// --------------------------------------------------------------------------
// Source code provided FOR REVIEW ONLY, as part of the submission entitled
// "Moving Level-of-Detail Surfaces".
//
// A proper version of this code will be released if the paper is accepted
// with the proper licence, documentation and bug fix.
// Currently, this material has to be considered confidential and shall not
// be used outside the review process.
//
// All right reserved. The Authors
// --------------------------------------------------------------------------

#pragma once

#include <unordered_map>
#include <array>
#include "octree.h"
#include "glm/glm.hpp"
#include <vector>

struct oriented_point
{
	oriented_point(){position = Eigen::Vector3f::Zero(); normal = Eigen::Vector3f::Zero();}
	oriented_point(glm::vec3 p, glm::vec3 n): position(Eigen::Vector3f(p.x, p.y, p.z)), normal(Eigen::Vector3f(n.x, n.y, n.z)){}
	oriented_point(Eigen::Vector3f p, Eigen::Vector3f n): position(p), normal(n){}

	Eigen::Vector3f position;
	Eigen::Vector3f normal;

	oriented_point operator +(oriented_point& p2){
		return oriented_point(position + p2.position, normal + p2.normal);}

	void operator +=(oriented_point& p2){
		position += p2.position;
		normal += p2.normal;}

	void operator /=(unsigned d){
		position /= float(d);
		normal /= float(d);}
};

struct vec3i_hash
{
	size_t operator()(const Eigen::Vector3i& id) const
	{
		return (size_t(id[0]) * size_t(73856093))
				^ (size_t(id[1]) * size_t(19349669))
				^ (size_t(id[2]) * size_t(39349693));
	}
};

template<typename Cell, typename Surface, typename Normal>
class dual_contour
{
public:
	dual_contour(const sorted_octree<Cell>& t, const Surface& s, const Normal& n) : tree(t), surf(s), normal_field(n) {}
	dual_contour(const sorted_octree<Cell>& t, const Surface& s, const Normal& n, std::unordered_map<Eigen::Vector3i, float, vec3i_hash>& c) : tree(t), surf(s), normal_field(n), corners(c) {}

	/// Extract the surface into a mesh
	void extract(std::vector<oriented_point>& vertices, std::vector<uint32_t>& indices)
	{
		// For each minimum edge, add 1 or 2 triangles (depending on whether it connects 3 or 4 nodes)
		std::vector<uint32_t> node_vertex(tree.size(), uint32_t(-1));
		auto get_vertex = [&](size_t n) -> uint32_t
		{
			uint32_t& vert = node_vertex[n];
			if(vert == uint32_t(-1))
			{
				vert = (uint32_t)vertices.size();
				vertices.push_back(compute_vertex(tree.id(n)));
			}
			return vert;
		};
		on_min_edges([&](const std::array<size_t, 4>& nodes, const std::array<unsigned int, 2>&, unsigned int, bool flip_tris)
		{
			auto add_triangle = [&](unsigned int i0, unsigned int i1, unsigned int i2)
			{
				if(i0 == i1 || i1 == i2 || i2 == i0) return;
				const uint32_t v[3] = { get_vertex(nodes[i0]), get_vertex(nodes[i1]), get_vertex(nodes[i2]) };
				if(v[0] == v[1] || v[1] == v[2] || v[2] == v[0]) return;
				indices.push_back(v[0]);
				indices.push_back(flip_tris ? v[2] : v[1]);
				indices.push_back(flip_tris ? v[1] : v[2]);
			};
			add_triangle(0, 1, 3);
			add_triangle(0, 3, 2);
		});
	}

protected:

	struct QEF
	{
		Eigen::Matrix3f A = Eigen::Matrix3f::Zero();
		Eigen::Vector3f b = Eigen::Vector3f::Zero();
		Eigen::Vector3f sum = Eigen::Vector3f::Zero();
		size_t count = 0;

		Eigen::Vector3f center() const { return sum / float(count); }
		void add(const Eigen::Vector3f& pos, const Eigen::Vector3f& normal)
		{
			A += normal * normal.transpose();
			b += pos.dot(normal) * normal;
			sum += pos;
			count++;
		}
		Eigen::Vector3f solve() const
		{
			const Eigen::Vector3f c = center();
			// Seems more precise than LDLT, which sometimes returns weird results for underconstrained systems
			return A.completeOrthogonalDecomposition().solve(b - A * c) + c;
		}
	};

protected:
	/// Compute the intersection of the surface and the given segment
	Eigen::Vector3f intersection(const Eigen::Vector3f& u, const Eigen::Vector3f& v, const float val_u, const float val_v) const
	{
		/// TODO: here we suppose linear interpolation...
		const float x = -val_u / (val_v - val_u);
		return x * v + (1 - x) * u;
	}

	/// Compute the position and orientation of the vertex of a node
	oriented_point compute_vertex(const octree_id& n)
	{
		// we only need to put the vertex at the center of the node, we will project the point on the surface with APSS later:
//		oriented_point vertex;
//		vertex.position = unit_octree::center(n);
//		return vertex;

		// QEF - based code :
		static const std::pair<uint8_t, uint8_t> node_edges[12] =
		{
			{ 0, 1 }, { 0, 2 }, { 1, 3 }, { 2, 3 },
			{ 4, 5 }, { 4, 6 }, { 5, 7 }, { 6, 7 },
			{ 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 }
		};
		const float w = 1.0f / (1 << tree.max_depth());
		const int delta = 1 << (tree.max_depth() - n.depth);
		oriented_point vertex;

		// Find edge intersections
		QEF qef;
		for(const std::pair<uint8_t, uint8_t>& e : node_edges)
		{
			const Eigen::Vector3i u = delta * n.corner(e.first);
			const Eigen::Vector3i v = delta * n.corner(e.second);
			const float val_u = get_corner(u);
			const float val_v = get_corner(v);
			if((val_u < 0) == (val_v < 0)) continue;

			// Intersection with implicit surface
			const Eigen::Vector3f p = intersection(u.cast<float>() * w, v.cast<float>() * w, val_u, val_v);
			const Eigen::Vector3f f = normal_field(p).normalized();
			vertex.normal += f;
			qef.add(p, f);
		}

		// Final position
		vertex.normal.normalize();
		vertex.position = qef.solve();
		// Special cases:
		//  1. Sometimes, due to the configuration of the normal field, the vertex might lie outside the node bounds,
		//     and this can create self-intersecting surfaces. In order to fix that, we discard normal data and just
		//     take the average of the intersection points.
		//  2. When a minimal edge is shared by 3 leaves only (i.e., there is a depth difference), then the minimal
		//     will only span a subpart of and edge of the biggest node, so it is possible that all corners of the
		//     biggest node are inside (or all outside) while the minimal edge definitely have a sign change, hence
		//     no points added to the qef. The simple fix here is to take the node center.
		if(!unit_octree::bounds(n).contains(vertex.position))
			vertex.position = qef.count > 0 ? qef.center() : unit_octree::center(n);
		if (std::isnan(vertex.position.x()) || std::isnan(vertex.position.y()) || std::isnan(vertex.position.z())) //NaN
			vertex.position = unit_octree::center(n);

		return vertex;
	}

	/// Return the value of the surface at the selected corner (and cache the result)
	float get_corner(const Eigen::Vector3i& c)
	{
		const float w = 1.0f / (1 << tree.max_depth());
		auto it = corners.find(c);
		if(it != corners.end()) return it->second;
		const float v = surf(c.cast<float>() * w);
		corners[c] = v;
		return v;
	}

	/// Call f(nodes, node_corners, selected_node, flip_triangles) on all minimum edges of the octree that have a sign change
	template<typename F>
	void on_min_edges(F f)
	{
		minimal_edges::traverse(tree, [&](size_t n1, size_t n2, size_t n3, size_t n4, direction dir)
		{
			ASSERT(n1 != size_t(-1) && n2 != size_t(-1) && n3 != size_t(-1) && n4 != size_t(-1), "Need valid nodes");
			ASSERT(n1 != n4 && n2 != n3, "Opposite nodes from a minimal edges cannot be the same");
			ASSERT((n1 == n2) + (n2 == n4) + (n4 == n3) + (n3 == n1) <= 1, "A minimal edge must be shared by at least 3 nodes");

			const std::array<size_t, 4> nodes = { { n1, n2, n3, n4 } };
			unsigned int selected = 0;
			int sel_depth = -1;
			for(unsigned i = 0; i < 4; i++)
				if(tree.id(nodes[i]).depth > sel_depth)
				{
					selected = i;
					sel_depth = tree.id(nodes[i]).depth;
				}

			static constexpr const unsigned int corner_indices[3][4][2] =
			{
				{ { 6, 7 }, { 4, 5 }, { 2, 3 }, { 0, 1 } },
				{ { 7, 5 }, { 6, 4 }, { 3, 1 }, { 2, 0 } },
				{ { 3, 7 }, { 2, 6 }, { 1, 5 }, { 0, 4 } }
			};
			const std::array<unsigned int, 2> c = { { corner_indices[(int)dir][selected][0], corner_indices[(int)dir][selected][1] } };
			const octree_id& n = tree.id(nodes[selected]);
			const int delta = 1 << (tree.max_depth() - n.depth);
			const Eigen::Vector3i vc0 = delta * n.corner(c[0]);
			const Eigen::Vector3i vc1 = delta * n.corner(c[1]);

			const float f0 = get_corner(vc0);
			const float f1 = get_corner(vc1);
			if(std::isnan(f0) || std::isnan(f1)) return;

			const bool m0 = f0 < 0;
			const bool m1 = f1 < 0;
			if(m0 != m1) f(nodes, c, selected, m0);
		});
	}

protected:
	const sorted_octree<Cell>& tree;
	const Surface& surf;
	const Normal& normal_field;
	std::unordered_map<Eigen::Vector3i, float, vec3i_hash> corners;
};

template<typename Cell, typename Surface, typename Normal>
dual_contour<Cell, Surface, Normal> dual_contouring(const sorted_octree<Cell>& tree, const Surface& surface, const Normal& normal_field)
{
	return { tree, surface, normal_field };
}
template<typename Cell, typename Surface, typename Normal>
dual_contour<Cell, Surface, Normal> dual_contouring(const sorted_octree<Cell>& tree, const Surface& surface, const Normal& normal_field, std::unordered_map<Eigen::Vector3i, float, vec3i_hash>& corners)
{
	return { tree, surface, normal_field, corners };
}
