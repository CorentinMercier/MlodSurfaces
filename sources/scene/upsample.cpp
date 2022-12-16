// --------------------------------------------------------------------------
// This file is part of the reference implementation for the paper
//    Moving Level-of-Detail Surfaces.
//    C. Mercier, T. Lescoat, P. Roussillon, T. Boubekeur, and J-M. Thiery
//    ACM Transaction On Graphics 2022
//    DOI: 10.1145/3528223.3530151
//
// All rights reserved. Use of this source code is governed by a
// MIT license that can be found in the LICENSE file.
// --------------------------------------------------------------------------

#include "upsample.h"

glm::vec3 to_glm(const Eigen::Vector3f& v) { return { v[0], v[1], v[2] }; }
Eigen::Vector3f to_eigen(const glm::vec3& v) { return { v[0], v[1], v[2] }; }


static constexpr const double PI = 3.14159265358979323846;

inline glm::vec3 fibonacci_point(size_t i, size_t n)
{
	const float alpha = PI * (3.0f - sqrt(5.0f));
	const float z = 1.0f - (2.0f * float(i) + 1.0f) / float(n);
	const float r = sqrt(1.0f - z * z);
	return glm::vec3(float(r * cos(alpha * i)), float(r * sin(alpha * i)), z);
}

void remove_duplicates(std::vector<octree_id>& vec)
{
	std::sort(vec.begin(), vec.end(), node_comparator);
	vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}

void remove_duplicates_reverse(std::vector<octree_id>& vec)
{
	std::sort(vec.rbegin(), vec.rend(), node_comparator);
	vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}

void regularize_up(std::vector<octree_id>& nodes)
{
	for(size_t t = 0; t < nodes.size(); t++)
		if(!nodes[t].is_root())
			nodes.push_back(nodes[t].parent());
	remove_duplicates(nodes);
}

void regularize_siblings(std::vector<octree_id>& nodes)
{
	const size_t s = nodes.size();
	auto add_siblings = [&](size_t t)
	{
		const octree_id& n = nodes[t];
		if(n.depth <= 0) return;
		const octree_id p = n.parent();
		for(unsigned int k = 0; k < 8; k++)
			nodes[s + 8 * t + k] = p.child(k);
	};

	nodes.resize(9 * s);
#pragma omp parallel for
	for(int t = 0; t < (int)s; t++)
		add_siblings(t);
	remove_duplicates(nodes);
}

void grade(std::vector<octree_id>& nodes)
{
	const std::vector<unsigned> corners = {0, 2, 6, 8, 18, 20, 24, 26};
	std::function<void(size_t)> add_neighbors;
	add_neighbors = [&nodes, &corners](size_t t)
	{
		if(nodes[t].depth <= 0) return;
		for (unsigned i : corners)
		{
			nodes.push_back(nodes[t].neighbor(i).parent());
		}
	};
	std::sort(nodes.rbegin(), nodes.rend(), node_comparator); //Reverse ordering, so that maxDepth values are first
	size_t i = 0;
	int currentDepth = nodes[0].depth;
	while(i<nodes.size())
	{
		if (currentDepth != nodes[i].depth)
		{
			remove_duplicates_reverse(nodes);
			currentDepth = nodes[i].depth;
		}
		add_neighbors(i);
		i++;
	}
	remove_duplicates(nodes);
}

void cut(std::vector<octree_id>& nodes, unsigned int max_depth)
{
	nodes.erase(std::remove_if(nodes.begin(), nodes.end(),
		[&](const octree_id& n) { return n.depth > int(max_depth); }),
		nodes.end());
}

void vector_difference_inplace(std::vector<octree_id>& vec, const std::vector<octree_id>& to_remove)
{
	auto should_remove = [&to_remove](const octree_id& n)
	{
		for(const octree_id& m : to_remove)
			if(n == m) return true;
		return false;
	};
	vec.erase(std::remove_if(vec.begin(), vec.end(), should_remove), vec.end());
}

void append_vector(std::vector<octree_id>& vec, const std::vector<octree_id>& to_add)
{
	vec.insert(vec.end(), to_add.begin(), to_add.end());
}

std::vector<octree_id> get_all_ancestors(const std::vector<octree_id>& nodes)
{
	std::vector<octree_id> res = nodes;
	for(size_t t = 0; t < res.size(); t++)
		if(res[t].depth > 0)
			res.push_back(res[t].parent());
	remove_duplicates(res);
	return res;
}

std::vector<octree_id> OctreeGenerator::generateFromPoints(unsigned int max_depth, const Pn& pts, const Eigen::AlignedBox3f &box) const
{
	if (pts.size() == 0) return {};
	aabb_octree octzone;
	if (box.isEmpty())
		octzone = {Eigen::AlignedBox3f(Eigen::Vector3f(-1.f, -1.f, -1.f), Eigen::Vector3f(1.f, 1.f, 1.f))};
	else
		octzone = {box};
	std::vector<octree_id> nodes(pts.size());
#pragma omp parallel for
	for (int i=0; i< (int)pts.size(); i++)
	{
		nodes[i] = octzone.enclosing(to_eigen(pts[i]), max_depth);
	}
	remove_duplicates(nodes);
	return nodes;
}

int getRelativeDepth(const std::vector<octree_id>& nodes, Eigen::Vector3i& max, Eigen::Vector3i& min)
{
	if (nodes.size() == 0) return 0;
	for(unsigned i=0; i<nodes.size(); i++)
	{
		octree_id ancestor = nodes[i].ancestor_at_depth(0);
		min = min.cwiseMin(ancestor.pos);
		max = max.cwiseMax(ancestor.pos);
	}

	const Eigen::Vector3i extent = max - min + Eigen::Vector3i(1, 1, 1);
	int delta_depth = 0;
	for(unsigned k = 0; k < 3; k++)
		delta_depth = std::max(delta_depth, (int)std::ceil(float(std::log2(extent[k]))));
	return delta_depth;
}

void onlySetDepthNodes(std::vector<octree_id>& nodes, unsigned int depth)
{
#pragma omp parallel for
	for (int i=0; i< (int)nodes.size(); i++)
	{
		if (unsigned(nodes[i].depth) > depth)
			nodes[i] = nodes[i].ancestor_at_depth(depth);
	}
	cut(nodes, depth);
	remove_duplicates(nodes);
}

bool intersectAABB_plane(const octree_id& node, const glm::vec3& pt, const glm::vec3& normal)
{
	const aabb_octree octzone =
	{
		Eigen::AlignedBox3f(Eigen::Vector3f( -1.0f, -1.0f, -1.0f), Eigen::Vector3f(1.0f, 1.0f, 1.0f))
	};
	glm::vec3 center = to_glm(octzone.center(node));
	glm::vec3 extent = glm::vec3(float(sqrt(3.)));
	float r = glm::dot(extent, glm::abs(normal));
	float s = glm::dot(normal, center) - glm::dot(normal, pt);
	return fabs(s) <= r;
}

void OctreeGenerator::expand_to_neighbors(std::vector<octree_id>& nodes, std::vector<octree_id>* total_tree) const
{
	std::cout << "Expanding to neighbors" << std::endl;
	const size_t originalSize = nodes.size();
	const aabb_octree octzone =
	{
		Eigen::AlignedBox3f(Eigen::Vector3f( -1.0f, -1.0f, -1.0f), Eigen::Vector3f(1.0f, 1.0f, 1.0f))
	};
	std::vector<glm::vec3> centers(nodes.size());
#pragma omp parallel for
	for (int i=0; i< (int)nodes.size(); i++)
	{
		Eigen::Vector3f c = octzone.center(nodes[i]);
		centers[i] = glm::vec3(c.x(), c.y(), c.z());
	}
	std::vector<glm::vec3> projectedNormals(nodes.size());
	apss->erasePointsFromGPU();
	apss->copyPointsToGPU(nodes.size(), centers.data());
	apss->project(nodes.size(), centers.data(), projectedNormals.data(), 1, kernel);

	for (unsigned i=0; i<originalSize; i++)
	{
		for(const octree_id& m : neighbors(nodes[i]))
		{
			if (intersectAABB_plane(m, centers[i], projectedNormals[i]))
				nodes.push_back(m);
		}
	}

	remove_duplicates(nodes);
	const size_t addedNodes = nodes.size() - originalSize;
	std::cout << addedNodes << " nodes were added" << std::endl;
	if (total_tree && addedNodes > 0)
	{
		remove_duplicates(*total_tree);
		append_vector(*total_tree, nodes);
		remove_duplicates(*total_tree);
	}
}

std::vector<octree_id> OctreeGenerator::generate(unsigned int max_depth, std::vector<octree_id>* total_tree) const
{
	std::vector<octree_id> depth0 = { octree_id() };
	std::vector<octree_id> old_depth0, nodes;
	unsigned int max_relative_depth = max_depth;
	Eigen::Vector3i min(INT_MAX, INT_MAX, INT_MAX);
	Eigen::Vector3i max = -min;
	while(!depth0.empty())
	{
		if(total_tree) append_vector(*total_tree, depth0);
		std::vector<octree_id> current = compute_leaves(depth0, max_relative_depth, total_tree);
		max_relative_depth = std::max(0, int(max_depth) - getRelativeDepth(current, max, min));//To comment for no adaptive max depth estimation
		append_vector(nodes, current);
		append_vector(old_depth0, depth0);
		depth0 = depth0_ancestors(current);
		expand(depth0);
		vector_difference_inplace(depth0, old_depth0);
	}
	remove_duplicates(nodes);
	onlySetDepthNodes(nodes, max_relative_depth);//To comment for no adaptive max depth estimation
	if(total_tree)
		onlySetDepthNodes(*total_tree, max_relative_depth);//To comment for no adaptive max depth estimation
	return nodes;
}

void uniformizeFromPoints(std::vector<octree_id> nodes, sorted_octree<char>& soct, glm::vec3& new_diag, glm::vec3& new_origin)
{
	grade(nodes);
	regularize_up(nodes);

	regularize_siblings(nodes);
	soct = sorted_octree<char>(nodes);
}

void uniformize(unsigned int max_depth, std::vector<octree_id> nodes,
	sorted_octree<char>& soct, glm::vec3& new_diag, glm::vec3& new_origin, bool doCut)
{
	append_vector(nodes, depth0_ancestors(nodes));
	retreiveCorrectRoot(nodes, new_diag, new_origin);
	grade(nodes);
	regularize_up(nodes);
	if (doCut) cut(nodes, max_depth);
	regularize_siblings(nodes);
	soct = sorted_octree<char>(nodes);
}

std::vector<octree_id> OctreeGenerator::compute_leaves(std::vector<octree_id> current_nodes, unsigned int max_depth, std::vector<octree_id>* total_tree/*, unsigned a, const Eigen::Matrix4f& transform*/) const
{
	std::vector<octree_id> toSubdivide;
	for(unsigned int depth = 0; depth < max_depth + 1; depth++)
	{
		find_stable_nodesGPU(current_nodes, toSubdivide);
		current_nodes.clear();
		add_children(toSubdivide, current_nodes);
		toSubdivide.clear();
		expand(current_nodes, depth + 1);//Explore neighbors
		if(total_tree) append_vector(*total_tree, current_nodes);
	}
	return current_nodes;
}

std::vector<octree_id> OctreeGenerator::project(const octree_id& octId) const
{
	const aabb_octree octzone =
	{
		Eigen::AlignedBox3f(Eigen::Vector3f( -1.0f, -1.0f, -1.0f), Eigen::Vector3f(1.0f, 1.0f, 1.0f))
	};
	const Eigen::Vector3f c = octzone.center(octId);
	glm::vec3 projectedPoint, projectedNormal;
	glm::vec3 center(c.x(), c.y(), c.z());
	apss->projectCPU(center, projectedPoint, projectedNormal, 1, kernel, std::numeric_limits<float>::max());
	for(unsigned k = 0; k < 3; k++)
		if(std::isinf(projectedPoint[k]) || std::isnan(projectedPoint[k]))
			return {};
	octree_id newId = octzone.enclosing(Eigen::Vector3f(projectedPoint.x, projectedPoint.y, projectedPoint.z), octId.depth);

	return { newId };// nodesToSubdivide;
}

void OctreeGenerator::find_stable_nodes(const std::vector<octree_id>& nodes, std::vector<octree_id>& out) const
{
	for(const octree_id& n : nodes)
	{
		const std::vector<octree_id> projected = project(n);
		out.insert(out.end(), projected.begin(), projected.end());
	}
	remove_duplicates(out);
}

void OctreeGenerator::find_stable_nodesGPU(const std::vector<octree_id>& nodes, std::vector<octree_id>& out) const
{
	const aabb_octree octzone =
	{
		Eigen::AlignedBox3f(Eigen::Vector3f( -1.0f, -1.0f, -1.0f), Eigen::Vector3f(1.0f, 1.0f, 1.0f))
	};
	std::vector<glm::vec3> centers(nodes.size());
#pragma omp parallel for
	for (int i=0; i< (int)nodes.size(); i++)
	{
		Eigen::Vector3f c = octzone.center(nodes[i]);
		centers[i] = glm::vec3(c.x(), c.y(), c.z());
	}
	std::vector<glm::vec3> projectedNormals(nodes.size());
	apss->erasePointsFromGPU();
	apss->copyPointsToGPU(nodes.size(), centers.data());
	apss->project(nodes.size(), centers.data(), projectedNormals.data(), 1, kernel, std::numeric_limits<float>::max());

	out.resize(nodes.size());
#pragma omp parallel for
	for (int i=0; i< (int)out.size(); i++)
	{
		if (glm::length2(projectedNormals[i]) < 0.5f)
			out[i] = octree_id();//Point was not projected, so we do not wish to add a bad octree_id. Root is always included.
		else
			out[i] = octzone.enclosing(Eigen::Vector3f(centers[i].x, centers[i].y, centers[i].z), nodes[i].depth);
	}

	remove_duplicates(out);
}

void OctreeGenerator::add_children(const std::vector<octree_id>& nodes, std::vector<octree_id>& out) const
{
	out.resize(nodes.size() * 8);
#pragma omp parallel for
	for (int i = 0; i < (int)nodes.size(); i++)
	{
		const octree_id& n = nodes[i];
		for(unsigned int c = 0; c < 8; c++)
			out[8 * i + c] = n.child(c);
	}
}

std::vector<octree_id> depth0_ancestors(std::vector<octree_id> nodes)
{
#pragma omp parallel for
	for(int i = 0; i < (int)nodes.size(); i++)
		nodes[i] = nodes[i].ancestor_at_depth(0);
	remove_duplicates(nodes);
	return nodes;
}

Eigen::Vector4i union_nodes(std::vector<octree_id>& nodes)
{
	Eigen::Vector3i min(INT_MAX, INT_MAX, INT_MAX);
	Eigen::Vector3i max = -min;
	bool has_depth0 = false;
	for(const octree_id& n : nodes)
	{
		if(n.depth != 0) continue;
		min = min.cwiseMin(n.pos);
		max = max.cwiseMax(n.pos);
		has_depth0 = true;
	}
	if(!has_depth0) return { -1, 0, 0, 0 };

	const Eigen::Vector3i extent = max - min + Eigen::Vector3i(1, 1, 1);
	int delta_depth = 0;
	for(unsigned k = 0; k < 3; k++)
		delta_depth = std::max(delta_depth, (int)std::ceil(float(std::log2(extent[k])))); /// TODO: use integer logarithm
	if(delta_depth == 0) return { -1, 0, 0, 0 };

#pragma omp parallel for
	for(int i=0; i< (int)nodes.size(); i++)
	{
		octree_id& n = nodes[i];
		n.pos -= min * (1 << n.depth);
		n.depth += delta_depth;
	}

	return { delta_depth, min[0], min[1], min[2] };
}

Eigen::Vector4i retreiveCorrectRoot(std::vector<octree_id>& nodes, glm::vec3& new_diag, glm::vec3& new_origin)
{
	const Eigen::Vector4i delta = union_nodes(nodes);
	const octree_id old_min_depth0(0, delta[1], delta[2], delta[3]);
	new_origin = to_glm(unit_octree::bounds(old_min_depth0).min());
	new_diag = glm::vec3(1 << delta[0]);
	return delta;
}

void expand(std::vector<octree_id>& nodes, int depth)
{
	const size_t n = nodes.size();
	for(size_t t = 0; t < n; t++)
	{
		if(depth >= 0 && nodes[t].depth != depth) continue;
		for(const octree_id& m : neighbors(nodes[t]))
			nodes.push_back(m);
	}
	remove_duplicates(nodes);
}
