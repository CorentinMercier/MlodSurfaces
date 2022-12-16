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

#pragma once

#include <unordered_map>
#include "base.h"
#include <vector>

struct octree_id
{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	octree_id(int d, int x, int y = 0, int z = 0) : depth(d), pos(x, y, z) {}
	octree_id(int d = 0, const Eigen::Vector3i& p = { 0, 0, 0 }) : depth(d), pos(p) {/* ASSERT(is_valid(), "Invalid octree id"); */}
	bool ascend();
	bool is_root() const { return depth == 0; }
	bool is_valid() const;
	octree_id parent() const;
	octree_id ancestor_at_depth(int d) const;
	octree_id child(unsigned int c) const;
	Eigen::Vector3i corner(unsigned int c) const;
	octree_id neighbor(unsigned int c) const;
	bool is_ancestor_of(octree_id id) const;
	bool is_interior() const;
	bool operator==(const octree_id& id) const;
	bool operator!=(const octree_id& id) const { return !(*this == id); }
	Eigen::AlignedBox3i ring(unsigned int size) const;

	static int width(int depth) { return 1 << depth; }
	static octree_id invalid();
	static Eigen::AlignedBox3i space(int depth);

	int depth;
	Eigen::Vector3i pos;
};

struct octree_id_neighbors
{
	struct iterator
	{
		const octree_id c;
		octree_id n;
		unsigned char t;

		void validate()
		{
			if (t==13) t++;
			if(t < 27) n = c.neighbor(t);
		}
		iterator(const octree_id& c, unsigned char t) : c(c), t(t) { validate(); }
		iterator& operator++() { t++; validate(); return *this; }
		bool operator!=(const iterator& i) const { return i.t != t; }
		octree_id operator*() const { return n; }
	};

	octree_id c;
	octree_id_neighbors(const octree_id& c) : c(c) {}
	iterator begin() { return iterator(c, 0); }
	iterator end() { return iterator(c, 27); }
};
inline octree_id_neighbors neighbors(const octree_id& c) { return c; }

/// An octree, stored in a flat hashed container
template<typename T>
class octree
{
public:
	bool add(const octree_id& id, T data);
	void remove(const octree_id& id);

	bool contains(const octree_id& id) const;
	T& data(const octree_id& id);
	const T& data(const octree_id& id) const;
	T& operator[](const octree_id& id);
	T* find(const octree_id& id);
	const T* find(const octree_id& id) const;

	size_t size() const { return nodes.size(); }
	size_t size(int depth) const;
	int max_depth() const;

	template<typename F>
	void for_each(F f) const;

	template<typename F>
	void for_each(F f)
	{
		for(typename node_storage::iterator it = nodes.begin(); it != nodes.end(); ++it)
			f(it->first, it->second);
	}

	std::vector<octree_id> all_nodes() const;
	std::vector<octree_id> all_nodes(int depth) const;
	std::vector<octree_id> all_min_depth_nodes(int depth) const;
	std::vector<octree_id> all_boundary_nodes() const;

	/// Ensure that for each node its hierarchy is present, and its sibling too
	void regularize(const T& def = T());

	/// Add all neighbors of all nodes at \p depth
	void expand(int depth, const T& def = T());

private:
	struct id_hash
	{
		size_t operator()(const octree_id& id) const
		{
			return (size_t(id.depth) * size_t(39349693))
				^ (size_t(id.pos[0]) * size_t(73856093))
				^ (size_t(id.pos[1]) * size_t(19349669))
				^ (size_t(id.pos[2]) * size_t(83492791));
		}
	};

private:
	using node_storage = std::unordered_map<octree_id, T, id_hash>;
	node_storage nodes;
};

struct unit_octree
{
	static Eigen::AlignedBox3f bounds();
	static Eigen::AlignedBox3f bounds(const octree_id& id);
	static Eigen::Vector3f center(const octree_id& id);
	static octree_id enclosing(const Eigen::Vector3f& pos, int depth);
};

struct aabb_octree
{
	Eigen::AlignedBox3f box;

	Eigen::AlignedBox3f bounds() const;
	Eigen::AlignedBox3f bounds(const octree_id& id) const;
	Eigen::Vector3f center(const octree_id& id) const;
	octree_id enclosing(const Eigen::Vector3f& pos, int depth) const;
};

inline bool node_comparator(const octree_id& a, const octree_id& b)
{
	if(a.depth != b.depth) return a.depth < b.depth;
	// Keep the order of children, so sort first on Z, then Y, then X
	if(a.pos[2] != b.pos[2]) return a.pos[2] < b.pos[2];
	if(a.pos[1] != b.pos[1]) return a.pos[1] < b.pos[1];
	return a.pos[0] < b.pos[0];
}

struct sorted_nodes
{
	std::vector<octree_id> nodes;
	std::vector<size_t> depth_layers;

	template<typename T>
	sorted_nodes(const octree<T>& tree)
	{
		nodes = tree.all_nodes();
		std::sort(nodes.begin(), nodes.end(), node_comparator);
		depth_layers.push_back(0);
		for(size_t t = 0; t < nodes.size(); t++)
			if(t == 0 || nodes[t].depth != nodes[t - 1].depth)
				while(depth_layers.size() <= unsigned(nodes[t].depth))
					depth_layers.push_back(t);
	}

	sorted_nodes(std::vector<octree_id> nodes) : nodes(std::move(nodes))
	{
		std::sort(nodes.begin(), nodes.end(), node_comparator);
		depth_layers.push_back(0);
		for(size_t t = 0; t < nodes.size(); t++)
			if(t == 0 || nodes[t].depth != nodes[t - 1].depth)
				while(depth_layers.size() <= unsigned(nodes[t].depth))
					depth_layers.push_back(t);
	}

	std::vector<octree_id>::const_iterator begin() const { return nodes.begin(); }
	std::vector<octree_id>::const_iterator end() const { return nodes.end(); }
	std::vector<octree_id>::iterator begin() { return nodes.begin(); }
	std::vector<octree_id>::iterator end() { return nodes.end(); }
	size_t begin(size_t depth) const { return depth_layers[depth]; }
	size_t end(size_t depth) const { return depth == depth_layers.size() - 1 ? nodes.size() : depth_layers[depth + 1]; }
	const octree_id& operator[](size_t t) const { return nodes[t]; }
	size_t size() const { return nodes.size(); }
	size_t relative_index(const octree_id& n) const
	{
		const std::vector<octree_id>::const_iterator s = nodes.begin() + begin(n.depth);
		const std::vector<octree_id>::const_iterator e = nodes.begin() + end(n.depth);
		const std::vector<octree_id>::const_iterator it = std::lower_bound(s, e, n, node_comparator);
		return it != e && *it == n ? std::distance(s, it) : -1;
	}
	size_t index(const octree_id& n) const
	{
		const size_t i = relative_index(n);
		return i != size_t(-1) ? begin(n.depth) + i : i;
	}
};

template<typename T>
class sorted_octree
{
public:
	struct range_node_ids
	{
		struct iterator
		{
			const sorted_octree<T>& tree;
			size_t t;
			bool operator!=(const iterator& it) const { return t != it.t; }
			iterator& operator++() { t++; return *this; }
			const octree_id& operator*() const { return tree.id(t); }
		};

		const sorted_octree<T>& tree;
		size_t start, stop;
		iterator begin() const { return { tree, start }; }
		iterator end() const { return { tree, stop }; }
	};

public:
	sorted_octree(){}
	explicit sorted_octree(const octree<T>& tree) : sorted_octree(tree.all_nodes())
	{
	}

	explicit sorted_octree(const std::vector<octree_id>& all_nodes, T val = T())
	{
		// Nodes
		nodes.reserve(all_nodes.size());
		for(const octree_id& n : all_nodes)
			nodes.emplace_back(n, val);
		std::sort(nodes.begin(), nodes.end());

		// Depths
		depths.push_back(0);
		for(size_t t = 0; t < nodes.size(); t++)
			if(t == 0 || nodes[t].id.depth != nodes[t - 1].id.depth)
				while(depths.size() <= unsigned(nodes[t].id.depth))
					depths.push_back(t);

		// Relations : children
		relations.resize(size());
		{
			struct parent_rel { uint32_t parent, child; };
			std::vector<parent_rel> rels;
			for(size_t t = 0; t < size(); t++)
			{
				if(id(t).is_root()) continue;
				const size_t i = index(id(t).parent());
				ASSERT(i != size_t(-1), "All nodes should have an explicit parent (except root)");
				rels.push_back(parent_rel { uint32_t(i), uint32_t(t) });
			}
			std::sort(rels.begin(), rels.end(), [](const parent_rel& a, const parent_rel& b)
			{
				return a.parent != b.parent ? a.parent < b.parent : a.child < b.child;
			});

			abs_children.reserve(rels.size());
			uint32_t last_parent = uint32_t(-1);
			for(const parent_rel& r : rels)
			{
				relations[r.parent].num_children++;
				if(last_parent != r.parent)
					relations[r.parent].children = (uint32_t)abs_children.size();
				last_parent = r.parent;
				abs_children.push_back(r.child);
			}
		}

		// Relations : one rings
		{
			struct ring1_rel { uint32_t center, n; };
			std::vector<ring1_rel> rels;
			for(size_t t = 0; t < size(); t++)
			{
				const octree_id& n = id(t);
				for(int dx = -1; dx <= 1; dx++)
				for(int dy = -1; dy <= 1; dy++)
				for(int dz = -1; dz <= 1; dz++)
				{
					const size_t i = dx == 0 && dy == 0 && dz == 0
						? t
						: index(octree_id(n.depth, n.pos[0] + dx, n.pos[1] + dy, n.pos[2] + dz));
					if(i == size_t(-1)) continue;
					rels.push_back(ring1_rel { uint32_t(t), uint32_t(i) });
				}
			}
			std::sort(rels.begin(), rels.end(), [](const ring1_rel& a, const ring1_rel& b)
			{
				return a.center != b.center ? a.center < b.center : a.n < b.n;
			});

			abs_ring1.reserve(rels.size());
			uint32_t last_node = uint32_t(-1);
			for(const ring1_rel& r : rels)
			{
				relations[r.center].num_ring1++;
				if(last_node != r.center)
					relations[r.center].ring1 = (uint32_t)abs_ring1.size();
				last_node = r.center;
				abs_ring1.push_back(r.n);
			}
		}
	}

	T& data(const octree_id& id) { return nodes[index(id)].data; }
	const T& data(const octree_id& id) const { return nodes[index(id)].data; }
	T& operator[](size_t t) { return nodes[t].data; }
	const T& operator[](size_t t) const { return nodes[t].data; }
	T& data(size_t t) { return nodes[t].data; }
	const T& data(size_t t) const { return nodes[t].data; }
	const octree_id& id(size_t t) const { return nodes[t].id; }

	pointer_range<const uint32_t> children(size_t t) const
	{
		const node_relations& r = relations[t];
		const uint32_t* begin = abs_children.data() + r.children;
		return { begin, begin + r.num_children };
	}

	pointer_range<const uint32_t> ring1(size_t t) const
	{
		const node_relations& r = relations[t];
		const uint32_t* begin = abs_ring1.data() + r.ring1;
		return { begin, begin + r.num_ring1 };
	}

	size_t begin(int depth) const { return depths[depth]; }
	size_t end(int depth) const { return depth == int(depths.size()) - 1 ? nodes.size() : depths[depth + 1]; }

	size_t size() const { return nodes.size(); }
	size_t size(int depth) const { return end(depth) - begin(depth); }
	int max_depth() const { return int(depths.size()) - 1; }

	unsigned num_children(size_t n) const { return relations[n].num_children; }
	bool is_leaf(size_t n) const { return relations[n].num_children == 0; }

	size_t relative_index(const octree_id& n) const
	{
		const typename std::vector<node>::const_iterator s = nodes.begin() + begin(n.depth);
		const typename std::vector<node>::const_iterator e = nodes.begin() + end(n.depth);
		const typename std::vector<node>::const_iterator it = std::lower_bound(s, e, n);
		return it != e && it->id == n ? std::distance(s, it) : -1;
	}
	size_t index(const octree_id& n) const
	{
		const size_t i = relative_index(n);
		return i != size_t(-1) ? begin(n.depth) + i : i;
	}

	template<typename F>
	void for_each(F f) const
	{
		for(size_t t = 0; t < nodes.size(); t++)
			f(nodes[t].id, nodes[t].data);
	}

	range_node_ids all_nodes() const { return { *this, 0, nodes.size() }; }
	range_node_ids all_nodes(int depth) const { return { *this, begin(depth), end(depth) }; }

private:
	struct node
	{
		T data;
		octree_id id;
		node(const octree_id& n, const T& d = T()) : data(d), id(n) {}
		bool operator<(const node& n) const { return node_comparator(id, n.id); }
	};
	struct node_relations
	{
		uint32_t children = 0;
		uint32_t ring1 = 0;
		uint8_t num_children = 0;
		uint8_t num_ring1 = 0;
		node_relations() {}
	};

private:
	std::vector<node> nodes;
	std::vector<size_t> depths;
	std::vector<node_relations> relations;
	std::vector<uint32_t> abs_children;
	std::vector<uint32_t> abs_ring1;
};

/** A *minimal edge* is one that does not contain any other edge in a neighbouring octree cell.
	In the case where every node has either 8 children (intermediate node) or 0 children (leaf
	node), minimal edges are only in leaf nodes: any edge of an intermediate node can be
	subdivided along the node's children, and hence is not minimal.
	The algorithm here is from "Dual Contouring of Hermite Data", and enumerate them in O(n)
	with n the number of nodes in the octree. A graphical illustration of it can be found in
	"Adaptive surface extraction from anisotropic volumetric data: contouring on generalized
	octrees", from Lobello, Denis and Dupont, page 8. */
struct minimal_edges
{
	/* Children:
		.     6        7
		.      +------+       Z
		.     /|     /|       | / Y
		.  4 / |   5/ |       |/
		.   +------+  |       +---- X
		.   |  +---|--+
		.   | /2   | / 3
		.   |/     |/
		.   +------+
		.  0        1
	*/
	template<typename T, typename F>
	static void traverse(const sorted_octree<T>& tree, F f)
	{
		const size_t root = tree.index(octree_id());
		if(root != size_t(-1)) traverse_node(root, tree, f);
	}

private:
	template<typename T, typename F>
	static void traverse_node(size_t n, const sorted_octree<T>& tree, F fn)
	{
		if(tree.is_leaf(n)) return;
		const uint32_t* children = tree.children(n).begin();
		ASSERT(tree.children(n).end() - children == 8, "Node must have exactly 8 children");

		// Inner nodes
		for(unsigned int t = 0; t < 8; t++)
			traverse_node(children[t], tree, fn);

		// Inner faces
		traverse_face(children[0], children[1], direction::X, tree, fn); // X
		traverse_face(children[2], children[3], direction::X, tree, fn);
		traverse_face(children[4], children[5], direction::X, tree, fn);
		traverse_face(children[6], children[7], direction::X, tree, fn);
		traverse_face(children[0], children[2], direction::Y, tree, fn); // Y
		traverse_face(children[1], children[3], direction::Y, tree, fn);
		traverse_face(children[4], children[6], direction::Y, tree, fn);
		traverse_face(children[5], children[7], direction::Y, tree, fn);
		traverse_face(children[0], children[4], direction::Z, tree, fn); // Z
		traverse_face(children[1], children[5], direction::Z, tree, fn);
		traverse_face(children[2], children[6], direction::Z, tree, fn);
		traverse_face(children[3], children[7], direction::Z, tree, fn);

		// Inner edges
		traverse_edge(children[0], children[2], children[4], children[6], direction::X, tree, fn); // X
		traverse_edge(children[1], children[3], children[5], children[7], direction::X, tree, fn);
		traverse_edge(children[0], children[1], children[4], children[5], direction::Y, tree, fn); // Y
		traverse_edge(children[2], children[3], children[6], children[7], direction::Y, tree, fn);
		traverse_edge(children[0], children[1], children[2], children[3], direction::Z, tree, fn); // Z
		traverse_edge(children[4], children[5], children[6], children[7], direction::Z, tree, fn);
	}

	template<typename T, typename F>
	static void traverse_face(size_t n1, size_t n2, direction dir, const sorted_octree<T>& tree, F fn)
	{
		ASSERT(tree.is_leaf(n1) || tree.num_children(n1) == 8, "Node must have exactly 8 children");
		ASSERT(tree.is_leaf(n2) || tree.num_children(n2) == 8, "Node must have exactly 8 children");
		const uint32_t* children1 = tree.is_leaf(n1) ? nullptr : tree.children(n1).begin();
		const uint32_t* children2 = tree.is_leaf(n2) ? nullptr : tree.children(n2).begin();
		if(!children1 && !children2) return;

		// Inner faces
		static constexpr const unsigned int subfaces[3][2][4] =
		{
			{ { 1, 3, 5, 7 }, { 0, 2, 4, 6 } }, // Depending on the direction (first index) and whether it is
			{ { 2, 3, 6, 7 }, { 0, 1, 4, 5 } }, // about n1 or n2 (second index), these are the index of the
			{ { 4, 5, 6, 7 }, { 0, 1, 2, 3 } }  // children to which we propagate the recursion.
		};
		for(unsigned int f = 0; f < 4; f++)
			traverse_face(children1 ? children1[subfaces[(int)dir][0][f]] : n1,
				children2 ? children2[subfaces[(int)dir][1][f]] : n2,
				dir, tree, fn);

		// Inner edges
		static constexpr const unsigned int subedges[3][4][4] =
		{
			{ { 1, 5, 0, 4 }, { 3, 7, 2, 6 }, { 1, 3, 0, 2 }, { 5, 7, 4, 6 } }, // Per direction, per inner edge
			{ { 3, 7, 1, 5 }, { 2, 6, 0, 4 }, { 2, 3, 0, 1 }, { 6, 7, 4, 5 } }, // children indices
			{ { 4, 5, 0, 1 }, { 6, 7, 2, 3 }, { 4, 6, 0, 2 }, { 5, 7, 1, 3 } }
		};
		static constexpr const direction subedgedir[3][4] =
		{
			{ direction::Y, direction::Y, direction::Z, direction::Z },
			{ direction::X, direction::X, direction::Z, direction::Z },
			{ direction::Y, direction::Y, direction::X, direction::X }
		};
		for(unsigned int e = 0; e < 4; e++)
		{
			const size_t m1 = children1 ? children1[subedges[(int)dir][e][0]] : n1;
			const size_t m2 = children1 ? children1[subedges[(int)dir][e][1]] : n1;
			const size_t m3 = children2 ? children2[subedges[(int)dir][e][2]] : n2;
			const size_t m4 = children2 ? children2[subedges[(int)dir][e][3]] : n2;
			// To keep the ordering of the nodes it is needed to flip them in some cases
			const bool flip = dir == direction::X || (dir == direction::Y && e < 2);
			if(!flip) traverse_edge(m1, m2, m3, m4, subedgedir[(int)dir][e], tree, fn);
			else      traverse_edge(m1, m3, m2, m4, subedgedir[(int)dir][e], tree, fn);
		}
	}

	template<typename T, typename F>
	static void traverse_edge(size_t n1, size_t n2, size_t n3, size_t n4, direction dir, const sorted_octree<T>& tree, F fn)
	{
		ASSERT(tree.is_leaf(n1) || tree.num_children(n1) == 8, "Node must have exactly 8 children");
		ASSERT(tree.is_leaf(n2) || tree.num_children(n2) == 8, "Node must have exactly 8 children");
		ASSERT(tree.is_leaf(n3) || tree.num_children(n3) == 8, "Node must have exactly 8 children");
		ASSERT(tree.is_leaf(n4) || tree.num_children(n4) == 8, "Node must have exactly 8 children");
		const uint32_t* children1 = tree.is_leaf(n1) ? nullptr : tree.children(n1).begin();
		const uint32_t* children2 = tree.is_leaf(n2) ? nullptr : tree.children(n2).begin();
		const uint32_t* children3 = tree.is_leaf(n3) ? nullptr : tree.children(n3).begin();
		const uint32_t* children4 = tree.is_leaf(n4) ? nullptr : tree.children(n4).begin();
		if(!children1 && !children2 && !children3 && !children4)
		{
			// This is a minimal edge
			fn(n1, n2, n3, n4, dir);
			return;
		}

		// Inner edges
		static constexpr const unsigned int subedges[3][2][4] =
		{
			{ { 6, 4, 2, 0 }, { 7, 5, 3, 1 } }, // Per direction, per edge, index
			{ { 5, 4, 1, 0 }, { 7, 6, 3, 2 } }, // of the children nodes sharing
			{ { 3, 2, 1, 0 }, { 7, 6, 5, 4 } }  // an inner edge.
		};
		for(unsigned e = 0; e < 2; e++)
			traverse_edge(
				children1 ? children1[subedges[(int)dir][e][0]] : n1,
				children2 ? children2[subedges[(int)dir][e][1]] : n2,
				children3 ? children3[subedges[(int)dir][e][2]] : n3,
				children4 ? children4[subedges[(int)dir][e][3]] : n4,
				dir, tree, fn);
	}
};

// Implementation
// ------------------------------------------------------------------------------------------------
inline bool octree_id::ascend()
{
	if(is_root()) return false;
	*this = parent();
	return true;
}

inline octree_id octree_id::parent() const
{
	ASSERT(depth > 0, "Should not get parent of root");
	octree_id id = *this;
	if(!id.is_root())
	{
		id.depth--;
		const Eigen::Vector3i d = ((id.pos.array() < 0) && ((id.pos.array() - (2 * (id.pos.array()/2))) != 0)).cast<int>();
		id.pos /= 2;
		id.pos -= d;
	}
	return id;
}

inline octree_id octree_id::ancestor_at_depth(int d) const
{
	octree_id id = *this;
	ASSERT (id.depth >= d, "Parent must be at a lower depth");
	while(id.depth > d)
	{
		id.depth--;
		const Eigen::Vector3i d2 = ((id.pos.array() < 0) && ((id.pos.array() - (2 * (id.pos.array()/2))) != 0)).cast<int>();
		id.pos /= 2;
		id.pos -= d2;
	}
	return id;
}

inline octree_id octree_id::child(unsigned int n) const
{
	ASSERT(n < 8, "Invalid octree node child id");

	return octree_id(depth + 1,
		2 * pos[0] + int((n & 1) >> 0),
		2 * pos[1] + int((n & 2) >> 1),
		2 * pos[2] + int((n & 4) >> 2));
}

inline Eigen::Vector3i octree_id::corner(unsigned int c) const
{
	ASSERT(c < 8, "Invalid octree node corner id");
	return pos + Eigen::Vector3i(int((c & 1) >> 0), int((c & 2) >> 1), int((c & 4) >> 2));
}

inline octree_id octree_id::neighbor(unsigned int n) const
{
	ASSERT(n < 27, "Invalid octree node neighbor id");
	//ASSERT(is_valid(), "Invalid octree id");

	const int z = n / 9;
	const int y = (n - 9 * z) / 3;
	const int x = n - 9 * z - 3 * y;
	octree_id c = *this;
	c.pos += Eigen::Vector3i(x - 1, y - 1, z - 1);
	return c;
}

inline bool octree_id::is_valid() const
{
	const int w = width(depth);
	return depth >= 0
		&& pos[0] >= 0 && pos[0] < w
		&& pos[1] >= 0 && pos[1] < w
		&& pos[2] >= 0 && pos[2] < w;
}

inline bool octree_id::is_ancestor_of(octree_id id) const
{
	if(id.depth <= depth) return false;
	while(id.depth > depth)
		id = id.parent();
	return id == *this;
}

inline bool octree_id::is_interior() const
{
	const int w = width(depth);
	return pos[0] > 0 && pos[0] < w - 1
		&& pos[1] > 0 && pos[1] < w - 1
		&& pos[2] > 0 && pos[2] < w - 1;
}

inline bool octree_id::operator==(const octree_id& id) const
{
	return depth == id.depth && pos == id.pos;
}

inline Eigen::AlignedBox3i octree_id::ring(unsigned int s) const
{
	return Eigen::AlignedBox3i(pos.array() - s, pos.array() + s);
}

inline octree_id octree_id::invalid()
{
	octree_id c;
	c.depth = c.pos[0] = c.pos[1] = c.pos[2] = -1;
	return c;
}

inline Eigen::AlignedBox3i octree_id::space(int depth)
{
	return Eigen::AlignedBox3i(Eigen::Vector3i::Zero(), Eigen::Vector3i::Constant(width(depth) - 1));
}

template<typename T>
inline bool octree<T>::add(const octree_id& id, T data)
{
	return nodes.insert({ id, std::move(data) }).second;
}

template<typename T>
inline void octree<T>::remove(const octree_id& id)
{
	for(typename node_storage::iterator it = nodes.begin(); it != nodes.end();)
		if(it->first == id || id.is_ancestor_of(it->first))
			it = nodes.erase(it);
		else
			++it;
}

template<typename T>
inline bool octree<T>::contains(const octree_id& id) const
{
	return nodes.find(id) != nodes.end();
}

template<typename T>
inline T& octree<T>::data(const octree_id& id)
{
	auto it = nodes.find(id);
	ASSERT(it != nodes.end(), "Accessing non-existing octree node");
	return it->second;
}

template<typename T>
inline const T& octree<T>::data(const octree_id& id) const
{
	auto it = nodes.find(id);
	ASSERT(it != nodes.end(), "Accessing non-existing octree node");
	return it->second;
}

template<typename T>
inline T* octree<T>::find(const octree_id& id)
{
	auto it = nodes.find(id);
	return it != nodes.end() ? &it->second : nullptr;
}

template<typename T>
inline const T* octree<T>::find(const octree_id& id) const
{
	auto it = nodes.find(id);
	return it != nodes.end() ? &it->second : nullptr;
}

template<typename T>
inline T& octree<T>::operator[](const octree_id& id)
{
	return nodes[id];
}

template<typename T>
inline size_t octree<T>::size(int depth) const
{
	return std::count_if(nodes.begin(), nodes.end(), [&](const std::pair<const octree_id, T>& n) { return n.first.depth == depth; });
}

template<typename T>
inline int octree<T>::max_depth() const
{
	/// TODO: faster implementation
	int d = -999999;
	for(typename node_storage::const_iterator it = nodes.begin(); it != nodes.end(); ++it)
		if(it->first.depth > d)
			d = it->first.depth;
	return d;
}

template<typename T>
template<typename F>
inline void octree<T>::for_each(F f) const
{
	for(typename node_storage::const_iterator it = nodes.begin(); it != nodes.end(); ++it)
		f(it->first, it->second);
}

template<typename T>
inline std::vector<octree_id> octree<T>::all_nodes() const
{
	std::vector<octree_id> cells;
	cells.reserve(nodes.size());
	for(typename node_storage::const_iterator it = nodes.begin(); it != nodes.end(); ++it)
		cells.push_back(it->first);
	return cells;
}

template<typename T>
inline std::vector<octree_id> octree<T>::all_nodes(int depth) const
{
	std::vector<octree_id> cells;
	for(typename node_storage::const_iterator it = nodes.begin(); it != nodes.end(); ++it)
		if(it->first.depth == depth)
			cells.push_back(it->first);
	return cells;
}

template<typename T>
inline std::vector<octree_id> octree<T>::all_min_depth_nodes(int depth) const
{
	std::vector<octree_id> cells;
	for(typename node_storage::const_iterator it = nodes.begin(); it != nodes.end(); ++it)
		if(it->first.depth >= depth)
			cells.push_back(it->first);
	return cells;
}

template<typename T>
inline std::vector<octree_id> octree<T>::all_boundary_nodes() const
{
	std::vector<octree_id> cells;
	for(typename node_storage::const_iterator it = nodes.begin(); it != nodes.end(); ++it)
		if(!it->first.is_interior())
			cells.push_back(it->first);
	return cells;
}

template<typename T>
inline void octree<T>::regularize(const T& def)
{
	sorted_nodes sn(*this);
	std::reverse(sn.begin(), sn.end()); // Important: this invalidates sn.depth_layers
	for(size_t t = 0; t < sn.size(); t++)
	{
		if(sn[t].is_root()) continue;
		const octree_id n = sn[t].parent();
		if(nodes.insert({ n, def }).second) sn.nodes.push_back(n);
		for(uint8_t i = 0; i < 8; i++)
			nodes.insert({ n.child(i), def });
	}
}

template<typename T>
inline void octree<T>::expand(int depth, const T& def)
{
	for(const octree_id& n : all_nodes(depth))
		for(const octree_id& m : neighbors(n))
			nodes.insert({ m, def });
}

inline Eigen::AlignedBox3f unit_octree::bounds()
{
	return Eigen::AlignedBox3f(Eigen::Vector3f::Zero(), Eigen::Vector3f::Ones());
}

inline Eigen::AlignedBox3f unit_octree::bounds(const octree_id& id)
{
	const float width = 1.0f / float(1 << id.depth);
	return Eigen::AlignedBox3f(id.pos.cast<float>() * width, (id.pos.cast<float>() + Eigen::Vector3f::Ones()) * width);
}

inline Eigen::Vector3f unit_octree::center(const octree_id& id)
{
	const float width = 1.0f / float(1 << id.depth);
	return (id.pos.cast<float>().array() + 0.5f) * width;
}

inline octree_id unit_octree::enclosing(const Eigen::Vector3f& pos, int depth)
{
	const int width = octree_id::width(depth);
	octree_id c;
	c.depth = depth;
	c.pos = (pos * float(width)).array().floor().matrix().cast<int>();
	return c;
}

inline Eigen::Vector3i node_color_by_depth(size_t, const octree_id& n)
{
	const Eigen::Vector3i cols[] =
	{

		{ 128, 128, 255 },
		{ 255, 128, 255 },
		{ 255, 0,   0   },
		{ 0,   255, 0   },
		{ 128, 255, 128 },
		{ 128, 255, 255 },
		{ 0,   0,   255 },
		{ 255, 128, 128 },
		{ 255, 255, 128 }
	};
	const size_t num_colors = (sizeof(cols) / sizeof(cols[0]));
	return cols[n.depth % num_colors];
};

template<typename Fc>
inline void save_octree_ply(std::ostream& out, const std::vector<octree_id>& nodes, Fc color,
	const Eigen::Matrix4f& transform = Eigen::Matrix4f::Identity(), bool use_empty_triangles = true)
{
	std::vector<Eigen::Vector3f> vertices;
	std::vector<size_t> vert_index;
	size_t lines = 0;
	auto add_line = [&](size_t index, const Eigen::Vector3f& a, const Eigen::Vector3f& b)
	{
		lines++;
		vertices.emplace_back(a[0], a[1], a[2]);
		vert_index.push_back(index);
		vertices.emplace_back(b[0], b[1], b[2]);
		vert_index.push_back(index);
		if(use_empty_triangles)
		{
			const Eigen::Vector3f m = (a + b) / 2.0f;
			vertices.emplace_back(m[0], m[1], m[2]);
			vert_index.push_back(index);
		}
	};
	int max_depth = -1;
	for(size_t t = 0; t < nodes.size(); t++)
	{
		const octree_id& c = nodes[t];
		max_depth = std::max(max_depth, c.depth);
		const Eigen::AlignedBox3f b = unit_octree::bounds(c);
		using corner = Eigen::AlignedBox3f::CornerType;
		add_line(t, b.corner(corner::BottomLeftFloor),  b.corner(corner::BottomRightFloor));
		add_line(t, b.corner(corner::BottomLeftFloor),  b.corner(corner::TopLeftFloor));
		add_line(t, b.corner(corner::TopLeftFloor),     b.corner(corner::TopRightFloor));
		add_line(t, b.corner(corner::BottomRightFloor), b.corner(corner::TopRightFloor));
		add_line(t, b.corner(corner::BottomLeftCeil),   b.corner(corner::BottomRightCeil));
		add_line(t, b.corner(corner::BottomLeftCeil),   b.corner(corner::TopLeftCeil));
		add_line(t, b.corner(corner::TopLeftCeil),      b.corner(corner::TopRightCeil));
		add_line(t, b.corner(corner::BottomRightCeil),  b.corner(corner::TopRightCeil));
		add_line(t, b.corner(corner::BottomLeftFloor),  b.corner(corner::BottomLeftCeil));
		add_line(t, b.corner(corner::BottomRightFloor), b.corner(corner::BottomRightCeil));
		add_line(t, b.corner(corner::TopLeftFloor),     b.corner(corner::TopLeftCeil));
		add_line(t, b.corner(corner::TopRightFloor),    b.corner(corner::TopRightCeil));
	}

	out << "ply\n"
		<< "format ascii 1.0\n"
		<< "comment Octree visualization (depth: " << max_depth << ")\n"
		<< "element vertex " << vertices.size() << "\n"
		<< "property float x\n"
		<< "property float y\n"
		<< "property float z\n"
		<< "property uchar red\n"
		<< "property uchar green\n"
		<< "property uchar blue\n";
	if(use_empty_triangles)
		out << "element face " << lines << "\n"
		<< "property list uchar int vertex_index\n";
	else
		out << "element edge " << lines << "\n"
		<< "property int vertex1\n"
		<< "property int vertex2\n";
	out << "end_header\n";

	for(size_t v = 0; v < vertices.size(); v++)
	{
		Eigen::Vector3f vert = vertices[v];
		Eigen::Vector4f vv = transform * Eigen::Vector4f(vert[0], vert[1], vert[2], 1.0f);
		vert = vv.block<3, 1>(0, 0) / vv[3];
		const Eigen::Vector3i col = color(vert_index[v], nodes[vert_index[v]]);
		out << vert[0] << ' ' << vert[1] << ' ' << vert[2] << ' '
			<< col[0] << ' ' << col[1] << ' ' << col[2] << '\n';
	}
	for(size_t l = 0; l < lines; l++)
	{
		const size_t c = use_empty_triangles ? 3 : 2;
		if(c != 2) out << c << ' ';
		for(size_t v = 0; v < c; v++)
		{
			if(v != 0) out << ' ';
			out << (l * c + v);
		}
		out << '\n';
	}
}

template<typename Fc>
inline void save_octree_ply(const std::string& filename, const std::vector<octree_id>& nodes, Fc color, const Eigen::Matrix4f& transform = Eigen::Matrix4f::Identity(), bool use_empty_triangles = true)
{
	std::ofstream out(filename);
	save_octree_ply(out, nodes, color, transform, use_empty_triangles);
}

template<typename T, typename Fc>
inline void save_octree_ply(const std::string& filename, const octree<T>& tree, Fc color, const Eigen::Matrix4f& transform = Eigen::Matrix4f::Identity(), bool use_empty_triangles = true)
{
	std::ofstream out(filename);
	save_octree_ply(out, tree.all_nodes(), color, transform, use_empty_triangles);
}

inline Eigen::AlignedBox3f aabb_octree::bounds() const
{
	return box;
}

inline Eigen::AlignedBox3f aabb_octree::bounds(const octree_id& id) const
{
	const Eigen::Array3f w = box.diagonal().array() / (1 << id.depth);
	const Eigen::Vector3f min = box.min().array() + w * id.pos.array().cast<float>();
	const Eigen::Vector3f max = box.min().array() + w * (id.pos.array() + 1).cast<float>();
	return { min, max };
}

inline Eigen::Vector3f aabb_octree::center(const octree_id& id) const
{
	return  bounds(id).center();
}

inline octree_id aabb_octree::enclosing(const Eigen::Vector3f& pos, int depth) const
{
	const float w = 1.0f / (1 << depth);
	const Eigen::Array3f P = (pos - box.min()).cwiseQuotient(box.diagonal()) / w;
	octree_id c;
	c.depth = depth;
	c.pos = P.floor().cast<int>();
	return c;
}
