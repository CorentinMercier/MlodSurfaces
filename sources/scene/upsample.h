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

#include <functional>
#include <random>
#include "../apss/Pn.h"

void regularize_up(std::vector<octree_id>& nodes);
void expand(std::vector<octree_id>& nodes, int depth = -1);
Eigen::Vector4i retreiveCorrectRoot(std::vector<octree_id>& nodes, glm::vec3 &transformDiag, glm::vec3 &transformVec) ;
Eigen::Vector4i union_nodes(std::vector<octree_id>& nodes);
void uniformize(unsigned int max_depth, std::vector<octree_id> leaf_nodes,
	sorted_octree<char>& soct, glm::vec3& new_diag, glm::vec3& new_origin, bool cut);
void uniformizeFromPoints(std::vector<octree_id> leaf_nodes, sorted_octree<char>& soct, glm::vec3& new_diag, glm::vec3& new_origin);
std::vector<octree_id> depth0_ancestors(std::vector<octree_id> nodes);

struct OctreeGenerator
{
	OctreeGenerator(APSS *apss, const Kernel& kernel) : apss(apss), kernel(kernel) {}
	std::vector<octree_id> generate(unsigned int max_depth,
		std::vector<octree_id>* total_nodes = nullptr) const;
	std::vector<octree_id> generateFromPoints(unsigned int max_depth,
		const Pn &pts, const Eigen::AlignedBox3f& box = Eigen::AlignedBox3f()) const;

private:
	std::vector<octree_id> compute_leaves(std::vector<octree_id> current_nodes, unsigned int max_depth, std::vector<octree_id>* total_nodes = nullptr) const;

	std::vector<octree_id> project(const octree_id& octId) const;

	void find_stable_nodes(const std::vector<octree_id>& nodes, std::vector<octree_id>& out) const;

	void find_stable_nodesGPU(const std::vector<octree_id>& nodes, std::vector<octree_id>& out) const;

	void add_children(const std::vector<octree_id>& nodes, std::vector<octree_id>& out) const;

	void expand_to_neighbors(std::vector<octree_id>& nodes, std::vector<octree_id>* total_tree) const;

private:
	APSS *apss;
	const Kernel& kernel;
};

namespace mesh
{
inline void smoothMesh(std::vector<glm::vec3> & iPoints, std::vector<uint32_t> & iTris)
{
	std::vector<glm::vec3> newPos(iPoints.size(), glm::vec3(0,0,0));
	std::vector<float> valence(iPoints.size(), 0);

	const unsigned int iTrisNb = iTris.size() / 3;
	for( unsigned int t = 0 ; t < iTrisNb ; ++t )
	{
		unsigned int a = iTris[3*t];
		unsigned int b = iTris[3*t + 1];
		unsigned int c = iTris[3*t + 2];

		newPos[a] += iPoints[b] + iPoints[c];
		newPos[b] += iPoints[a] + iPoints[c];
		newPos[c] += iPoints[b] + iPoints[a];
		valence[a] += 2;
		valence[b] += 2;
		valence[c] += 2;
	}
#pragma omp parallel for
	for( unsigned int v = 0 ; v < iPoints.size() ; ++v )
		iPoints[v] = newPos[v] / valence[v];
}

inline void subdivideMesh(std::vector<glm::vec3>& iPoints, std::vector<uint32_t>& iTris)
{
	struct Edge{
		unsigned int a , b;
		Edge(){}
		Edge(unsigned int c , unsigned int d) : a( std::min(c,d) ) , b( std::max(c,d) ) {}
		bool operator < (Edge const & o) const { return a < o.a || (a == o.a  &&  b < o.b) ; }
	};
	std::map< Edge , unsigned int > newVertOnEdge;
	const unsigned int iTrisNb = iTris.size() / 3;
	for( unsigned int t = 0 ; t < iTrisNb ; ++t )
	{
		unsigned int a = iTris[3*t];
		unsigned int b = iTris[3*t + 1];
		unsigned int c = iTris[3*t + 2];
		unsigned int ab , bc , ca;
		std::map< Edge , unsigned int >::iterator ab_vert = newVertOnEdge.find( Edge(a,b) );
		if( ab_vert ==  newVertOnEdge.end() ) {
			iPoints.push_back( (iPoints[a] + iPoints[b])/2.f );
			ab = iPoints.size() - 1;
			newVertOnEdge[ Edge(a,b) ] = ab;
		}
		else ab = ab_vert->second;
		std::map< Edge , unsigned int >::iterator bc_vert = newVertOnEdge.find( Edge(c,b) );
		if( bc_vert ==  newVertOnEdge.end() ) {
			iPoints.push_back( (iPoints[c] + iPoints[b])/2.f );
			bc = iPoints.size() - 1;
			newVertOnEdge[ Edge(c,b) ] = bc;
		}
		else bc = bc_vert->second;
		std::map< Edge , unsigned int >::iterator ca_vert = newVertOnEdge.find( Edge(c,a) );
		if( ca_vert ==  newVertOnEdge.end() ) {
			iPoints.push_back( (iPoints[c] + iPoints[a])/2.f );
			ca = iPoints.size() - 1;
			newVertOnEdge[ Edge(c,a) ] = ca;
		}
		else ca = ca_vert->second;

		iTris[3*t] = a;
		iTris[3*t+1] = ab;
		iTris[3*t+2] = ca;

		iTris.push_back( ab );
		iTris.push_back( b );
		iTris.push_back( bc );

		iTris.push_back( ab );
		iTris.push_back( bc );
		iTris.push_back( ca );

		iTris.push_back( ca );
		iTris.push_back( bc );
		iTris.push_back( c );
	}
}
}
