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

#include "Kernel.h"

#include <set>
#include <queue>
#include "../sources/dual-contouring/kdTreeNN/kdTreeNN.h"

#define MAX_DEPTH 12

//#define SNDORDER

///
/// CPU structs
///

// all points in the node are assumed to be given the same weight w:
struct apssNodeStats{
	glm::vec3 s_ai_pi , s_ai_ni;
	float s_ai_pi_ni , s_ai_pi_pi;
	float s_ai;
#ifdef SNDORDER
	glm::mat3 s_ai_pi_piT, s_ai_ni_piT;
	glm::vec3 s_ai_pi_pi_pi, s_ai_pi_ni_pi;

	apssNodeStats() : s_ai_pi(0,0,0) , s_ai_ni(0,0,0) , s_ai_pi_ni(0) , s_ai_pi_pi(0) , s_ai(0) , s_ai_pi_piT(0) , s_ai_ni_piT(0) , s_ai_pi_pi_pi(0) , s_ai_pi_ni_pi(0){}
	apssNodeStats( glm::vec3 const & p , glm::vec3 const & n ) : s_ai_pi(p) , s_ai_ni(n) , s_ai_pi_ni( glm::dot(p,n) ), s_ai_pi_pi( glm::dot(p,p) ) , s_ai(1) ,
		s_ai_pi_piT(glm::outerProduct(p,p)) , s_ai_ni_piT(glm::outerProduct(n,p)) , s_ai_pi_pi_pi(glm::dot(p,p)*p) , s_ai_pi_ni_pi(glm::dot(p,n)*n) {}
	apssNodeStats( glm::vec3 const & p , glm::vec3 const & n , float a ) : s_ai_pi(a*p) , s_ai_ni(a*n) , s_ai_pi_ni( a*glm::dot(p,n) ), s_ai_pi_pi( a*glm::dot(p,p) ) , s_ai(a) ,
		s_ai_pi_piT(a*glm::outerProduct(p,p)) , s_ai_ni_piT(a*glm::outerProduct(n,p)) , s_ai_pi_pi_pi(a*glm::dot(p,p)*p) , s_ai_pi_ni_pi(a*glm::dot(p,n)*n) {}
#else
	apssNodeStats() : s_ai_pi(0,0,0) , s_ai_ni(0,0,0) , s_ai_pi_ni(0) , s_ai_pi_pi(0) , s_ai(0) {}
	apssNodeStats( glm::vec3 const & p , glm::vec3 const & n ) : s_ai_pi(p) , s_ai_ni(n) , s_ai_pi_ni( glm::dot(p,n) ), s_ai_pi_pi( glm::dot(p,p) ) , s_ai(1) {}
	apssNodeStats( glm::vec3 const & p , glm::vec3 const & n , float a ) : s_ai_pi(a*p) , s_ai_ni(a*n) , s_ai_pi_ni( a*glm::dot(p,n) ), s_ai_pi_pi( a*glm::dot(p,p) ) , s_ai(a) {}
#endif

	void operator = (apssNodeStats const & o);
};

struct apssStats{
	glm::vec3 s_ai_wi_pi , s_ai_wi_ni;
	float s_ai_wi, s_ai_wi_pi_ni , s_ai_wi_pi_pi;

	apssStats() : s_ai_wi_pi(0,0,0) , s_ai_wi_ni(0,0,0) , s_ai_wi(0) , s_ai_wi_pi_ni(0) , s_ai_wi_pi_pi(0) {}

	void operator += (apssStats const & o);
	apssStats operator + (apssStats const & o);
	void operator /= (float const & f);
	apssStats operator * (float const & f);
};

#ifdef SNDORDER
inline apssStats operator*(const apssNodeStats& ns, const glm::vec4& gradWAndW);
#endif

inline apssStats operator*(const apssNodeStats& ns, float w);


struct BBOX // To move somewhere else
{
	glm::vec3 bb,BB;

	BBOX() { clear(); }

	void clear();
	bool isCleared() const { return bb[0] <= BB[0]; }
	void set( const glm::vec3 & p );
	void set( const glm::vec3 & pbb, const glm::vec3 & pBB );
	void add( const glm::vec3 & p );
	void add( const BBOX & b );

	float squareDiagonal() const { return glm::distance2(BB, bb); }
	inline float diagonal() const { return sqrt( glm::distance2(BB, bb)); }
	inline float radius() const { return diagonal() / 2.f; }
	inline glm::vec3 center() const { return (bb + BB) / 2.f; }
	inline float squareRadius() const { return squareDiagonal() / 4.f; }

	char getLargestExtent() const;
	inline float getExtentValue(char i) const { return BB[i] - bb[i]; }
	inline float getLargestExtentValue() const { return getExtentValue( getLargestExtent() ); }

	inline float getPseudoExtentInDirection( glm::vec3 const & dir ) const {
		return getExtentValue(0) * fabs(dir[0]) +getExtentValue(1) * fabs(dir[1]) +getExtentValue(2) * fabs(dir[2]);
	}

	inline void splitAlongAxis( char axis , double value , BBOX & bbox1 , BBOX & bbox2 );
	glm::vec3 corner(unsigned int i) const;
};

struct relations
{
	int father = -1;
	int firstChild = -1;
	int nextBrother = -1;
};

struct octreeNode {
	// breadth-first indexing:
	unsigned int breadth_first_index;

	octreeNode * children[8];
	unsigned int numberOfChildren;
	std::vector< unsigned int > indices;
	unsigned int depth;
	BBOX boundingBox;

	void setBoundingBox( const BBOX & bbox ) { boundingBox = bbox; }
	BBOX const & getBoundingBox() const { return boundingBox; }

	apssNodeStats nodeapssNodeStats;

	octreeNode();

	unsigned int nbOfNodes() const;

	void copy(std::vector<octreeNode const *> &nodes, std::vector<relations> &related, int & nodePos) const;

	void getRelations(std::vector<octreeNode const *> &nodes, std::vector<relations> &related) const;

	void estimateAreasRecursive(std::vector<float> & areas);

	unsigned int computeDepth() const;

	void deleteRecursive();

	void buildRecursive(
			const std::vector< glm::vec3 > & i_points,
			const std::vector< unsigned int > & i_indices,
			unsigned int currentDepth,
			unsigned int maxDepth,
			unsigned int maxNumberOfPointsPerLeaf);

	const std::vector<const octreeNode *> getLeaves() const;
	const std::vector<const octreeNode *> all_nodes() const;

	bool isLeaf() const { return indices.size() > 0 ; }

	void buildapssNodeStatsRecursive( std::vector< apssNodeStats > const & leavesapssNodeStats );

	float blendingFunc(float x)const;

	void gatherTraversedNodes( glm::vec3 const & q , std::set<unsigned int> & traversed_nodes , unsigned int minimal_depth , float scalingProtectionSphere ) const;

	// this will build the statistics required to fit an algebraic sphere to the octree:
	apssStats approximateAPSSStatsRecursive(glm::vec3 const & q , Kernel const & kernel , unsigned int minimal_depth , float scalingProtectionSphere ,
											 std::vector< apssNodeStats > const & leavesapssNodeStats) const;

	apssStats APSSStats(unsigned nbOfLeaves, glm::vec3 const & q, Kernel const & kernel,
						std::vector< apssNodeStats > const & leavesapssNodeStats) const;
};

class GlobalAPSS;
class Fastapss;

class APSSOctree {
	friend class GlobalAPSS;
	friend class Fastapss;

public:
	APSSOctree(){m_root = new octreeNode();}
	virtual ~APSSOctree(){m_root->deleteRecursive();}

	void printState() const;

	void build( const std::vector< glm::vec3 > & i_points ,
				std::vector< glm::vec3 > const & normals ,
				bool updateRootBoundingBox = true , // sometimes you don t have all the points when you build, but you want to setup a fixed bounding box
				unsigned int maxDepth = MAX_DEPTH,
				unsigned int maxNumberOfPointsPerLeaf = 1);

	BBOX const & getBoundingBox() const { return m_root->getBoundingBox();}
	const octreeNode* get_root() const { return m_root; }
	const std::vector< apssNodeStats >& getPointsApssNodeStats() const { return m_pointsapssNodeStats; }
	float getScalingProtectionSphere() const {return m_scalingProtectionSphere;}
	unsigned getNbOfLeaves() const{return m_numberOfLeaves;}
	unsigned getMinDepth(){return m_minDepth;}
	void setScalingProtectionSphere(float val){m_scalingProtectionSphere = val;}
	void setMinDepth(unsigned val){m_minDepth = val;}

	void knnMode(bool setMode, unsigned knn = 20);
	void ballMode(bool setMode, float ballRadius = 0.1);
	bool isAPSS(){return (m_knnMode || m_ballMode);}

	void cpuSphere(glm::vec3 const & qStart, float &u0, glm::vec3 &u123, float &u4, Kernel const & kernel = Kernel()) const;
	void getKnn(glm::vec3 const & x, unsigned k, std::vector<glm::vec3> & neighbors) const;
protected:
	void projectCPU(glm::vec3 const & qStart, glm::vec3 & outputPoint, glm::vec3 & outputNormal,
					unsigned int n_iterations, bool fast, Kernel const & kernel = Kernel() , float stepMaxSize = -1.f) const;

private:
	octreeNode *m_root;
	std::vector<glm::vec3> m_points;
	std::vector< apssNodeStats > m_pointsapssNodeStats;
	unsigned m_numberOfLeaves;
	float m_scalingProtectionSphere = 2.0;
	unsigned m_minDepth = 0;

	//APSS knn
	bool m_knnMode = false;
	unsigned m_knn = 20;
	//APSS in ball
	bool m_ballMode = false;
	float m_ballRadius = 0.1;
	kdTreeNN::kdTree<glm::vec3> m_kdTree;
	std::vector<glm::vec3> m_normals;

	void computeBoundingBox();
	void estimateAreas(std::vector<float>& areas);
	void buildapssNodeStats( std::vector< glm::vec3 > const & normals , std::vector< float > const & areas );
	void index_nodes_by_breadth_first_search();

	void oneStepProjection(glm::vec3 const & q, glm::vec3 & outputPoint , glm::vec3 & outputNormal , bool fast,
							Kernel const & kernel , unsigned int minimal_depth , float scalingProtectionSphere , float step_max_size ) const;
	void projectOnAPSSStats(glm::vec3 const & q,  glm::vec3 & outputPoint , glm::vec3 & outputNormal , apssStats const & treeStats , float step_max_size ) const;

	void projectAPSSCPU(glm::vec3 const & qStart, glm::vec3 & outputPoint, glm::vec3 & outputNormal,
						unsigned int n_iterations, Kernel const & kernel = Kernel() , float stepMaxSize = -1.f) const;
};
