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

#include "CPUUtils.h"

///
/// struct apssNodeStats
///

void apssNodeStats::operator = (apssNodeStats const & o) {
	s_ai_pi = o.s_ai_pi;
	s_ai_ni = o.s_ai_ni;
	s_ai_pi_ni = o.s_ai_pi_ni;
	s_ai_pi_pi = o.s_ai_pi_pi;
	s_ai = o.s_ai;
#ifdef SNDORDER
	s_ai_pi_piT = o.s_ai_pi_piT;
	s_ai_ni_piT = o.s_ai_ni_piT;
	s_ai_pi_pi_pi = o.s_ai_pi_pi_pi;
	s_ai_pi_ni_pi = o.s_ai_pi_ni_pi;
#endif
}

///
/// struct apssStats
///

void apssStats::operator += (apssStats const & o) {
	s_ai_wi_pi += o.s_ai_wi_pi;
	s_ai_wi_ni += o.s_ai_wi_ni;
	s_ai_wi += o.s_ai_wi;
	s_ai_wi_pi_ni += o.s_ai_wi_pi_ni;
	s_ai_wi_pi_pi += o.s_ai_wi_pi_pi;
}

apssStats apssStats::operator + (apssStats const & o) {
	apssStats stats;
	stats.s_ai_wi_pi = s_ai_wi_pi + o.s_ai_wi_pi;
	stats.s_ai_wi_ni = s_ai_wi_ni + o.s_ai_wi_ni;
	stats.s_ai_wi = s_ai_wi + o.s_ai_wi;
	stats.s_ai_wi_pi_ni = s_ai_wi_pi_ni + o.s_ai_wi_pi_ni;
	stats.s_ai_wi_pi_pi = s_ai_wi_pi_pi + o.s_ai_wi_pi_pi;
	return stats;
}

void apssStats::operator /= (float const & f) {
	s_ai_wi_pi /= f;
	s_ai_wi_ni /= f;
	s_ai_wi /= f;
	s_ai_wi_pi_ni /= f;
	s_ai_wi_pi_pi /= f;
}

apssStats apssStats::operator * (float const & f) {
	apssStats stats;
	stats.s_ai_wi_pi = f * s_ai_wi_pi;
	stats.s_ai_wi_ni = f * s_ai_wi_ni;
	stats.s_ai_wi = s_ai_wi * f;
	stats.s_ai_wi_pi_ni = s_ai_wi_pi_ni * f;
	stats.s_ai_wi_pi_pi = s_ai_wi_pi_pi * f;
	return stats;
}

#ifdef SNDORDER
inline apssStats operator*(const apssNodeStats& ns, const glm::vec4& gradWAndW)
{
	apssStats ts;
	glm::vec3 p = ns.s_ai_pi / ns.s_ai;
	ts.s_ai_wi = gradWAndW.w * ns.s_ai;
	ts.s_ai_wi_ni = gradWAndW.w * ns.s_ai_ni + glm::vec3(gradWAndW) * (ns.s_ai_ni_piT - glm::outerProduct(ns.s_ai_ni ,p));
	ts.s_ai_wi_pi = gradWAndW.w * ns.s_ai_pi + glm::vec3(gradWAndW) * (ns.s_ai_pi_piT - glm::outerProduct(ns.s_ai_pi,p));
	ts.s_ai_wi_pi_ni = gradWAndW.w * ns.s_ai_pi_ni + glm::dot(glm::vec3(gradWAndW), (ns.s_ai_pi_ni_pi - ns.s_ai_pi_ni * p));
	ts.s_ai_wi_pi_pi = gradWAndW.w * ns.s_ai_pi_pi + glm::dot(glm::vec3(gradWAndW), (ns.s_ai_pi_pi_pi - ns.s_ai_pi_pi * p));
	return ts;
}
#endif

inline apssStats operator*(const apssNodeStats& ns, float w)
{
	apssStats ts;
	ts.s_ai_wi = w * ns.s_ai;
	ts.s_ai_wi_ni = w * ns.s_ai_ni;
	ts.s_ai_wi_pi = w * ns.s_ai_pi;
	ts.s_ai_wi_pi_ni = w * ns.s_ai_pi_ni;
	ts.s_ai_wi_pi_pi = w * ns.s_ai_pi_pi;
	return ts;
}

///
/// struct BBOX
///

void BBOX::clear()
{
	bb[0] = FLT_MAX;
	BB[0] = -FLT_MAX;
}

void BBOX::set( const glm::vec3 & p )
{
	bb = glm::vec3(p[0],p[1],p[2]);
	BB = glm::vec3(p[0],p[1],p[2]);
}

void BBOX::set( const glm::vec3 & pbb, const glm::vec3 & pBB )
{
	bb = glm::vec3(pbb[0],pbb[1],pbb[2]);
	BB = glm::vec3(pBB[0],pBB[1],pBB[2]);
}

void BBOX::add( const glm::vec3 & p )
{
	for( unsigned int c = 0 ; c < 3 ; ++c ) {
		bb[c] = std::min( bb[c],p[c] );
		BB[c] = std::max( BB[c],p[c] );
	}
}

void BBOX::add( const BBOX & b )
{
	for( unsigned int c = 0 ; c < 3 ; ++c ) {
		bb[c] = std::min( bb[c],b.bb[c] );
		BB[c] = std::max( BB[c],b.BB[c] );
	}
}

inline char BBOX::getLargestExtent() const {
	if(  BB[0] - bb[0]  >  BB[1] - bb[1]  ) {
		if( BB[0] - bb[0] > BB[2] - bb[2] )
			return 0;
		return 2;
	}
	else {
		if( BB[1] - bb[1] > BB[2] - bb[2] )
			return 1;
		return 2;
	}
}

inline void BBOX::splitAlongAxis( char axis , double value , BBOX & bbox1 , BBOX & bbox2 ) {
	// for safety:
	value = std::max( std::min( (float)value , BB[axis] ) , bb[axis] );
	glm::vec3 BB1 = BB;
	BB1[axis] = value;
	glm::vec3 bb2 = bb;
	bb2[axis] = value;
	bbox1.set( bb , BB1 );
	bbox2.set( bb2 , BB );
}

glm::vec3 BBOX::corner(unsigned int i) const
{
	glm::vec3 c = bb;
	if((i >> 0) % 2 == 1) c[0] = BB[0];
	if((i >> 1) % 2 == 1) c[1] = BB[1];
	if((i >> 2) % 2 == 1) c[2] = BB[2];
	return c;
}

///
/// struct octreeNode
///

octreeNode::octreeNode() {
	for(unsigned int c = 0 ; c < 8 ; ++c) children[c] = NULL;
	numberOfChildren = 0;
}

unsigned int octreeNode::nbOfNodes() const
{
	unsigned int num = 1;
	for (unsigned int c=0; c<8; c++)
		if (children[c])
			num += children[c]->nbOfNodes();
	return num;
}

void octreeNode::copy(std::vector<octreeNode const *> &nodes, std::vector<relations> &related, int & nodePos) const
{
	nodePos++;
	nodes.push_back(this);
	related[nodePos].firstChild = -1;
	related[nodePos].nextBrother = -1;
	bool firstChild = true;
	int precedentChild = -1;
	int fatherPos = nodePos;
	for (unsigned int c=0; c<8; c++)
		if (children[c])
		{
			if (firstChild)
			{
				related[nodePos].firstChild = nodePos + 1;
				related[nodePos+1].father = nodePos;
				firstChild = false;
			}
			else {
				related[precedentChild].nextBrother = nodePos + 1;
				related[nodePos+1].father = fatherPos;
			}
			precedentChild = nodePos + 1;
			children[c]->copy(nodes, related, nodePos);
		}
}

void octreeNode::getRelations(std::vector<octreeNode const *> &nodes, std::vector<relations> &related) const
{
	nodes[breadth_first_index] = this;
	bool firstChild = true;
	int precedentChild = -1;
	for (unsigned int c=0; c<8; c++)
		if (children[c] != nullptr)
		{
			int childIndex = int(children[c]->breadth_first_index);
			related[childIndex].father = int(breadth_first_index);
			if (firstChild)
			{
				related[breadth_first_index].firstChild = childIndex;
				firstChild = false;
			}
			else {
				related[precedentChild].nextBrother = childIndex;
			}
			precedentChild = childIndex;
			children[c]->getRelations(nodes, related);
		}
}

void octreeNode::estimateAreasRecursive(std::vector<float> & areas) {
	// traverse the octree, consider only the leaves, and for each leaf:
	//    compute an area for the leaf, and split the area among the point samples inside it
	if( isLeaf() ) {
		// compute the area:
		float leafArea = boundingBox.diagonal() * boundingBox.diagonal();
		float sampleArea = leafArea / indices.size();
		for( unsigned int idx : indices ) {
			areas[ idx ] = sampleArea;
		}
	}
	else {
		for (unsigned int c=0; c<8; c++)
			if (children[c])
				children[c]->estimateAreasRecursive(areas);
	}
}

unsigned int octreeNode::computeDepth() const
{
	unsigned d = depth;
	for (unsigned c=0; c<8; c++)
		if (children[c])
			d = std::max(d, children[c]->computeDepth());
	return d;
}

void octreeNode::deleteRecursive()
{
	for (unsigned int c=0; c<8; c++)
		if (children[c])
			children[c]->deleteRecursive();
	delete(this);
}

void octreeNode::buildRecursive(
		const std::vector< glm::vec3 > & i_points,
		const std::vector< unsigned int > & i_indices,
		unsigned int currentDepth,
		unsigned int maxDepth,
		unsigned int maxNumberOfPointsPerLeaf) {
	depth = currentDepth;

	if( i_indices.size() <= maxNumberOfPointsPerLeaf   ||   currentDepth == maxDepth ) {
		indices = i_indices;
		assert(i_indices.size() > 0);
		return;
	}

	std::vector< unsigned int > childrenIndices [8];
	double xSplit = (boundingBox.bb[0] + boundingBox.BB[0])/2.0;
	double ySplit = (boundingBox.bb[1] + boundingBox.BB[1])/2.0;
	double zSplit = (boundingBox.bb[2] + boundingBox.BB[2])/2.0;

	double xVals [3] = {boundingBox.bb[0] , xSplit , boundingBox.BB[0]};
	double yVals [3] = {boundingBox.bb[1] , ySplit , boundingBox.BB[1]};
	double zVals [3] = {boundingBox.bb[2] , zSplit , boundingBox.BB[2]};

	for( unsigned int pIt = 0 ; pIt < i_indices.size() ; ++pIt ) {
		unsigned int pIdx = i_indices[pIt];
		glm::vec3 const & p = i_points[pIdx];

		unsigned int xi = 0;
		if( p[0] == xSplit ) if( rand() > RAND_MAX / 2 ) xi = 1;
		if( p[0] > xSplit ) xi = 1;

		unsigned int yi = 0;
		if( p[1] == ySplit ) if( rand() > RAND_MAX / 2 ) yi = 1;
		if( p[1] > ySplit ) yi = 1;

		unsigned int zi = 0;
		if( p[2] == zSplit ) if( rand() > RAND_MAX / 2 ) zi = 1;
		if( p[2] > zSplit ) zi = 1;

		childrenIndices[ xi + 2 * yi + 4 * zi ].push_back( pIdx );
	}

	numberOfChildren = 0;
	for( unsigned int xi = 0 ; xi < 2 ; ++xi ) {
		for( unsigned int yi = 0 ; yi < 2 ; ++yi ) {
			for( unsigned int zi = 0 ; zi < 2 ; ++zi ) {
				unsigned int c = xi + 2 * yi + 4 * zi;
				if( childrenIndices[c].size() > 0 ) {
					++numberOfChildren;
					children[c] = new octreeNode;
					BBOX bbox;
					bbox.set( glm::vec3( xVals[xi] , yVals[yi] , zVals[zi] )  ,  glm::vec3( xVals[xi+1] , yVals[yi+1] , zVals[zi+1] ) );
					children[c]->setBoundingBox(bbox);
					children[c]->buildRecursive(i_points , childrenIndices[c] , currentDepth + 1 , maxDepth , maxNumberOfPointsPerLeaf );
				}
			}
		}
	}
}

void append_vector(std::vector<const octreeNode*>& vec, const std::vector<const octreeNode*>& to_add)
{
	vec.insert(vec.end(), to_add.begin(), to_add.end());
}

const std::vector<const octreeNode*> octreeNode::getLeaves() const
{
	std::vector<const octreeNode*> leaves;
	if (isLeaf())
		leaves.push_back(this);
	else
		for (unsigned c=0; c<8; c++)
			if (children[c])
				append_vector(leaves, children[c]->getLeaves());
	return leaves;
}

const std::vector<const octreeNode*> octreeNode::all_nodes() const
{
	std::vector<const octreeNode*> nodes;
	nodes.push_back(this);
	for (unsigned c=0; c<8; c++)
		if (children[c])
			append_vector(nodes, children[c]->all_nodes());
	return nodes;
}

void octreeNode::buildapssNodeStatsRecursive( std::vector< apssNodeStats > const & leavesapssNodeStats ) {
	nodeapssNodeStats = apssNodeStats(); // make sure it's cleared
	if( isLeaf() ) {
		for( unsigned int lIt = 0 ; lIt < indices.size() ; ++lIt ) {
			unsigned int leafIdx = indices[lIt];
			apssNodeStats const & lapssNodeStats = leavesapssNodeStats[ leafIdx ];

			nodeapssNodeStats.s_ai += lapssNodeStats.s_ai;
			nodeapssNodeStats.s_ai_ni += lapssNodeStats.s_ai_ni;
			nodeapssNodeStats.s_ai_pi += lapssNodeStats.s_ai_pi;
			nodeapssNodeStats.s_ai_pi_ni += lapssNodeStats.s_ai_pi_ni;
			nodeapssNodeStats.s_ai_pi_pi += lapssNodeStats.s_ai_pi_pi;
#ifdef SNDORDER
			nodeapssNodeStats.s_ai_ni_piT += lapssNodeStats.s_ai_ni_piT;
			nodeapssNodeStats.s_ai_pi_piT += lapssNodeStats.s_ai_pi_piT;
			nodeapssNodeStats.s_ai_pi_ni_pi += lapssNodeStats.s_ai_pi_ni_pi;
			nodeapssNodeStats.s_ai_pi_pi_pi += lapssNodeStats.s_ai_pi_pi_pi;
#endif
		}
	}
	else{
		for( unsigned int c = 0 ; c < 8 ; ++c ) {
			octreeNode * child = children[c];
			if( child != NULL ) {
				child->buildapssNodeStatsRecursive( leavesapssNodeStats );
				apssNodeStats const & lapssNodeStats = child->nodeapssNodeStats;
				nodeapssNodeStats.s_ai += lapssNodeStats.s_ai;
				nodeapssNodeStats.s_ai_ni += lapssNodeStats.s_ai_ni;
				nodeapssNodeStats.s_ai_pi += lapssNodeStats.s_ai_pi;
				nodeapssNodeStats.s_ai_pi_ni += lapssNodeStats.s_ai_pi_ni;
				nodeapssNodeStats.s_ai_pi_pi += lapssNodeStats.s_ai_pi_pi;
#ifdef SNDORDER
			nodeapssNodeStats.s_ai_ni_piT += lapssNodeStats.s_ai_ni_piT;
			nodeapssNodeStats.s_ai_pi_piT += lapssNodeStats.s_ai_pi_piT;
			nodeapssNodeStats.s_ai_pi_ni_pi += lapssNodeStats.s_ai_pi_ni_pi;
			nodeapssNodeStats.s_ai_pi_pi_pi += lapssNodeStats.s_ai_pi_pi_pi;
#endif
			}
		}
	}
}

float octreeNode::blendingFunc(float x)const
{
#ifdef EXPFUNC
	if (x>=1) return 1;
	return exp(-(exp(1 / (x-1))) / (x*x));
#else
	return (1 - cos(3.14159265 * x))/2.0;
#endif
}

void octreeNode::gatherTraversedNodes( glm::vec3 const & q , std::set<unsigned int> & traversed_nodes , unsigned int minimal_depth , float scalingProtectionSphere ) const
{
	traversed_nodes.insert( this->breadth_first_index );
	if( isLeaf() ) {
		return ;
	}
	// otherwise, check if you're sufficiently far to allow for simple approximation:
	if( depth >= minimal_depth )
	{
		if( glm::length(q - boundingBox.center())  > scalingProtectionSphere * boundingBox.radius() )
		{
			return ;
		}
		else { // do blending between levels
			for( unsigned int c = 0 ; c < 8 ; ++c ) {
				octreeNode * child = children[c];
				if( child != NULL ) {
					child->gatherTraversedNodes( q , traversed_nodes , minimal_depth , scalingProtectionSphere );
				}
			}
			return ;
		}
	}
	// if not, then just accumulate the wn of the children:
	for( unsigned int c = 0 ; c < 8 ; ++c ) {
		octreeNode * child = children[c];
		if( child != NULL ) {
			child->gatherTraversedNodes( q , traversed_nodes , minimal_depth , scalingProtectionSphere );
		}
	}
}

// this will build the statistics required to fit an algebraic sphere to the octree:
apssStats octreeNode::approximateAPSSStatsRecursive( glm::vec3 const & q , Kernel const & kernel , unsigned int minimal_depth , float scalingProtectionSphere ,
										 std::vector< apssNodeStats > const & leavesapssNodeStats ) const {
	apssStats treeStats; // initialize with nothing, we will accumulate along the tree

	if( isLeaf() ) {
		for( unsigned int lIt = 0 ; lIt < indices.size() ; ++lIt ) {
			unsigned int leafIdx = indices[lIt];
			apssNodeStats const & lapssNodeStats = leavesapssNodeStats[ leafIdx ];
#ifdef SNDORDER
			glm::vec4 w = kernel.gradW( lapssNodeStats.s_ai_pi / lapssNodeStats.s_ai , q ); // evaluate weight wrt node average point
			treeStats = lapssNodeStats * w;
#else
			float w = kernel.w( lapssNodeStats.s_ai_pi / lapssNodeStats.s_ai , q ); // evaluate weight wrt node average point
			treeStats = lapssNodeStats * w;
#endif
		}
		return treeStats;
	}

	// otherwise, check if you're sufficiently far to allow for simple approximation:
	if( depth >= minimal_depth )
	{
		apssNodeStats const & lapssNodeStats = this->nodeapssNodeStats;
		glm::vec3 p = lapssNodeStats.s_ai_pi / lapssNodeStats.s_ai;
		if( glm::length(q - boundingBox.center())  > scalingProtectionSphere * boundingBox.radius() )
		{
#ifdef SNDORDER
			glm::vec4 w = kernel.gradW( p , q ); // evaluate weight wrt node average point
			treeStats = lapssNodeStats * w;
#else
			float w = kernel.w( p , q ); // evaluate weight wrt node average point
			treeStats = lapssNodeStats * w;
#endif
			return treeStats;
		}
		else {
			// do blending between levels
			apssStats parentStats;
#ifdef SNDORDER
			glm::vec4 w = kernel.gradW( p , q ) / float(numberOfChildren); // evaluate weight wrt node average point
			parentStats = lapssNodeStats * w;
#else
			float w = kernel.w( p , q ) / float(numberOfChildren);
			parentStats = lapssNodeStats * w;
#endif

			float distToParent = fabs(glm::length(q - boundingBox.center()) - (boundingBox.radius() * scalingProtectionSphere));
			for( unsigned int c = 0 ; c < 8 ; ++c ) {
				octreeNode * child = children[c];
				if( child != NULL ) {
					float wChild = distToParent / (distToParent + std::max<float>(glm::length(q - child->boundingBox.center()) - (child->boundingBox.radius() * scalingProtectionSphere), 0.f)); //kernel.w( p , q );
					apssStats childStats;
					childStats = child->approximateAPSSStatsRecursive(q , kernel , minimal_depth , scalingProtectionSphere , leavesapssNodeStats);
					float valChild = blendingFunc(wChild);
					float valParent = 1 - valChild;
					treeStats += childStats * valChild + parentStats * valParent;
				}
			}
			return treeStats;
		}
	}

	// if not, then just accumulate the wn of the children:
	for( unsigned int c = 0 ; c < 8 ; ++c ) {
		octreeNode * child = children[c];
		if( child != NULL ) {
			treeStats += child->approximateAPSSStatsRecursive(q , kernel , minimal_depth , scalingProtectionSphere , leavesapssNodeStats);
		}
	}
	return treeStats;
}

apssStats octreeNode::APSSStats(unsigned nbOfLeaves, const glm::vec3 &q, const Kernel &kernel, const std::vector<apssNodeStats> &leavesapssNodeStats) const
{
	apssStats treeStats;
	for (unsigned i=0; i<nbOfLeaves; i++)
	{
		float w = kernel.w(leavesapssNodeStats[i].s_ai_pi / leavesapssNodeStats[i].s_ai, q);
		treeStats.s_ai_wi = w * leavesapssNodeStats[i].s_ai + treeStats.s_ai_wi;
		treeStats.s_ai_wi_ni = treeStats.s_ai_wi_ni + w * leavesapssNodeStats[i].s_ai_ni;
		treeStats.s_ai_wi_pi = treeStats.s_ai_wi_pi + w * leavesapssNodeStats[i].s_ai_pi;
		treeStats.s_ai_wi_pi_ni = w * leavesapssNodeStats[i].s_ai_pi_ni + treeStats.s_ai_wi_pi_ni;
		treeStats.s_ai_wi_pi_pi = w * leavesapssNodeStats[i].s_ai_pi_pi + treeStats.s_ai_wi_pi_pi;
	}
	return treeStats;
}

///
/// class APSSOctree
///

void APSSOctree::printState() const
{
	std::cout << "Size of pointset : " << m_points.size() << std::endl;
	std::cout << "Size of octree : " << m_root->nbOfNodes() << std::endl;
}

void APSSOctree::build( const std::vector< glm::vec3 > & i_points ,
					  std::vector< glm::vec3 > const & normals ,
					  bool updateRootBoundingBox, // sometimes you don t have all the points when you build, but you want to setup a fixed bounding box
					  unsigned int maxDepth,
					  unsigned int maxNumberOfPointsPerLeaf)
{
	assert( i_points.size() > 0  &&  "Come ooooon" );
	m_points = i_points;
	if(updateRootBoundingBox) computeBoundingBox();
	std::vector< unsigned int > indices( m_points.size() );
	for( unsigned int t = 0 ; t < m_points.size() ; ++t )
		indices[t] = t;

	//std::cout << "Building tree" << std::endl;
	m_root->buildRecursive(m_points, indices , 0 , maxDepth , maxNumberOfPointsPerLeaf);
	//std::cout << "Tree built" << std::endl;

	std::vector< float > sampleAreas(m_points.size());
	estimateAreas( sampleAreas );

	buildapssNodeStats(normals , sampleAreas);

	index_nodes_by_breadth_first_search(); // Is that useful ?

	//Kd tree for APSS
	//std::cout << "Building kd tree for APSS knn" << std::endl;
	m_normals = normals;
	m_kdTree.build(i_points);
	//std::cout << "Kd tree built" << std::endl;
}

void APSSOctree::computeBoundingBox()
{
	assert( m_points.size() > 0  &&  "Come ooooon" );
	BBOX boundingBox;
	boundingBox.set(m_points[0]);
	for( unsigned int t = 1 ; t < m_points.size() ; ++t )
		boundingBox.add(m_points[t]);
	m_root->setBoundingBox(boundingBox);
}

void APSSOctree::estimateAreas(std::vector<float> & areas)
{
	// traverse the octree, consider only the leaves, and for each leaf:
	//    compute an area for the leaf, and split the area among the point samples inside it
	areas.resize( m_points.size() );
	m_root->estimateAreasRecursive(areas);
}

void APSSOctree::buildapssNodeStats( std::vector< glm::vec3 > const & normals , std::vector< float > const & areas )
{
	assert( m_points.size() == normals.size() );

	m_pointsapssNodeStats.resize( m_points.size() );
	for( unsigned int pIt = 0 ; pIt < m_points.size() ; ++pIt ) {
		m_pointsapssNodeStats[pIt] = apssNodeStats(m_points[pIt] , normals[pIt] , areas[pIt]);
	}
	m_numberOfLeaves = m_points.size();
	m_root->buildapssNodeStatsRecursive( m_pointsapssNodeStats );
}

void APSSOctree::index_nodes_by_breadth_first_search()
{
	unsigned int shared_index = 0;
	std::queue< octreeNode * > nodesToVisit;
	nodesToVisit.push( m_root );
	while ( ! nodesToVisit.empty() ) {
		octreeNode * n = nodesToVisit.front();
		n->breadth_first_index = shared_index;
		++shared_index;

		for( unsigned int c = 0 ; c < 8 ; ++c ) {
			octreeNode * child = n->children[c];
			if( child != NULL ) {
				nodesToVisit.push(child);
			}
		}
		nodesToVisit.pop();
	}
}

void APSSOctree::projectCPU(glm::vec3 const & qStart , glm::vec3 & outputPoint , glm::vec3 & outputNormal , unsigned int n_iterations , bool fast,
			  Kernel const& kernel, float stepMaxSize) const
{
	if (!m_knnMode && !m_ballMode)
	{
		if (stepMaxSize < 0)
			stepMaxSize = float(m_root->getBoundingBox().radius());
		glm::vec3 q = qStart;
		for( unsigned int i = 0 ; i < n_iterations ; ++i ) {
			oneStepProjection(q , outputPoint , outputNormal , fast, kernel , m_minDepth , m_scalingProtectionSphere , stepMaxSize );
			q = outputPoint;
		}
	}
	else
		projectAPSSCPU(qStart, outputPoint, outputNormal, n_iterations, kernel, stepMaxSize);
}

void APSSOctree::cpuSphere(glm::vec3 const & qStart, float & u0, glm::vec3 & u123, float & u4, Kernel const & kernel) const
{
	unsigned k = m_knn; //20 nearest neighbors (default), or m_knn
	// find k nearest neighbour
	kdTreeNN::searchResults knn;
	m_kdTree.knearest( qStart , k , knn , false );

	apssStats treeStats;

	for(unsigned int i=0;i<k;i++){
		const glm::vec3 &p = m_points[knn[i].idx];
		const glm::vec3 &n = m_normals[knn[i].idx];

		float wi = kernel.w(p, qStart);
		treeStats.s_ai_wi_pi_ni += wi * glm::dot(p, n);
		treeStats.s_ai_wi_pi_pi += wi * glm::dot(p, p);
		treeStats.s_ai_wi_pi += wi * p;
		treeStats.s_ai_wi_ni += wi * n;
		treeStats.s_ai_wi += wi ;
	}
	u4 = 0.5*( treeStats.s_ai_wi_pi_ni/*/treeStats.s_ai_wi*/ - glm::dot(treeStats.s_ai_wi_pi/*/treeStats.s_ai_wi*/,treeStats.s_ai_wi_ni/treeStats.s_ai_wi) ) /
			( treeStats.s_ai_wi_pi_pi/*/treeStats.s_ai_wi*/ - glm::dot(treeStats.s_ai_wi_pi/*/treeStats.s_ai_wi*/,treeStats.s_ai_wi_pi/treeStats.s_ai_wi) );
	u123= (treeStats.s_ai_wi_ni - 2*u4*treeStats.s_ai_wi_pi)/treeStats.s_ai_wi;
	u0= -(glm::dot(treeStats.s_ai_wi_pi,u123) + u4*treeStats.s_ai_wi_pi_pi)/treeStats.s_ai_wi;

}

void APSSOctree::getKnn(glm::vec3 const & x, unsigned k, std::vector<glm::vec3> & neighbors) const
{
	kdTreeNN::searchResults knn;
	m_kdTree.knearest( x , k , knn , false );
	neighbors.resize(k);
	for (unsigned i=0; i<k; i++)
		neighbors[i] = m_points[knn[i].idx];
}

void APSSOctree::knnMode(bool setMode, unsigned knn)
{
	m_knnMode = setMode;
	if (m_knnMode)
		m_ballMode = !setMode;
	m_knn = knn;
}

void APSSOctree::ballMode(bool setMode, float ballRadius)
{
	m_ballMode = setMode;
	if (m_ballMode)
		m_knnMode = !setMode;
	m_ballRadius = ballRadius;
}

void APSSOctree::projectAPSSCPU(glm::vec3 const & qStart , glm::vec3 & outputPoint , glm::vec3 & outputNormal , unsigned int n_iterations,
			  Kernel const& kernel, float stepMaxSize) const
{
	unsigned k = m_knn; //20 nearest neighbors (default), or m_knn
	outputPoint = qStart;
	for(unsigned int it=0;it<n_iterations;it++)
	{
		unsigned nb_elements = 0;
		kdTreeNN::searchResults knn;
		if (m_knnMode)// find k nearest neighbour
		{
			m_kdTree.knearest( outputPoint , k , knn , false );
			nb_elements = k;
		}
		else// find point in ball of radius m_ballRadius
		{
			if (m_ballRadius < 0.f)
				m_kdTree.getInBall( outputPoint, kernel.sigmas[0] * 3, knn, false);
			else
				m_kdTree.getInBall( outputPoint, m_ballRadius, knn, false);
			nb_elements = knn.size();
			if (nb_elements == 0)//No points inside the ball
			{
				outputPoint = qStart;
				outputNormal = glm::vec3(0, 0, 0);
				return;
			}
		}

		apssStats treeStats;

		for(unsigned int i=0;i<nb_elements;i++){
			const glm::vec3 &p = m_points[knn[i].idx];
			const glm::vec3 &n = m_normals[knn[i].idx];

			float wi = kernel.w(p, outputPoint);
			treeStats.s_ai_wi_pi_ni += wi * glm::dot(p, n);
			treeStats.s_ai_wi_pi_pi += wi * glm::dot(p, p);
			treeStats.s_ai_wi_pi += wi * p;
			treeStats.s_ai_wi_ni += wi * n;
			treeStats.s_ai_wi += wi ;
		}

		if (stepMaxSize < 0)
			stepMaxSize = float(m_root->getBoundingBox().radius());
		projectOnAPSSStats(outputPoint, outputPoint, outputNormal, treeStats, stepMaxSize);
	}
}

void APSSOctree::projectOnAPSSStats( glm::vec3 const & q ,  glm::vec3 & outputPoint , glm::vec3 & outputNormal , apssStats const & treeStats , float step_max_size ) const
{
	// algebraic sphere: u4.||X||^2 + u123.X + u0 = 0_ai
	// geometric sphere: ||X-C||^2 - r^2 = 0
	// geometric plane:  (X-C).n = 0
	float u4 = 0.5*( treeStats.s_ai_wi_pi_ni/*/treeStats.s_ai_wi*/ - glm::dot(treeStats.s_ai_wi_pi/*/treeStats.s_ai_wi*/,treeStats.s_ai_wi_ni/treeStats.s_ai_wi) ) /
			( treeStats.s_ai_wi_pi_pi/*/treeStats.s_ai_wi*/ - glm::dot(treeStats.s_ai_wi_pi/*/treeStats.s_ai_wi*/,treeStats.s_ai_wi_pi/treeStats.s_ai_wi) );
	glm::vec3 u123= (treeStats.s_ai_wi_ni - 2*u4*treeStats.s_ai_wi_pi)/treeStats.s_ai_wi;
	float u0= -(glm::dot(treeStats.s_ai_wi_pi,u123) + u4*treeStats.s_ai_wi_pi_pi)/treeStats.s_ai_wi;

	outputPoint = q;

	if( fabs(u4) < 0.000000000001 ) {
		// then project on a plane (it's a degenerate sphere)
		glm::vec3 n = -u123;
		float lambda = ( u0 - glm::dot( outputPoint , n ) ) / glm::length2(n);
		outputPoint += lambda * n;
		outputNormal = glm::normalize(treeStats.s_ai_wi_ni);
	}
	else {
		glm::vec3 sphere_center = u123/(-2*u4);
		float sphere_radius = sqrt( std::max<float>(0.0 , glm::length2(sphere_center) - u0/u4 ) );

		// projection of the inputpoint onto the sphere:
		glm::vec3 pc= glm::normalize(outputPoint-sphere_center);
		outputPoint = sphere_center + sphere_radius*pc;

		// compute normal by looking at the gradient there:
		outputNormal = glm::normalize(u123 + 2*u4*outputPoint);
	}

	//Invalid points have a normal equals to (0, 0, 0) and their position is kept for debug
	if (std::isnan(outputPoint.x) || std::isnan(outputPoint.y) || std::isnan(outputPoint.z))
	{
		outputPoint = q;
		outputNormal = glm::vec3(0, 0, 0);
		return;
	}

	float step_size = glm::length(outputPoint - q);
	if( step_size > step_max_size ) {
		outputPoint = q + (step_max_size / step_size) * (outputPoint - q) ;
	}
}

void APSSOctree::oneStepProjection(glm::vec3 const & q , glm::vec3 & outputPoint , glm::vec3 & outputNormal , bool fast,
						Kernel const & kernel , unsigned int minimal_depth , float scalingProtectionSphere , float step_max_size ) const
{
	apssStats treeStats;
	if (fast)
		treeStats = m_root->approximateAPSSStatsRecursive(q , kernel , minimal_depth , scalingProtectionSphere , m_pointsapssNodeStats);
	else
		treeStats = m_root->APSSStats(m_numberOfLeaves, q , kernel, m_pointsapssNodeStats);
	projectOnAPSSStats(q, outputPoint , outputNormal, treeStats,step_max_size);
}
