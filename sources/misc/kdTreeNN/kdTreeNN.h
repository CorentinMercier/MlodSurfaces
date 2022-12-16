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

#ifndef kdTreeNN_H
#define kdTreeNN_H


#include <vector>
#include <set>
#include <cassert>
#include <algorithm>
#include <climits>

#define ALLOW_BREADTHFIRST_INDEXING_IN_KDTREENN

#ifdef ALLOW_BREADTHFIRST_INDEXING_IN_KDTREENN
#include <queue>
#endif


namespace kdTreeNN
{

template< class point_t >
class kdTreeAABB
{
public:
	point_t bb,BB;

	kdTreeAABB() { clear(); }

	void clear()
	{
		bb[0] = FLT_MAX;
		BB[0] = -FLT_MAX;
	}
	bool isCleared() const { return bb[0] <= BB[0]; }

	void set( const point_t & p )
	{
		bb = point_t(p[0],p[1],p[2]);
		BB = point_t(p[0],p[1],p[2]);
	}

	void set( const point_t & pbb, const point_t & pBB )
	{
		bb = point_t(pbb[0],pbb[1],pbb[2]);
		BB = point_t(pBB[0],pBB[1],pBB[2]);
	}

	void add( const point_t & p )
	{
		for( unsigned int c = 0 ; c < 3 ; ++c ) {
			bb[c] = std::min( bb[c],p[c] );
			BB[c] = std::max( BB[c],p[c] );
		}
	}

	void add( const kdTreeAABB<point_t> & b )
	{
		for( unsigned int c = 0 ; c < 3 ; ++c ) {
			bb[c] = std::min( bb[c],b.bb[c] );
			BB[c] = std::max( BB[c],b.BB[c] );
		}
	}

	double squareDiagonal() const { return (BB - bb).sqrnorm(); }
	inline double diagonal() const { return sqrt( (BB - bb).sqrnorm() ); }
	inline double radius() const { return diagonal() / 2.0; }
	inline double squareRadius() const { return squareDiagonal() / 4.0; }

	inline char getLargestExtent() const {
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
	inline double getExtentValue(char i) const { return BB[i] - bb[i]; }
	inline double getLargestExtentValue() const { return getExtentValue( getLargestExtent() ); }

	inline double getPseudoExtentInDirection( point_t const & dir ) const {
		return getExtentValue(0) * fabs(dir[0]) +getExtentValue(1) * fabs(dir[1]) +getExtentValue(2) * fabs(dir[2]);
	}

	inline void splitAlongAxis( char axis , double value , kdTreeAABB< point_t > & bbox1 , kdTreeAABB< point_t > & bbox2 ) {
		// for safety:
		value = std::max< double >( std::min< double >( value , BB[axis] ) , bb[axis] );
		point_t BB1 = BB;
		BB1[axis] = value;
		point_t bb2 = bb;
		bb2[axis] = value;
		bbox1.set( bb , BB1 );
		bbox2.set( bb2 , BB );
	}
};




struct SplittingPlane {
	char axis; // 0: along x , 1: along y , 2: along z
	double value;

	SplittingPlane() : axis(0){}
};




struct searchResult {
	unsigned int idx;
	double sqrDist;

	searchResult(unsigned int i , double sd) : idx(i) , sqrDist(sd) {}

	bool operator < (searchResult const & o) const { return sqrDist > o.sqrDist; } // an element is worse if it is further away from the query point
	bool operator > (searchResult const & o) const { return sqrDist < o.sqrDist; }
};




struct searchResults {
	std::vector< searchResult > points;
	double sqrDistMin , sqrDistMax;

	searchResults() : sqrDistMin(FLT_MAX) , sqrDistMax(-FLT_MAX) {}
	void sort() {
		std::sort( points.begin() , points.end() , std::greater< searchResult >() );
		assert( points.front().sqrDist <= points.back().sqrDist );
	}
	void clear() {
		sqrDistMin = FLT_MAX;
		sqrDistMax = -FLT_MAX;
		points.clear();
	}

	bool isCleared() const { return sqrDistMin > sqrDistMax  &&  points.size() == 0; }

	unsigned int size() const { return points.size(); }
	searchResult & operator [] (unsigned int p) { return points[p]; }
	searchResult const & operator [] (unsigned int p) const { return points[p]; }
};




template< class point_t >
class kdTreeNode {
	kdTreeAABB<point_t> boundingBox;
	std::vector< unsigned int > indices; // only the leaves have points indices stored
	unsigned int depth;

	SplittingPlane splittingPlane;
	kdTreeNode * leftChild , * rightChild;

#ifdef ALLOW_BREADTHFIRST_INDEXING_IN_KDTREENN
	unsigned int idx;
#endif

public:
	kdTreeNode() : leftChild(NULL) , rightChild(NULL) {
#ifdef ALLOW_BREADTHFIRST_INDEXING_IN_KDTREENN
		idx = UINT_MAX;
#endif
	}

#ifdef ALLOW_BREADTHFIRST_INDEXING_IN_KDTREENN
	void setIdx( unsigned int i ) { idx = i; }
#endif

	~kdTreeNode() {
		clear();
	}
	kdTreeAABB< point_t > getBoundingBox() const { return boundingBox; }
	SplittingPlane getSplittingPlane() const { return splittingPlane; }
	kdTreeNode * getLeftChild() const { return leftChild; }
	kdTreeNode * getRightChild() const { return rightChild; }

	void setBoundingBox( const kdTreeAABB< point_t > & bbox ) { boundingBox = bbox; }
	void clear() {
		indices.clear();
		if( leftChild != NULL ) {
			delete leftChild;
			leftChild = NULL;
		}
		if( rightChild != NULL ) {
			delete rightChild;
			rightChild = NULL;
		}
	}

	bool isLeaf() const { return leftChild==NULL  &&  rightChild==NULL; }

	void buildRecursive(
			const std::vector< point_t > & i_points,
			const std::vector< unsigned int > & i_indices,
			const SplittingPlane & parentSplittingPlane,
			unsigned int currentDepth,
			unsigned int maxDepth,
			unsigned int maxNumberOfPointsPerLeaf) {
		depth = currentDepth;

		if( i_indices.size() <= maxNumberOfPointsPerLeaf   ||   currentDepth == maxDepth ) {
			indices = i_indices;
			return;
		}

		leftChild = new kdTreeNode;
		rightChild = new kdTreeNode;

		// Decide of a splitting plane, use the round robin rule for now (alternate between x, y, and z):
		splittingPlane.axis = (parentSplittingPlane.axis + 1)%3;
		// for now, go through the median:

		std::vector< double > pointsAlongChosenAxis;
		for( unsigned int p_it = 0 ; p_it < i_indices.size() ; ++p_it ) {
			unsigned int p_idx = i_indices[p_it];
			pointsAlongChosenAxis.push_back(  i_points[p_idx][splittingPlane.axis]  );
		}
		std::sort( pointsAlongChosenAxis.begin() , pointsAlongChosenAxis.end() );
		splittingPlane.value = pointsAlongChosenAxis[ pointsAlongChosenAxis.size() / 2 ];

		std::vector< unsigned int >  indicesLeft;
		std::vector< unsigned int >  indicesRight;
		for( unsigned int p_it = 0 ; p_it < i_indices.size() ; ++p_it ) {
			unsigned int p_idx = i_indices[p_it];
			double x = i_points[p_idx][splittingPlane.axis] ;

			if( x <= splittingPlane.value )
				indicesLeft.push_back(p_idx);
			else
				indicesRight.push_back(p_idx);
		}

		kdTreeAABB<point_t> bbLeft , bbRight;
		boundingBox.splitAlongAxis( splittingPlane.axis , splittingPlane.value , bbLeft , bbRight );
		leftChild->setBoundingBox(bbLeft);
		rightChild->setBoundingBox(bbRight);

		leftChild->buildRecursive(i_points , indicesLeft , splittingPlane , currentDepth + 1 , maxDepth , maxNumberOfPointsPerLeaf );
		rightChild->buildRecursive(i_points , indicesRight , splittingPlane , currentDepth + 1 , maxDepth , maxNumberOfPointsPerLeaf );
	}

	void addPoint( std::vector< point_t > const & i_points ,
				   unsigned int newPointIndex ,
				   const SplittingPlane & parentSplittingPlane,
				   unsigned int currentDepth ,
				   unsigned int maxDepth ,
				   unsigned int maxNumberOfPointsPerLeaf ) {
		// enlarge bounding box if necessary:
		boundingBox.add(i_points[newPointIndex]);

		if( isLeaf() ) {
			if( currentDepth == maxDepth  ||  indices.size() < maxNumberOfPointsPerLeaf ) {
				// then this can stay a leaf, just push the point:
				indices.push_back(newPointIndex);
				return;
			}
			else {
				// then this node should not stay a leaf anymore:
				// create two children, decide of a splitting plane, and split the points in two sets:
				leftChild = new kdTreeNode;
				rightChild = new kdTreeNode;

				// Decide of a splitting plane, use the round robin rule for now (alternate between x, y, and z):
				splittingPlane.axis = (parentSplittingPlane.axis + 1)%3;

				std::vector< double > pointsAlongChosenAxis;
				for( unsigned int t_it = 0 ; t_it < indices.size() ; ++t_it ) {
					unsigned int t = indices[t_it];
					pointsAlongChosenAxis.push_back(  i_points[t][splittingPlane.axis]  );
				}
				std::sort( pointsAlongChosenAxis.begin() , pointsAlongChosenAxis.end() );
				splittingPlane.value = pointsAlongChosenAxis[ pointsAlongChosenAxis.size() / 2 ];

				std::vector< unsigned int >  indicesLeft;
				std::vector< unsigned int >  indicesRight;
				for( unsigned int t_it = 0 ; t_it < indices.size() ; ++t_it ) {
					unsigned int t = indices[t_it];
					double x = i_points[t][splittingPlane.axis] ;

					if( x <= splittingPlane.value )
						indicesLeft.push_back(t);
					else
						indicesRight.push_back(t);
				}

				kdTreeAABB<point_t> bbLeft , bbRight;
				boundingBox.splitAlongAxis( splittingPlane.axis , splittingPlane.value , bbLeft , bbRight );
				leftChild->setBoundingBox(bbLeft);
				rightChild->setBoundingBox(bbRight);

				leftChild->buildRecursive(i_points , indicesLeft , splittingPlane , currentDepth + 1 , maxDepth , maxNumberOfPointsPerLeaf );
				rightChild->buildRecursive(i_points , indicesRight , splittingPlane , currentDepth + 1 , maxDepth , maxNumberOfPointsPerLeaf );

				indices.clear(); // we keep track of indices in the leaves only
				return;
			}
		}
		else {
			// then, just add the point in the correct child:
			double x = i_points[newPointIndex][splittingPlane.axis] ;

			if( x <= splittingPlane.value )
				leftChild->addPoint(i_points,newPointIndex,splittingPlane,currentDepth+1,maxDepth,maxNumberOfPointsPerLeaf);
			else
				rightChild->addPoint(i_points,newPointIndex,splittingPlane,currentDepth+1,maxDepth,maxNumberOfPointsPerLeaf);
		}
	}


	void nearest( point_t const & q ,
				  std::vector< point_t > const & i_points ,
				  searchResult & bestCandidate ) const {
		if(isLeaf()) {
			for( unsigned int lp_it = 0 ; lp_it < indices.size() ; ++lp_it ) {
				unsigned int p_index = indices[lp_it];
				point_t const & p = i_points[p_index];
				double sqrDist = (p[0]-q[0])*(p[0]-q[0]) + (p[1]-q[1])*(p[1]-q[1]) + (p[2]-q[2])*(p[2]-q[2]) ;
				if( sqrDist < bestCandidate.sqrDist ) {
					bestCandidate.idx = p_index;
					bestCandidate.sqrDist = sqrDist;
				}
			}
		}
		else {
			// traverse children in promising order, if the ball search intersects them:
			double x = q[splittingPlane.axis] ;

			if( x <= splittingPlane.value ) {
				leftChild->nearest(q,i_points,bestCandidate);
				if( (x-splittingPlane.value)*(x-splittingPlane.value)  <=  bestCandidate.sqrDist ) {
					rightChild->nearest(q,i_points,bestCandidate);
				}
			}
			else {
				rightChild->nearest(q,i_points,bestCandidate);
				if( (x-splittingPlane.value)*(x-splittingPlane.value)  <=  bestCandidate.sqrDist ) {
					leftChild->nearest(q,i_points,bestCandidate);
				}
			}
		}
	}



#ifdef ALLOW_BREADTHFIRST_INDEXING_IN_KDTREENN
	void gatherNodesForKnearest( point_t const & q ,
								 unsigned int k ,
								 std::set< unsigned int > & visitedNodes ,
								 searchResults & res ,
								 std::vector< point_t > const & i_points  ) const {
		visitedNodes.insert( idx );
		if(isLeaf()) {
			for( unsigned int lp_it = 0 ; lp_it < indices.size() ; ++lp_it ) {
				unsigned int p_index = indices[lp_it];
				point_t const & p = i_points[p_index];
				double sqrDist = (p[0]-q[0])*(p[0]-q[0]) + (p[1]-q[1])*(p[1]-q[1]) + (p[2]-q[2])*(p[2]-q[2]) ;
				if( res.points.size() < k ) {
					res.points.push_back(searchResult(p_index,sqrDist));
					std::push_heap( res.points.begin(),res.points.end() , std::greater<kdTreeNN::searchResult>() );

					res.sqrDistMax = res.points.front().sqrDist;
					res.sqrDistMin = std::min< double >( res.sqrDistMin , sqrDist );
				}
				else if( sqrDist < res.sqrDistMax ) {
					std::pop_heap( res.points.begin(),res.points.end() , std::greater<kdTreeNN::searchResult>() ); // put the worst at the back
					res.points.pop_back(); // get rid of the worst element

					res.points.push_back(searchResult(p_index,sqrDist));
					std::push_heap(res.points.begin(),res.points.end() , std::greater<kdTreeNN::searchResult>()); // put the worst in the front

					res.sqrDistMax = res.points.front().sqrDist;
					res.sqrDistMin = std::min< double >( res.sqrDistMin , sqrDist );
				}
			}
		}
		else {
			// traverse children in promising order, if the ball search intersects them:
			double x = q[splittingPlane.axis] ;

			if( x <= splittingPlane.value ) {
				leftChild->gatherNodesForKnearest(q , k , visitedNodes , res , i_points);
				if( res.points.size() < k || (x-splittingPlane.value)*(x-splittingPlane.value)  <=  res.sqrDistMax ) {
					rightChild->gatherNodesForKnearest(q , k , visitedNodes , res , i_points);
				}
			}
			else {
				rightChild->gatherNodesForKnearest(q , k , visitedNodes , res , i_points);
				if( res.points.size() < k || (x-splittingPlane.value)*(x-splittingPlane.value)  <=  res.sqrDistMax ) {
					leftChild->gatherNodesForKnearest(q , k , visitedNodes , res , i_points);
				}
			}
		}
	}
#endif


	void knearest( point_t const & q ,
				   unsigned int k ,
				   searchResults & res ,
				   std::vector< point_t > const & i_points  ) const {
		if(isLeaf()) {
			for( unsigned int lp_it = 0 ; lp_it < indices.size() ; ++lp_it ) {
				unsigned int p_index = indices[lp_it];
				point_t const & p = i_points[p_index];
				double sqrDist = (p[0]-q[0])*(p[0]-q[0]) + (p[1]-q[1])*(p[1]-q[1]) + (p[2]-q[2])*(p[2]-q[2]) ;
				if( res.points.size() < k ) {
					res.points.push_back(searchResult(p_index,sqrDist));
					std::push_heap( res.points.begin(),res.points.end() , std::greater<kdTreeNN::searchResult>() );

					res.sqrDistMax = res.points.front().sqrDist;
					res.sqrDistMin = std::min< double >( res.sqrDistMin , sqrDist );
				}
				else if( sqrDist < res.sqrDistMax ) {
					std::pop_heap( res.points.begin(),res.points.end() , std::greater<kdTreeNN::searchResult>() ); // put the worst at the back
					res.points.pop_back(); // get rid of the worst element

					res.points.push_back(searchResult(p_index,sqrDist));
					std::push_heap(res.points.begin(),res.points.end() , std::greater<kdTreeNN::searchResult>()); // put the worst in the front

					res.sqrDistMax = res.points.front().sqrDist;
					res.sqrDistMin = std::min< double >( res.sqrDistMin , sqrDist );
				}
			}
		}
		else {
			// traverse children in promising order, if the ball search intersects them:
			double x = q[splittingPlane.axis] ;

			if( x <= splittingPlane.value ) {
				leftChild->knearest(q , k , res , i_points);
				if( res.points.size() < k || (x-splittingPlane.value)*(x-splittingPlane.value)  <=  res.sqrDistMax ) {
					rightChild->knearest(q , k , res , i_points);
				}
			}
			else {
				rightChild->knearest(q , k , res , i_points);
				if( res.points.size() < k || (x-splittingPlane.value)*(x-splittingPlane.value)  <=  res.sqrDistMax ) {
					leftChild->knearest(q , k , res , i_points);
				}
			}
		}
	}

	void getInBall( point_t const & q , double r , searchResults & res , std::vector< point_t > const & i_points ) const {
		if(isLeaf()) {
			for( unsigned int lp_it = 0 ; lp_it < indices.size() ; ++lp_it ) {
				unsigned int p_index = indices[lp_it];
				point_t const & p = i_points[p_index];
				double sqrDist = (p[0]-q[0])*(p[0]-q[0]) + (p[1]-q[1])*(p[1]-q[1]) + (p[2]-q[2])*(p[2]-q[2]) ;
				if( sqrDist <= r*r ) {
					res.points.push_back( searchResult(p_index,sqrDist) );
					res.sqrDistMin = std::min< double >( res.sqrDistMin , sqrDist );
					res.sqrDistMax = std::max< double >( res.sqrDistMax , sqrDist );
				}
			}
		}
		else {
			// traverse children, if the ball search intersects them:
			double x = q[splittingPlane.axis] ;

			if( x <= splittingPlane.value ) {
				leftChild->getInBall(q , r , res , i_points);
				if( x + r >= splittingPlane.value ) {
					rightChild->getInBall(q , r , res , i_points);
				}
			}
			else {
				rightChild->getInBall(q , r , res , i_points);
				if( x - r <= splittingPlane.value ) {
					leftChild->getInBall(q , r , res , i_points);
				}
			}
		}
	}
	bool ballIsEmpty( point_t const & q , double r , std::vector< point_t > const & i_points ) const {
		if(isLeaf()) {
			for( unsigned int lp_it = 0 ; lp_it < indices.size() ; ++lp_it ) {
				unsigned int p_index = indices[lp_it];
				point_t const & p = i_points[p_index];
				double sqrDist = (p[0]-q[0])*(p[0]-q[0]) + (p[1]-q[1])*(p[1]-q[1]) + (p[2]-q[2])*(p[2]-q[2]) ;
				if( sqrDist <= r*r ) {
					return false;
				}
			}
			return true;
		}
		else {
			// traverse children, if the ball search intersects them:
			double x = q[splittingPlane.axis] ;

			if( x <= splittingPlane.value ) {
				if( ! leftChild->ballIsEmpty(q , r , i_points) )
					return false;
				if( x + r >= splittingPlane.value ) {
					if( ! rightChild->ballIsEmpty(q , r , i_points) )
						return false;
				}
			}
			else {
				if( ! rightChild->ballIsEmpty(q , r , i_points) )
					return false;
				if( x - r <= splittingPlane.value ) {
					if( ! leftChild->ballIsEmpty(q , r , i_points) )
						return false;
				}
			}
		}
		return true;
	}
};






template< class point_t >
class kdTree{
public:
	kdTree()  {}
	~kdTree()  { clear(); }
	void clear() { root.clear(); }
	std::vector< point_t > const & getPoints() const { return points; }

	// Without a mesh access structure:
	void build( const std::vector< point_t > & i_points ,
				bool updateRootBoundingBox = true , // sometimes you don t have all the points when you build, but you want to setup a fixed bounding box
				unsigned int maxDepth = 12,
				unsigned int maxNumberOfPointsPerLeaf = 1)
	{
		assert( i_points.size() > 0  &&  "Come ooooon" );
		points = i_points;
		if(updateRootBoundingBox) computeBoundingBox();
		std::vector< unsigned int > indices( points.size() );
		for( unsigned int p_idx = 0 ; p_idx < points.size() ; ++p_idx )
			indices[p_idx] = p_idx;
		initSplittingPlane.axis = getBoundingBox().getLargestExtent();
		root.buildRecursive(points,indices , initSplittingPlane , 0 , maxDepth , maxNumberOfPointsPerLeaf);
	}


	kdTreeAABB<point_t> getBoundingBox() const { return root.getBoundingBox(); }

	bool empty() const { return points.size() == 0; }
	unsigned int nPoints() const { return points.size(); }

	void addPoint( point_t const & p ,
				   unsigned int maxDepth = 12,
				   unsigned int maxNumberOfPointsPerLeaf = 1
			) {
		points.push_back(p);
		// add point in the hierarchy:
		root.addPoint( points , points.size() - 1 , initSplittingPlane , 0 , maxDepth , maxNumberOfPointsPerLeaf);
	}

	searchResult nearest( point_t const & q ) const {
		searchResult bestCandidate(UINT_MAX,FLT_MAX);
		root.nearest(q , points , bestCandidate);
		return bestCandidate;
	}

	void getInBall( point_t const & q , double r , searchResults & res , bool sortThem = false ) const {
		root.getInBall(q , r , res , points);
		if(sortThem) {
			res.sort();
		}
	}
	bool ballIsEmpty( point_t const & q , double r ) const {
		return root.ballIsEmpty(q , r , points);
	}
	void knearest( point_t const & q , unsigned int k , searchResults & res , bool sortThem = false ) const {
		assert( res.isCleared() );
		std::make_heap(res.points.begin() , res.points.end() , std::greater<kdTreeNN::searchResult>() );
		res.sqrDistMax = FLT_MAX;
		root.knearest(q , k , res , points);
		if( res.points.size() != std::min<unsigned int>(k , nPoints()) ) {
			std::cout << "Error : requested " << k << " NNs, got " << res.points.size() << std::endl;
			assert( res.points.size() == std::min<unsigned int>(k , nPoints()) );
		}
		if(sortThem) {
			std::sort_heap( res.points.begin(),res.points.end() , std::greater<kdTreeNN::searchResult>() );
		}
	}



#ifdef ALLOW_BREADTHFIRST_INDEXING_IN_KDTREENN
	void index_nodes_by_breadth_first_search() {
		unsigned int shared_index = 0;
		std::queue< kdTreeNode< point_t > * > nodesToVisit;
		nodesToVisit.push( &root );
		while ( ! nodesToVisit.empty() ) {
			kdTreeNode< point_t > * n = nodesToVisit.front();
			n->setIdx( shared_index );
			++shared_index;

			{
				kdTreeNode< point_t > * child = n->getLeftChild();
				if( child != NULL ) {
					nodesToVisit.push(child);
				}
			}
			{
				kdTreeNode< point_t > * child = n->getRightChild();
				if( child != NULL ) {
					nodesToVisit.push(child);
				}
			}

			nodesToVisit.pop();
		}
	}


	void gatherNodesForKnearest( point_t const & q , unsigned int k , std::set< unsigned int > & visitedNodes , searchResults & res , bool sortThem = false ) const {
		assert( res.isCleared() );
		std::make_heap(res.points.begin() , res.points.end() , std::greater<kdTreeNN::searchResult>() );
		res.sqrDistMax = FLT_MAX;
		root.gatherNodesForKnearest(q , k , visitedNodes , res , points);
		if( res.points.size() != std::min<unsigned int>(k , nPoints()) ) {
			std::cout << "Error : requested " << k << " NNs, got " << res.points.size() << std::endl;
			assert( res.points.size() == std::min<unsigned int>(k , nPoints()) );
		}
		if(sortThem) {
			std::sort_heap( res.points.begin(),res.points.end() , std::greater<kdTreeNN::searchResult>() );
		}
	}
#endif

private:
	kdTreeNode< point_t > root;
	std::vector< point_t > points;
	SplittingPlane initSplittingPlane;

	void computeBoundingBox( ) {
		assert( points.size() > 0  &&  "Come ooooon" );
		kdTreeAABB<point_t> boundingBox;
		boundingBox.set(points[0]);
		for( unsigned int t = 1 ; t < points.size() ; ++t )
			boundingBox.add(points[t]);
		root.setBoundingBox(boundingBox);
	}
};



} // namespace kdTreeNN


#endif // kdTreeNN_H
