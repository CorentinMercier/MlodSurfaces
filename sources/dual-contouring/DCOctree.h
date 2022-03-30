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

#include <vector>
#include <queue>
#include <cassert>
#include <algorithm>

#include "src/Vec3.h"

// ------------------------------------------------------------------------------------------------------------
// Octree for Dual Contouring
// ------------------------------------------------------------------------------------------------------------
namespace DCOctree {
	using BBOX = BBOX2;

struct QEM{
	Mat3 A; Vec3 b; double c;
	// quadric(q) = q^T A q - 2 b^T q + c
	// minimizer: q^* = A^{-1}b

	QEM() : A(0,0,0,0,0,0,0,0,0) , b(0,0,0) , c(0) {}
	QEM( Mat3 const _A , Vec3 const & _b , double _c ) : A(_A) , b(_b) , c(_c) {}

	static
	QEM fromPlan( Vec3 const & n , Vec3 const & p ) {
		return QEM( tensor(n,n) , dot( p , n ) * n , dot( p , n )*dot( p , n ) );
	}

	void operator += (QEM const & o) {
		A += o.A;
		b += o.b;
		c += o.c;
	}

	Vec3 minimizer(Vec3 const & pInit) {
		// if det(A) != 0 then returns A^{-1}b, otherwise minimize quadric(q) for q = pInit + s * b
		const float det = determinant(A);
		if(fabs(det) >= 0.000001f)
			return inverse(A) * b;
		else {
			return pInit + ( ( length2(b) - dot(pInit , A*b)) / (dot(b , A*b)) ) * b;
		}
	}
};


// An octree is "regular" if a node either is a leaf, or has exactly 8 children (either all or none children exist)
class RegularOctree{
public:
	struct RegularOctreeNode {
		RegularOctreeNode * children[8]; // these are organized as zyx (bit-wise)
		bool is_leaf;

		unsigned int depth;
		BBOX boundingBox;
		void setBoundingBox( const BBOX & bbox ) { boundingBox = bbox; }
		BBOX const & getBoundingBox() const { return boundingBox; }

		unsigned int idx; // unique identifier of the node, can come in handy
		unsigned int idx_point_in_mesh; // used for the dual contouring

		QEM qem;

		RegularOctreeNode() {
			idx_point_in_mesh = UINT_MAX;
			for(unsigned int c = 0 ; c < 8 ; ++c) children[c] = NULL;
			is_leaf = true;
			idx = UINT_MAX;
		}

		bool isLeaf() const { return is_leaf ; }

		void subdivideIfLeaf( ) {
			if( is_leaf ) {
				is_leaf = false;

				double xSplit = (boundingBox.bmin[0] + boundingBox.bmax[0])/2.0;
				double ySplit = (boundingBox.bmin[1] + boundingBox.bmax[1])/2.0;
				double zSplit = (boundingBox.bmin[2] + boundingBox.bmax[2])/2.0;

				double xVals [3] = {boundingBox.bmin[0] , xSplit , boundingBox.bmax[0]};
				double yVals [3] = {boundingBox.bmin[1] , ySplit , boundingBox.bmax[1]};
				double zVals [3] = {boundingBox.bmin[2] , zSplit , boundingBox.bmax[2]};

				for( unsigned int xi = 0 ; xi < 2 ; ++xi ) {
					for( unsigned int yi = 0 ; yi < 2 ; ++yi ) {
						for( unsigned int zi = 0 ; zi < 2 ; ++zi ) {
							unsigned int c = xi + 2 * yi + 4 * zi;
							children[c] = new RegularOctreeNode;
							BBOX bbox;
							bbox.set( Vec3( xVals[xi] , yVals[yi] , zVals[zi] )  ,  Vec3( xVals[xi+1] , yVals[yi+1] , zVals[zi+1] ) );
							children[c]->setBoundingBox(bbox);
							children[c]->depth = depth + 1;
							children[c]->is_leaf = true;
						}
					}
				}
			}
		}

		void buildRecursive(
				const std::vector< Vec3 > & i_points,
				const std::vector< unsigned int > & i_indices,
				unsigned int maxDepth,
				unsigned int maxNumberOfPointsPerLeaf) {
			if( i_indices.size() <= maxNumberOfPointsPerLeaf   ||   depth == maxDepth ) {
				is_leaf = true;
				return;
			}

			// otherwise:
			subdivideIfLeaf();

			std::vector< unsigned int > childrenIndices [8];
			double xSplit = (boundingBox.bmin[0] + boundingBox.bmax[0])/2.0;
			double ySplit = (boundingBox.bmin[1] + boundingBox.bmax[1])/2.0;
			double zSplit = (boundingBox.bmin[2] + boundingBox.bmax[2])/2.0;


			for( unsigned int pIt = 0 ; pIt < i_indices.size() ; ++pIt ) {
				unsigned int pIdx = i_indices[pIt];
				Vec3 const & p = i_points[pIdx];

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

			for( unsigned int xi = 0 ; xi < 2 ; ++xi ) {
				for( unsigned int yi = 0 ; yi < 2 ; ++yi ) {
					for( unsigned int zi = 0 ; zi < 2 ; ++zi ) {
						unsigned int c = xi + 2 * yi + 4 * zi;
						if( childrenIndices[c].size() > 0 ) {
							children[c]->buildRecursive(i_points , childrenIndices[c] , maxDepth , maxNumberOfPointsPerLeaf );
						}
					}
				}
			}
		}

		void enforceMinimalDepth( unsigned int min_depth_to_reach ) {
			if( depth >= min_depth_to_reach ) {
				return;
			}

			// otherwise, go deeper:
			if( is_leaf ) {
				// then first allocate:
				subdivideIfLeaf();
			}

			for( unsigned int xi = 0 ; xi < 2 ; ++xi ) {
				for( unsigned int yi = 0 ; yi < 2 ; ++yi ) {
					for( unsigned int zi = 0 ; zi < 2 ; ++zi ) {
						unsigned int c = xi + 2 * yi + 4 * zi;
						children[c]->enforceMinimalDepth( min_depth_to_reach );
					}
				}
			}
		}


		//--------------------------------------------------------------------------------------------------------------------------//
		//--------------------------------------------------------------------------------------------------------------------------//
		// A good explanation of the minimal edges gathering can be found in
		// Adaptive surface extraction from anisotropic volumetricdata: contouring on generalized octrees
		// Ricardo Uribe Lobello, Florence Denis, Florent Dupont
		// section 6
		//--------------------------------------------------------------------------------------------------------------------------//
		//--------------------------------------------------------------------------------------------------------------------------//
		template< class F >
		static
		void EdgeProc_traverse_minimal_edges( RegularOctreeNode * n1 , RegularOctreeNode * n2 , RegularOctreeNode * n3 , RegularOctreeNode * n4 , unsigned char axis , F & f ) {
			// n1  n2
			// n3  n4
			// if all the nodes are leaves, we have a minimal edge
			if( n1->isLeaf()  &&  n2->isLeaf()  &&  n3->isLeaf()  &&  n4->isLeaf() ) {
				// bingo, we have a minimal edge
				f( n1,n2,n3,n4,axis );
			}
			else {
				if( axis == 0 ) {
					{
						RegularOctreeNode * _n1 = n1->isLeaf() ? n1 : n1->children[2];
						RegularOctreeNode * _n2 = n2->isLeaf() ? n2 : n2->children[6];
						RegularOctreeNode * _n3 = n3->isLeaf() ? n3 : n3->children[0];
						RegularOctreeNode * _n4 = n4->isLeaf() ? n4 : n4->children[4];
						EdgeProc_traverse_minimal_edges( _n1 , _n2 , _n3 , _n4 , axis , f );
					}
					{
						RegularOctreeNode * _n1 = n1->isLeaf() ? n1 : n1->children[3];
						RegularOctreeNode * _n2 = n2->isLeaf() ? n2 : n2->children[7];
						RegularOctreeNode * _n3 = n3->isLeaf() ? n3 : n3->children[1];
						RegularOctreeNode * _n4 = n4->isLeaf() ? n4 : n4->children[5];
						EdgeProc_traverse_minimal_edges( _n1 , _n2 , _n3 , _n4 , axis , f );
					}
				}
				else if( axis == 1 ) {
					{
						RegularOctreeNode * _n1 = n1->isLeaf() ? n1 : n1->children[1];
						RegularOctreeNode * _n2 = n2->isLeaf() ? n2 : n2->children[0];
						RegularOctreeNode * _n3 = n3->isLeaf() ? n3 : n3->children[5];
						RegularOctreeNode * _n4 = n4->isLeaf() ? n4 : n4->children[4];
						EdgeProc_traverse_minimal_edges( _n1 , _n2 , _n3 , _n4 , axis , f );
					}
					{
						RegularOctreeNode * _n1 = n1->isLeaf() ? n1 : n1->children[3];
						RegularOctreeNode * _n2 = n2->isLeaf() ? n2 : n2->children[2];
						RegularOctreeNode * _n3 = n3->isLeaf() ? n3 : n3->children[7];
						RegularOctreeNode * _n4 = n4->isLeaf() ? n4 : n4->children[6];
						EdgeProc_traverse_minimal_edges( _n1 , _n2 , _n3 , _n4 , axis , f );
					}
				}
				else if( axis == 2 ) {
					{
						RegularOctreeNode * _n1 = n1->isLeaf() ? n1 : n1->children[3];
						RegularOctreeNode * _n2 = n2->isLeaf() ? n2 : n2->children[2];
						RegularOctreeNode * _n3 = n3->isLeaf() ? n3 : n3->children[1];
						RegularOctreeNode * _n4 = n4->isLeaf() ? n4 : n4->children[0];
						EdgeProc_traverse_minimal_edges( _n1 , _n2 , _n3 , _n4 , axis , f );
					}
					{
						RegularOctreeNode * _n1 = n1->isLeaf() ? n1 : n1->children[7];
						RegularOctreeNode * _n2 = n2->isLeaf() ? n2 : n2->children[6];
						RegularOctreeNode * _n3 = n3->isLeaf() ? n3 : n3->children[5];
						RegularOctreeNode * _n4 = n4->isLeaf() ? n4 : n4->children[4];
						EdgeProc_traverse_minimal_edges( _n1 , _n2 , _n3 , _n4 , axis , f );
					}
				}
			}
		}

		template< class F >
		static
		void FaceProc_traverse_minimal_edges( RegularOctreeNode * n1 , RegularOctreeNode * n2 , unsigned char axis , F & f ) {
			if( n1->isLeaf()  &&  n2->isLeaf() )
				return;

			// call FaceProc
			if( axis == 0 ) {
				{
					RegularOctreeNode * _n1 = n1->isLeaf() ? n1 : n1->children[1];
					RegularOctreeNode * _n2 = n2->isLeaf() ? n2 : n2->children[0];
					FaceProc_traverse_minimal_edges( _n1 , _n2 , axis , f );
				}
				{
					RegularOctreeNode * _n1 = n1->isLeaf() ? n1 : n1->children[5];
					RegularOctreeNode * _n2 = n2->isLeaf() ? n2 : n2->children[4];
					FaceProc_traverse_minimal_edges( _n1 , _n2 , axis , f );
				}
				{
					RegularOctreeNode * _n1 = n1->isLeaf() ? n1 : n1->children[3];
					RegularOctreeNode * _n2 = n2->isLeaf() ? n2 : n2->children[2];
					FaceProc_traverse_minimal_edges( _n1 , _n2 , axis , f );
				}
				{
					RegularOctreeNode * _n1 = n1->isLeaf() ? n1 : n1->children[7];
					RegularOctreeNode * _n2 = n2->isLeaf() ? n2 : n2->children[6];
					FaceProc_traverse_minimal_edges( _n1 , _n2 , axis , f );
				}
			}
			else if( axis == 1 ) {
				{
					RegularOctreeNode * _n1 = n1->isLeaf() ? n1 : n1->children[2];
					RegularOctreeNode * _n2 = n2->isLeaf() ? n2 : n2->children[0];
					FaceProc_traverse_minimal_edges( _n1 , _n2 , axis , f );
				}
				{
					RegularOctreeNode * _n1 = n1->isLeaf() ? n1 : n1->children[7];
					RegularOctreeNode * _n2 = n2->isLeaf() ? n2 : n2->children[5];
					FaceProc_traverse_minimal_edges( _n1 , _n2 , axis , f );
				}
				{
					RegularOctreeNode * _n1 = n1->isLeaf() ? n1 : n1->children[3];
					RegularOctreeNode * _n2 = n2->isLeaf() ? n2 : n2->children[1];
					FaceProc_traverse_minimal_edges( _n1 , _n2 , axis , f );
				}
				{
					RegularOctreeNode * _n1 = n1->isLeaf() ? n1 : n1->children[6];
					RegularOctreeNode * _n2 = n2->isLeaf() ? n2 : n2->children[4];
					FaceProc_traverse_minimal_edges( _n1 , _n2 , axis , f );
				}
			}
			else if( axis == 2 ) {
				{
					RegularOctreeNode * _n1 = n1->isLeaf() ? n1 : n1->children[4];
					RegularOctreeNode * _n2 = n2->isLeaf() ? n2 : n2->children[0];
					FaceProc_traverse_minimal_edges( _n1 , _n2 , axis , f );
				}
				{
					RegularOctreeNode * _n1 = n1->isLeaf() ? n1 : n1->children[5];
					RegularOctreeNode * _n2 = n2->isLeaf() ? n2 : n2->children[1];
					FaceProc_traverse_minimal_edges( _n1 , _n2 , axis , f );
				}
				{
					RegularOctreeNode * _n1 = n1->isLeaf() ? n1 : n1->children[6];
					RegularOctreeNode * _n2 = n2->isLeaf() ? n2 : n2->children[2];
					FaceProc_traverse_minimal_edges( _n1 , _n2 , axis , f );
				}
				{
					RegularOctreeNode * _n1 = n1->isLeaf() ? n1 : n1->children[7];
					RegularOctreeNode * _n2 = n2->isLeaf() ? n2 : n2->children[3];
					FaceProc_traverse_minimal_edges( _n1 , _n2 , axis , f );
				}
			}
			else {
				assert( 0  &&  "invalid axis");
			}

			// call EdgeProc:
			if( axis == 0 ) {
				// first consider edges along the y axis
				{
					// ok
					RegularOctreeNode * _n1 = n1->isLeaf() ? n1 : n1->children[5];
					RegularOctreeNode * _n2 = n2->isLeaf() ? n2 : n2->children[4];
					RegularOctreeNode * _n3 = n1->isLeaf() ? n1 : n1->children[1];
					RegularOctreeNode * _n4 = n2->isLeaf() ? n2 : n2->children[0];
					EdgeProc_traverse_minimal_edges( _n1 , _n2 , _n3 , _n4 , 1 , f );
				}
				{
					// ok
					RegularOctreeNode * _n1 = n1->isLeaf() ? n1 : n1->children[7];
					RegularOctreeNode * _n2 = n2->isLeaf() ? n2 : n2->children[6];
					RegularOctreeNode * _n3 = n1->isLeaf() ? n1 : n1->children[3];
					RegularOctreeNode * _n4 = n2->isLeaf() ? n2 : n2->children[2];
					EdgeProc_traverse_minimal_edges( _n1 , _n2 , _n3 , _n4 , 1 , f );
				}
				// then consider edges along the z axis
				{
					// ok
					RegularOctreeNode * _n1 = n1->isLeaf() ? n1 : n1->children[1];
					RegularOctreeNode * _n2 = n2->isLeaf() ? n2 : n2->children[0];
					RegularOctreeNode * _n3 = n1->isLeaf() ? n1 : n1->children[3];
					RegularOctreeNode * _n4 = n2->isLeaf() ? n2 : n2->children[2];
					EdgeProc_traverse_minimal_edges( _n1 , _n2 , _n3 , _n4 , 2 , f );
				}
				{
					// ok
					RegularOctreeNode * _n1 = n1->isLeaf() ? n1 : n1->children[5];
					RegularOctreeNode * _n2 = n2->isLeaf() ? n2 : n2->children[4];
					RegularOctreeNode * _n3 = n1->isLeaf() ? n1 : n1->children[7];
					RegularOctreeNode * _n4 = n2->isLeaf() ? n2 : n2->children[6];
					EdgeProc_traverse_minimal_edges( _n1 , _n2 , _n3 , _n4 , 2 , f );
				}
			}
			else if( axis == 1 ) {
				// first consider edges along the x axis
				{
					// ok
					RegularOctreeNode * _n1 = n1->isLeaf() ? n1 : n1->children[6];
					RegularOctreeNode * _n2 = n1->isLeaf() ? n1 : n1->children[2];
					RegularOctreeNode * _n3 = n2->isLeaf() ? n2 : n2->children[4];
					RegularOctreeNode * _n4 = n2->isLeaf() ? n2 : n2->children[0];
					EdgeProc_traverse_minimal_edges( _n1 , _n2 , _n3 , _n4 , 0 , f );
				}
				{
					// ok
					RegularOctreeNode * _n1 = n1->isLeaf() ? n1 : n1->children[7];
					RegularOctreeNode * _n2 = n1->isLeaf() ? n1 : n1->children[3];
					RegularOctreeNode * _n3 = n2->isLeaf() ? n2 : n2->children[5];
					RegularOctreeNode * _n4 = n2->isLeaf() ? n2 : n2->children[1];
					EdgeProc_traverse_minimal_edges( _n1 , _n2 , _n3 , _n4 , 0 , f );
				}
				// then consider edges along the z axis
				{
					// ok
					RegularOctreeNode * _n1 = n1->isLeaf() ? n1 : n1->children[2];
					RegularOctreeNode * _n2 = n1->isLeaf() ? n1 : n1->children[3];
					RegularOctreeNode * _n3 = n2->isLeaf() ? n2 : n2->children[0];
					RegularOctreeNode * _n4 = n2->isLeaf() ? n2 : n2->children[1];
					EdgeProc_traverse_minimal_edges( _n1 , _n2 , _n3 , _n4 , 2 , f );
				}
				{
					// ok
					RegularOctreeNode * _n1 = n1->isLeaf() ? n1 : n1->children[6];
					RegularOctreeNode * _n2 = n1->isLeaf() ? n1 : n1->children[7];
					RegularOctreeNode * _n3 = n2->isLeaf() ? n2 : n2->children[4];
					RegularOctreeNode * _n4 = n2->isLeaf() ? n2 : n2->children[5];
					EdgeProc_traverse_minimal_edges( _n1 , _n2 , _n3 , _n4 , 2 , f );
				}
			}
			else if( axis == 2 ) {
				// first consider edges along the x axis
				{
					// ok
					RegularOctreeNode * _n1 = n2->isLeaf() ? n2 : n2->children[0];
					RegularOctreeNode * _n2 = n1->isLeaf() ? n1 : n1->children[4];
					RegularOctreeNode * _n3 = n2->isLeaf() ? n2 : n2->children[2];
					RegularOctreeNode * _n4 = n1->isLeaf() ? n1 : n1->children[6];
					EdgeProc_traverse_minimal_edges( _n1 , _n2 , _n3 , _n4 , 0 , f );
				}
				{
					// ok
					RegularOctreeNode * _n1 = n2->isLeaf() ? n2 : n2->children[1];
					RegularOctreeNode * _n2 = n1->isLeaf() ? n1 : n1->children[5];
					RegularOctreeNode * _n3 = n2->isLeaf() ? n2 : n2->children[3];
					RegularOctreeNode * _n4 = n1->isLeaf() ? n1 : n1->children[7];
					EdgeProc_traverse_minimal_edges( _n1 , _n2 , _n3 , _n4 , 0 , f );
				}
				// then consider edges along the y axis
				{
					// ok
					RegularOctreeNode * _n1 = n2->isLeaf() ? n2 : n2->children[0];
					RegularOctreeNode * _n2 = n2->isLeaf() ? n2 : n2->children[1];
					RegularOctreeNode * _n3 = n1->isLeaf() ? n1 : n1->children[4];
					RegularOctreeNode * _n4 = n1->isLeaf() ? n1 : n1->children[5];
					EdgeProc_traverse_minimal_edges( _n1 , _n2 , _n3 , _n4 , 1 , f );
				}
				{
					// ok
					RegularOctreeNode * _n1 = n2->isLeaf() ? n2 : n2->children[2];
					RegularOctreeNode * _n2 = n2->isLeaf() ? n2 : n2->children[3];
					RegularOctreeNode * _n3 = n1->isLeaf() ? n1 : n1->children[6];
					RegularOctreeNode * _n4 = n1->isLeaf() ? n1 : n1->children[7];
					EdgeProc_traverse_minimal_edges( _n1 , _n2 , _n3 , _n4 , 1 , f );
				}
			}
			else {
				assert( 0  &&  "invalid axis");
			}
		}

		template< class F >
		static
		void CellProc_traverse_minimal_edges( RegularOctreeNode * n , F & f ) {
			if( ! n->isLeaf() ) {
				for( int c = 0 ; c < 8 ; ++c )
					CellProc_traverse_minimal_edges( n->children[c] , f );

				// call FaceProc
				FaceProc_traverse_minimal_edges( n->children[0] , n->children[1] , 0 , f );
				FaceProc_traverse_minimal_edges( n->children[2] , n->children[3] , 0 , f );
				FaceProc_traverse_minimal_edges( n->children[4] , n->children[5] , 0 , f );
				FaceProc_traverse_minimal_edges( n->children[6] , n->children[7] , 0 , f );

				FaceProc_traverse_minimal_edges( n->children[0] , n->children[2] , 1 , f );
				FaceProc_traverse_minimal_edges( n->children[1] , n->children[3] , 1 , f );
				FaceProc_traverse_minimal_edges( n->children[4] , n->children[6] , 1 , f );
				FaceProc_traverse_minimal_edges( n->children[5] , n->children[7] , 1 , f );

				FaceProc_traverse_minimal_edges( n->children[0] , n->children[4] , 2 , f );
				FaceProc_traverse_minimal_edges( n->children[1] , n->children[5] , 2 , f );
				FaceProc_traverse_minimal_edges( n->children[2] , n->children[6] , 2 , f );
				FaceProc_traverse_minimal_edges( n->children[3] , n->children[7] , 2 , f );

				//call EdgeProc
				EdgeProc_traverse_minimal_edges( n->children[4] , n->children[0] , n->children[6] , n->children[2] , 0 , f );
				EdgeProc_traverse_minimal_edges( n->children[5] , n->children[1] , n->children[7] , n->children[3] , 0 , f );

				EdgeProc_traverse_minimal_edges( n->children[4] , n->children[5] , n->children[0] , n->children[1] , 1 , f );
				EdgeProc_traverse_minimal_edges( n->children[6] , n->children[7] , n->children[2] , n->children[3] , 1 , f );

				EdgeProc_traverse_minimal_edges( n->children[0] , n->children[1] , n->children[2] , n->children[3] , 2 , f );
				EdgeProc_traverse_minimal_edges( n->children[4] , n->children[5] , n->children[6] , n->children[7] , 2 , f );
			}
		}
	};

private:
	RegularOctreeNode root;
	std::vector< Vec3 > points;

	void computeBoundingBox( ) {
		assert( points.size() > 0  &&  "Come ooooon" );
		BBOX boundingBox;
		boundingBox.set(points[0]);
		for( unsigned int t = 1 ; t < points.size() ; ++t )
			boundingBox.add(points[t]);
		root.setBoundingBox(boundingBox);
	}




public:
	RegularOctree() {}

	BBOX const & getBoundingBox() const { return root.getBoundingBox(); }

	void setBoundingBox( Vec3 const & bb , Vec3 const & BB ) { root.setBoundingBox(BBOX(bb,BB)); }

	bool empty() const { return points.size() == 0; }
	unsigned int nPoints() const { return points.size(); }


	void index_nodes_by_breadth_first_search() {
		unsigned int shared_index = 0;
		std::queue< RegularOctree::RegularOctreeNode * > nodesToVisit;
		nodesToVisit.push( &root );
		while ( ! nodesToVisit.empty() ) {
			RegularOctree::RegularOctreeNode * n = nodesToVisit.front();
			n->idx = shared_index;
			++shared_index;

			for( unsigned int c = 0 ; c < 8 ; ++c ) {
				RegularOctree::RegularOctreeNode * child = n->children[c];
				if( child != NULL ) {
					nodesToVisit.push(child);
				}
			}

			nodesToVisit.pop();
		}
	}

	// Without a mesh access structure:
	void build( const std::vector< Vec3 > & i_points ,
				bool updateRootBoundingBox = true , // sometimes you don t have all the points when you build, but you want to setup a fixed bounding box
				unsigned int maxDepth = 12,
				unsigned int maxNumberOfPointsPerLeaf = 1) {
		assert( i_points.size() > 0  &&  "Come ooooon" );
		points = i_points;
		if(updateRootBoundingBox) computeBoundingBox();
		std::vector< unsigned int > indices( points.size() );
		for( unsigned int t = 0 ; t < points.size() ; ++t )
			indices[t] = t;

		root.depth = 0;
		root.buildRecursive(points, indices , maxDepth , maxNumberOfPointsPerLeaf);
		index_nodes_by_breadth_first_search();
	}


	void enforceMinimalDepth( unsigned int min_depth_to_reach ) {
		root.enforceMinimalDepth( min_depth_to_reach );
	}

	void grade() {
		assert( 0  &&  "todo : octree::grade()");
	}


	struct Edge {
		RegularOctree::RegularOctreeNode * n1 , * n2 , * n3 ,  * n4;
		unsigned char axis;
		Edge(
				RegularOctree::RegularOctreeNode * _n1 , RegularOctree::RegularOctreeNode * _n2 ,
				RegularOctree::RegularOctreeNode * _n3 , RegularOctree::RegularOctreeNode * _n4,
				unsigned char _axis ) {
			n1 = _n1;
			n2 = _n2;
			n3 = _n3;
			n4 = _n4;
			axis = _axis;
		}

	public:
		RegularOctree::RegularOctreeNode * getN1() { return n1; }
		RegularOctree::RegularOctreeNode * getN2() { return n2; }
		RegularOctree::RegularOctreeNode * getN3() { return n3; }
		RegularOctree::RegularOctreeNode * getN4() { return n4; }
	};

	void locate_minimal_edge( Edge const & min_edge , Vec3 & p1 , Vec3 & p2 ) {
		if( min_edge.axis == 0 ) {
			if( min_edge.n1->depth >= min_edge.n2->depth   &&  min_edge.n1->depth >= min_edge.n3->depth   &&   min_edge.n1->depth >= min_edge.n4->depth ) {
				Vec3 b = min_edge.n1->boundingBox.bmin , B = min_edge.n1->boundingBox.bmax;
				p1 = Vec3( b[0] , B[1] , b[2] );        p2 = Vec3( B[0] , B[1] , b[2] );
			}
			else if( min_edge.n2->depth >= min_edge.n1->depth   &&  min_edge.n2->depth >= min_edge.n3->depth   &&   min_edge.n2->depth >= min_edge.n4->depth ) {
				Vec3 b = min_edge.n2->boundingBox.bmin , B = min_edge.n2->boundingBox.bmax;
				p1 = Vec3( b[0] , B[1] , B[2] );        p2 = Vec3( B[0] , B[1] , B[2] );
			}
			else if( min_edge.n3->depth >= min_edge.n1->depth   &&  min_edge.n3->depth >= min_edge.n2->depth   &&   min_edge.n3->depth >= min_edge.n4->depth ) {
				Vec3 b = min_edge.n3->boundingBox.bmin , B = min_edge.n3->boundingBox.bmax;
				p1 = Vec3( b[0] , b[1] , b[2] );        p2 = Vec3( B[0] , b[1] , b[2] );
			}
			else if( min_edge.n4->depth >= min_edge.n1->depth   &&  min_edge.n4->depth >= min_edge.n2->depth   &&   min_edge.n4->depth >= min_edge.n3->depth ) {
				Vec3 b = min_edge.n4->boundingBox.bmin , B = min_edge.n4->boundingBox.bmax;
				p1 = Vec3( b[0] , b[1] , B[2] );        p2 = Vec3( B[0] , b[1] , B[2] );
			}
		}
		else if( min_edge.axis == 1 ) {
			if( min_edge.n1->depth >= min_edge.n2->depth   &&  min_edge.n1->depth >= min_edge.n3->depth   &&   min_edge.n1->depth >= min_edge.n4->depth ) {
				Vec3 b = min_edge.n1->boundingBox.bmin , B = min_edge.n1->boundingBox.bmax;
				p1 = Vec3( B[0] , b[1] , b[2] );        p2 = Vec3( B[0] , B[1] , b[2] );
			}
			else if( min_edge.n2->depth >= min_edge.n1->depth   &&  min_edge.n2->depth >= min_edge.n3->depth   &&   min_edge.n2->depth >= min_edge.n4->depth ) {
				Vec3 b = min_edge.n2->boundingBox.bmin , B = min_edge.n2->boundingBox.bmax;
				p1 = Vec3( b[0] , b[1] , b[2] );        p2 = Vec3( b[0] , B[1] , b[2] );
			}
			else if( min_edge.n3->depth >= min_edge.n1->depth   &&  min_edge.n3->depth >= min_edge.n2->depth   &&   min_edge.n3->depth >= min_edge.n4->depth ) {
				Vec3 b = min_edge.n3->boundingBox.bmin , B = min_edge.n3->boundingBox.bmax;
				p1 = Vec3( B[0] , b[1] , B[2] );        p2 = Vec3( B[0] , B[1] , B[2] );
			}
			else if( min_edge.n4->depth >= min_edge.n1->depth   &&  min_edge.n4->depth >= min_edge.n2->depth   &&   min_edge.n4->depth >= min_edge.n3->depth ) {
				Vec3 b = min_edge.n4->boundingBox.bmin , B = min_edge.n4->boundingBox.bmax;
				p1 = Vec3( b[0] , b[1] , B[2] );        p2 = Vec3( b[0] , B[1] , B[2] );
			}
		}
		else if( min_edge.axis == 2 ) {
			if( min_edge.n1->depth >= min_edge.n2->depth   &&  min_edge.n1->depth >= min_edge.n3->depth   &&   min_edge.n1->depth >= min_edge.n4->depth ) {
				Vec3 b = min_edge.n1->boundingBox.bmin , B = min_edge.n1->boundingBox.bmax;
				p1 = Vec3( B[0] , B[1] , b[2] );        p2 = Vec3( B[0] , B[1] , B[2] );
			}
			else if( min_edge.n2->depth >= min_edge.n1->depth   &&  min_edge.n2->depth >= min_edge.n3->depth   &&   min_edge.n2->depth >= min_edge.n4->depth ) {
				Vec3 b = min_edge.n2->boundingBox.bmin , B = min_edge.n2->boundingBox.bmax;
				p1 = Vec3( b[0] , B[1] , b[2] );        p2 = Vec3( b[0] , B[1] , B[2] );
			}
			else if( min_edge.n3->depth >= min_edge.n1->depth   &&  min_edge.n3->depth >= min_edge.n2->depth   &&   min_edge.n3->depth >= min_edge.n4->depth ) {
				Vec3 b = min_edge.n3->boundingBox.bmin , B = min_edge.n3->boundingBox.bmax;
				p1 = Vec3( B[0] , b[1] , b[2] );        p2 = Vec3( B[0] , b[1] , B[2] );
			}
			else if( min_edge.n4->depth >= min_edge.n1->depth   &&  min_edge.n4->depth >= min_edge.n2->depth   &&   min_edge.n4->depth >= min_edge.n3->depth ) {
				Vec3 b = min_edge.n4->boundingBox.bmin , B = min_edge.n4->boundingBox.bmax;
				p1 = Vec3( b[0] , b[1] , b[2] );        p2 = Vec3( b[0] , b[1] , B[2] );
			}
		}
	}

	template< class F >
	void execute_on_minimal_edges( F f ) {
		RegularOctreeNode::CellProc_traverse_minimal_edges( &root , f );
	}

	void gather_minimal_edges( std::vector< Edge > & minimal_edges ) {
		struct record_minimal_edges_callback {
			std::vector< Edge > & min_edges;
			void operator()(
						RegularOctreeNode * _n1 , RegularOctreeNode * _n2 ,
						RegularOctreeNode * _n3 , RegularOctreeNode * _n4,
						unsigned char _axis  ) {
				min_edges.push_back( Edge(_n1,_n2,_n3,_n4,_axis) );
			}
			record_minimal_edges_callback( std::vector< Edge > & _min_edges ) : min_edges(_min_edges) {  }
		};
		record_minimal_edges_callback f( minimal_edges );
		execute_on_minimal_edges( f );
	}

	template< class F >
	void dual_contouring( F scalar_field , std::vector< Vec3 > & mesh_vertices , std::vector< unsigned int > & mesh_triangles
						  , std::vector< Vec3 > & dc_minimal_edges // useful only for debug
						  , int normalEstimationMode // 0 : p2-p1   ,  1 : linear interp of normal1 and normal2   ,  2 : normal at pInterp
						  ) {
		std::vector< Edge > minimal_edges;
		gather_minimal_edges( minimal_edges );

		mesh_vertices.clear();
		mesh_triangles.clear();
		std::vector< RegularOctreeNode * > node_corresponding_to_mesh_vertex;

		// for every minimal edge:
		for( Edge min_edge : minimal_edges ) {
			Vec3 p1 , p2;
			locate_minimal_edge( min_edge , p1 , p2 );

			// check if f changes sign
			double f1 = scalar_field(p1);
			double f2 = scalar_field(p2);
			Vec3 normal1 , normal2;
			if( normalEstimationMode == 1 ){
				f1 = scalar_field(p1 , normal1);
				f2 = scalar_field(p2 , normal2);
			}
			else {
				f1 = scalar_field(p1);
				f2 = scalar_field(p2);
			}

			// if f changes sign, create a new point
			if( f1*f2 < 0.0 ){
				// record traversed minimal edge for debug:
				dc_minimal_edges.push_back( p1 ); dc_minimal_edges.push_back( p2 );

				bool flipPatch = f2 < f1;

				DCOctree::RegularOctree::RegularOctreeNode * n1 = min_edge.getN1();
				DCOctree::RegularOctree::RegularOctreeNode * n2 = min_edge.getN2();
				DCOctree::RegularOctree::RegularOctreeNode * n3 = min_edge.getN3();
				DCOctree::RegularOctree::RegularOctreeNode * n4 = min_edge.getN4();
				Vec3 plane_center_to_fit = (f2 * p1 - f1 * p2) / (f2 - f1);
				// Vec3 plane_center_to_fit = (p1 + p2) / (2.0);


				Vec3 normal_to_fit = normalize(p2 - p1);
				if(normalEstimationMode == 1) normal_to_fit = normalize(f2 * normal1 - f1 * normal2);
				else if(normalEstimationMode == 2) {
					scalar_field(plane_center_to_fit , normal_to_fit);
				}

				QEM toAdd = QEM::fromPlan( normal_to_fit * length(p2 - p1), plane_center_to_fit );
				n1->qem += toAdd;
				n2->qem += toAdd;
				n3->qem += toAdd;
				n4->qem += toAdd;

				// add tris
				// 1 2 3    &    2 4 3    (if ! flipPatch)
				unsigned int i1 , i2 , i3 , i4;

				if( n1->idx_point_in_mesh == UINT_MAX ) {
					mesh_vertices.push_back( n1->boundingBox.center() );
					node_corresponding_to_mesh_vertex.push_back( n1 );
					assert( node_corresponding_to_mesh_vertex.size() == mesh_vertices.size() );
					i1 = mesh_vertices.size() - 1;
					n1->idx_point_in_mesh = i1;
				}
				else i1 = n1->idx_point_in_mesh;

				if( n2->idx_point_in_mesh == UINT_MAX ) {
					mesh_vertices.push_back( n2->boundingBox.center() );
					node_corresponding_to_mesh_vertex.push_back( n2 );
					assert( node_corresponding_to_mesh_vertex.size() == mesh_vertices.size() );
					i2 = mesh_vertices.size() - 1;
					n2->idx_point_in_mesh = i2;
				}
				else i2 = n2->idx_point_in_mesh;

				if( n3->idx_point_in_mesh == UINT_MAX ) {
					mesh_vertices.push_back( n3->boundingBox.center() );
					node_corresponding_to_mesh_vertex.push_back( n3 );
					assert( node_corresponding_to_mesh_vertex.size() == mesh_vertices.size() );
					i3 = mesh_vertices.size() - 1;
					n3->idx_point_in_mesh = i3;
				}
				else i3 = n3->idx_point_in_mesh;

				if( n4->idx_point_in_mesh == UINT_MAX ) {
					mesh_vertices.push_back( n4->boundingBox.center() );
					node_corresponding_to_mesh_vertex.push_back( n4 );
					assert( node_corresponding_to_mesh_vertex.size() == mesh_vertices.size() );
					i4 = mesh_vertices.size() - 1;
					n4->idx_point_in_mesh = i4;
				}
				else i4 = n4->idx_point_in_mesh;

				if( i1 != i2 && i1 != i3 && i3 != i2 ){
					if( flipPatch ) {
						mesh_triangles.push_back( i1 );
						mesh_triangles.push_back( i3 );
						mesh_triangles.push_back( i2 );
					}
					else {
						mesh_triangles.push_back( i1 );
						mesh_triangles.push_back( i2 );
						mesh_triangles.push_back( i3 );
					}
				}

				if( i4 != i2 && i2 != i3 && i3 != i4 ){
					if( flipPatch ) {
						mesh_triangles.push_back( i2 );
						mesh_triangles.push_back( i3 );
						mesh_triangles.push_back( i4 );
					}
					else {
						mesh_triangles.push_back( i2 );
						mesh_triangles.push_back( i4 );
						mesh_triangles.push_back( i3 );
					}
				}
			}
		}

		for( unsigned int v = 0 ; v < mesh_vertices.size() ; ++v ) {
			mesh_vertices[v] = node_corresponding_to_mesh_vertex[v]->qem.minimizer( mesh_vertices[v] );
		}
	}
};
}


