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

#include "src/Vec3.h"

#include "../apss/Fastapss.h"

extern std::vector< Vec3 > mesh_vertices;

extern bool drawDCME;
extern std::vector< Vec3 > dc_minimal_edges;


// ------------------------------------------------------------------------------------------------------------
// Grid
// ------------------------------------------------------------------------------------------------------------
struct GridPoint {
	double f;
};
struct GridCell {
	unsigned int idx_in_mesh_vertices;
};


struct Grid {
	// size of grid
	unsigned int nX,nY,nZ;

	// bounding box
	Vec3 bb , BB;

	// points (not centers of cells)
	std::vector< GridPoint > points;
	std::vector< GridCell > cells;


	Grid() {}
	Grid(unsigned int x,unsigned int y,unsigned int z , Vec3 _bb , Vec3 _bbb){
		assert(x > 0 && y > 0 && z > 0);
		nX=x;
		nY=y;
		nZ=z;
		bb=_bb;
		BB=_bbb;
		points.resize(nX*nY*nZ);
		cells.resize((nX-1)*(nY-1)*(nZ-1));
	}

	size_t pointSerializedIndex( unsigned int x , unsigned int y , unsigned int z ) const {
		return x + y*nX + z*nX*nY;
	}
	size_t cellSerializedIndex( unsigned int x , unsigned int y , unsigned int z ) const {
		return x + y*(nX-1) + z*(nX-1)*(nY-1);
	}

	GridPoint & point(unsigned int x,unsigned int y,unsigned int z){
		return points[ pointSerializedIndex(x,y,z) ];
	}
	GridCell & cell(unsigned int x,unsigned int y,unsigned int z){
		return cells[ cellSerializedIndex(x,y,z) ];
	}

	// return position of a point in the grid
	Vec3 position(float x,float y,float z){ // x \in [0, nX] , ...
		return Vec3( bb[0] + (x/nX) * (BB[0]-bb[0]) , bb[1] + (y/nY) * (BB[1]-bb[1]) , bb[2] + (z/nZ) * (BB[2]-bb[2]));
	}
};



// ------------------------------------------------------------------------------------------------------------
// Dual contouring on Grid
// ------------------------------------------------------------------------------------------------------------

// mode = 1 : new point at center of cell
// mode = 2 : new point at center projected to MLS surface
void dualContouring( Grid & grid , std::vector< Vec3 > &  new_positions , std::vector<unsigned int> & new_triangles);

void smoothMesh( std::vector< Vec3 > & iPoints , std::vector< uint32_t > & iTris);

void subdivideMesh( std::vector< Vec3 > & iPoints , std::vector< unsigned int > & iTris);

void dualContour(APSS * apss, Kernel * ker, std::vector<glm::vec3> & points, std::vector<glm::vec3> & vertices, std::vector<glm::vec3> & normals, std::vector<uint32_t> & mesh_triangles);
