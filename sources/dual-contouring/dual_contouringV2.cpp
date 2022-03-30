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

#include "dual_contouringV2.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <queue>
#include <map>
#include <algorithm>
#include <string>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <float.h>
#include "src/Vec3.h"
#include "kdTreeNN/kdTreeNN.h"

#include "DCOctree.h"
#include "../apss/Fastapss.h"

std::vector< Vec3 > mesh_vertices;


bool drawDCME = true;
std::vector< Vec3 > dc_minimal_edges;


// ----

// ------------------------------------------------------------------------------------------------------------
// Dual contouring on Grid
// ------------------------------------------------------------------------------------------------------------

// mode = 1 : new point at center of cell
// mode = 2 : new point at center projected to MLS surface
void dualContouring( Grid & grid , std::vector< Vec3 > &  new_positions , std::vector<unsigned int> & new_triangles ){

	std::cout<<"Dual contouring start."<<std::endl;
	// triangle mesh
	new_positions.clear();
	new_triangles.clear();

	// for every cell
	for(unsigned int x=0;x < grid.nX -1;x++){
		for(unsigned int y=0;y < grid.nY -1;y++){
			for(unsigned int z=0;z < grid.nZ -1;z++){

				// check if f changes sign
				float f_min=std::min(std::min(std::min(grid.point(x,y,z).f,grid.point(x+1,y,z).f),std::min(grid.point(x,y+1,z).f,grid.point(x,y,z+1).f)) ,
									 std::min(std::min(grid.point(x+1,y+1,z).f,grid.point(x+1,y,z+1).f),std::min(grid.point(x,y+1,z+1).f,grid.point(x+1,y+1,z+1).f)));

				float f_max=std::max(std::max(std::max(grid.point(x,y,z).f,grid.point(x+1,y,z).f),std::max(grid.point(x,y+1,z).f,grid.point(x,y,z+1).f)) ,
									 std::max(std::max(grid.point(x+1,y+1,z).f,grid.point(x+1,y,z+1).f),std::max(grid.point(x,y+1,z+1).f,grid.point(x+1,y+1,z+1).f))); // check THE CORNERS

				// if f changes sign, create a new point
				if(f_min*f_max<0){

					// new point at center of cell
					{
						new_positions.push_back(grid.position(x+0.5,y+0.5,z+0.5));
					}

					// record the created point's index
					grid.cell(x,y,z).idx_in_mesh_vertices = new_positions.size()-1;
				}
			}
		}
	}
	std::cout << "DC has detected " << new_positions.size() << " mesh vertices (among " << grid.nX*grid.nY*grid.nZ << " possible)" << std::endl;

	// for every edge in direction x
	for(unsigned int x=0 ; x < grid.nX -1 ; x++){
		for(unsigned int y=1 ; y < grid.nY -1 ; y++){
			for(unsigned int z=1 ; z < grid.nZ -1 ; z++){

				// if f changes sign, append 2 triangles to the triangle list
				if(grid.point(x,y,z).f * grid.point(x+1,y,z).f < 0.0)
				{
					unsigned int v1,v2,v3,v4;
					v1=grid.cell(x,y-1,z-1).idx_in_mesh_vertices;
					v2=grid.cell(x,y-1,z).idx_in_mesh_vertices;
					v3=grid.cell(x,y,z-1).idx_in_mesh_vertices;
					v4=grid.cell(x,y,z).idx_in_mesh_vertices;

					if( v1 >= new_positions.size() ) { std::cout << v1 << std::endl; assert(0); }
					if( v2 >= new_positions.size() ) { std::cout << v2 << std::endl; assert(0); }
					if( v3 >= new_positions.size() ) { std::cout << v3 << std::endl; assert(0); }
					if( v4 >= new_positions.size() ) { std::cout << v4 << std::endl; assert(0); }

					// decide orientation of triangle
					if(grid.point(x,y,z).f < 0.0)
					{
						new_triangles.push_back(v1);
						new_triangles.push_back(v3);
						new_triangles.push_back(v4);

						new_triangles.push_back(v1);
						new_triangles.push_back(v4);
						new_triangles.push_back(v2);
					}
					else
					{
						new_triangles.push_back(v1);
						new_triangles.push_back(v4);
						new_triangles.push_back(v3);

						new_triangles.push_back(v1);
						new_triangles.push_back(v2);
						new_triangles.push_back(v4);
					}
				}
			}
		}
	}

	// same for edges in direction y
	for(unsigned int x=1;x<grid.nX -1;x++){
		for(unsigned int y=0;y<grid.nY -1;y++){
			for(unsigned int z=1;z<grid.nZ -1;z++){
				if(grid.point(x,y,z).f * grid.point(x,y+1,z).f < 0.0)
				{
					unsigned int v1,v2,v3,v4;
					v1=grid.cell(x-1,y,z-1).idx_in_mesh_vertices;
					v2=grid.cell(x,y,z-1).idx_in_mesh_vertices;
					v3=grid.cell(x-1,y,z).idx_in_mesh_vertices;
					v4=grid.cell(x,y,z).idx_in_mesh_vertices;

					if( v1 >= new_positions.size() ) { std::cout << v1 << std::endl; assert(0); }
					if( v2 >= new_positions.size() ) { std::cout << v2 << std::endl; assert(0); }
					if( v3 >= new_positions.size() ) { std::cout << v3 << std::endl; assert(0); }
					if( v4 >= new_positions.size() ) { std::cout << v4 << std::endl; assert(0); }

					if(grid.point(x,y,z).f<0.0)
					{
						new_triangles.push_back(v1);
						new_triangles.push_back(v3);
						new_triangles.push_back(v4);
						new_triangles.push_back(v1);
						new_triangles.push_back(v4);
						new_triangles.push_back(v2);
					}
					else
					{
						new_triangles.push_back(v1);
						new_triangles.push_back(v4);
						new_triangles.push_back(v3);
						new_triangles.push_back(v1);
						new_triangles.push_back(v2);
						new_triangles.push_back(v4);
					}
				}
			}
		}
	}

	// same for edges in direction z
	for(unsigned int x=1;x<grid.nX -1;x++){
		for(unsigned int y=1;y<grid.nY -1;y++){
			for(unsigned int z=0;z<grid.nZ -1;z++){
				if(grid.point(x,y,z).f * grid.point(x,y,z+1).f < 0.0)
				{
					unsigned int v1,v2,v3,v4;
					v1=grid.cell(x-1,y-1,z).idx_in_mesh_vertices;
					v2=grid.cell(x-1,y,z).idx_in_mesh_vertices;
					v3=grid.cell(x,y-1,z).idx_in_mesh_vertices;
					v4=grid.cell(x,y,z).idx_in_mesh_vertices;

					if( v1 >= new_positions.size() ) { std::cout << v1 << std::endl; assert(0); }
					if( v2 >= new_positions.size() ) { std::cout << v2 << std::endl; assert(0); }
					if( v3 >= new_positions.size() ) { std::cout << v3 << std::endl; assert(0); }
					if( v4 >= new_positions.size() ) { std::cout << v4 << std::endl; assert(0); }

					if(grid.point(x,y,z).f<0.0)
					{
						new_triangles.push_back(v1);
						new_triangles.push_back(v3);
						new_triangles.push_back(v4);
						new_triangles.push_back(v1);
						new_triangles.push_back(v4);
						new_triangles.push_back(v2);
					}
					else
					{
						new_triangles.push_back(v1);
						new_triangles.push_back(v4);
						new_triangles.push_back(v3);
						new_triangles.push_back(v1);
						new_triangles.push_back(v2);
						new_triangles.push_back(v4);
					}
				}
			}
		}
	}

	std::cout<<"Dual contouring finished."<<std::endl;
}

void smoothMesh( std::vector< Vec3 > & iPoints , std::vector< uint32_t > & iTris ) {
	std::vector< Vec3 > newPos( iPoints.size() , Vec3(0,0,0));
	std::vector< double > valence( iPoints.size() , 0 );

	const unsigned int iTrisNb = iTris.size() / 3;
	for( unsigned int t = 0 ; t < iTrisNb ; ++t ) {
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
	for( unsigned int v = 0 ; v < iPoints.size() ; ++v ) {
		iPoints[v] = newPos[v] / valence[v];
	}
}

void subdivideMesh( std::vector< Vec3 > & iPoints , std::vector< unsigned int > & iTris ) {
	struct Edge{
		unsigned int a , b;
		Edge(){}
		Edge(unsigned int c , unsigned int d) : a( std::min(c,d) ) , b( std::max(c,d) ) {}
		bool operator < (Edge const & o) const { return a < o.a || (a == o.a  &&  b < o.b) ; }
	};
	std::map< Edge , unsigned int > newVertOnEdge;
	const unsigned int iTrisNb = iTris.size() / 3;
	for( unsigned int t = 0 ; t < iTrisNb ; ++t ) {
		unsigned int a = iTris[3*t];
		unsigned int b = iTris[3*t + 1];
		unsigned int c = iTris[3*t + 2];
		unsigned int ab , bc , ca;
		std::map< Edge , unsigned int >::iterator ab_vert = newVertOnEdge.find( Edge(a,b) );
		if( ab_vert ==  newVertOnEdge.end() ) {
			iPoints.push_back( (iPoints[a] + iPoints[b])/2.0 );
			ab = iPoints.size() - 1;
			newVertOnEdge[ Edge(a,b) ] = ab;
		}
		else ab = ab_vert->second;
		std::map< Edge , unsigned int >::iterator bc_vert = newVertOnEdge.find( Edge(c,b) );
		if( bc_vert ==  newVertOnEdge.end() ) {
			iPoints.push_back( (iPoints[c] + iPoints[b])/2.0 );
			bc = iPoints.size() - 1;
			newVertOnEdge[ Edge(c,b) ] = bc;
		}
		else bc = bc_vert->second;
		std::map< Edge , unsigned int >::iterator ca_vert = newVertOnEdge.find( Edge(c,a) );
		if( ca_vert ==  newVertOnEdge.end() ) {
			iPoints.push_back( (iPoints[c] + iPoints[a])/2.0 );
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


void dualContour(APSS * apss, Kernel * ker, std::vector<glm::vec3> & points, std::vector<glm::vec3> & vertices, std::vector<glm::vec3> & normals, std::vector<uint32_t> & mesh_triangles) {
	DCOctree::RegularOctree dc_octree;

	unsigned int dc_octree_max_depth_for_leaves = 6;

	std::vector< Vec3 > points_guiding_dual_contouring(points.size());
	for (unsigned i=0; i<points.size(); i++)
		points_guiding_dual_contouring[i] = Vec3(points[i].x, points[i].y, points[i].z);

	Vec3 gridbb = points_guiding_dual_contouring[0];
	Vec3 gridBB = points_guiding_dual_contouring[0];
	for( unsigned int pIt = 0 ; pIt < points_guiding_dual_contouring.size() ; ++pIt ) {
		for( unsigned int c = 0 ; c < 3 ; ++c ) {
			gridbb[c] = std::min(gridbb[c] , points_guiding_dual_contouring[pIt][c]);
			gridBB[c] = std::max(gridBB[c] , points_guiding_dual_contouring[pIt][c]);
		}
	}

	Vec3 gridCenter = (gridbb + gridBB) / 2.0;
	gridbb += 0.1f * (gridbb - gridCenter);
	gridBB += 0.1f * (gridBB - gridCenter);

	dc_octree.setBoundingBox( gridbb , gridBB ); // MANDATORY

	std::cout << "DCOctree::RegularOctree::build" << std::endl;
	dc_octree.build( points_guiding_dual_contouring , false , dc_octree_max_depth_for_leaves ); // MANDATORY
	std::cout << "DCOctree::RegularOctree::enforceMinimalDepth" << std::endl;
	dc_octree.enforceMinimalDepth(4); // OPTIONAL
	std::cout << "DCOctree::RegularOctree::index_nodes_by_breadth_first_search" << std::endl;
	dc_octree.index_nodes_by_breadth_first_search(); // MANDATORY


	struct scalar_field_apss_two_steps {
		// it will need to access all of that:
		std::vector< DCOctree::RegularOctree::Edge > minimal_edges;

		std::vector< Vec3 > minimal_edges_extremities;
		std::vector< Vec3 > minimal_edges_extremities_projection;
		std::vector< Vec3 > minimal_edges_extremities_projection_normal;
	};


	mesh_vertices.clear();
	mesh_triangles.clear();
	std::vector< DCOctree::RegularOctree::RegularOctreeNode * > node_corresponding_to_mesh_vertex;

	scalar_field_apss_two_steps fImpl2steps;
	dc_octree.gather_minimal_edges( fImpl2steps.minimal_edges );
	fImpl2steps.minimal_edges_extremities.resize( 2 * fImpl2steps.minimal_edges.size() );
	for( unsigned int eIt = 0 ; eIt < fImpl2steps.minimal_edges.size() ; ++eIt ) {
		DCOctree::RegularOctree::Edge const & e = fImpl2steps.minimal_edges[eIt];
		dc_octree.locate_minimal_edge( e , fImpl2steps.minimal_edges_extremities[2 * eIt] , fImpl2steps.minimal_edges_extremities[2 * eIt + 1] );
	}

	fImpl2steps.minimal_edges_extremities_projection.resize( 2 * fImpl2steps.minimal_edges.size() );
	fImpl2steps.minimal_edges_extremities_projection_normal.resize( 2 * fImpl2steps.minimal_edges.size() );

	// HERE, PROJECT ALL POINTS minimal_edges_extremities IN PARALLEL:
	std::vector<glm::vec3> outputPoints(fImpl2steps.minimal_edges_extremities.size()), outputNormals(fImpl2steps.minimal_edges_extremities.size());
	std::vector<glm::vec3> inputPoints(fImpl2steps.minimal_edges_extremities.size());
	for (unsigned i=0; i<fImpl2steps.minimal_edges_extremities.size(); i++)
		inputPoints[i] = glm::vec3(fImpl2steps.minimal_edges_extremities[i][0], fImpl2steps.minimal_edges_extremities[i][1], fImpl2steps.minimal_edges_extremities[i][2]);
	apss->erasePointsFromGPU();
	apss->copyPointsToGPU(fImpl2steps.minimal_edges_extremities.size(), inputPoints.data());
	apss->project(fImpl2steps.minimal_edges_extremities.size() , outputPoints.data(), outputNormals.data(), 1, *ker);
	for (unsigned i=0; i<fImpl2steps.minimal_edges_extremities.size(); i++)
	{
		fImpl2steps.minimal_edges_extremities_projection[i] = Vec3(outputPoints[i].x, outputPoints[i].y, outputPoints[i].z);
		fImpl2steps.minimal_edges_extremities_projection_normal[i] = Vec3(outputNormals[i].x, outputNormals[i].y, outputNormals[i].z);
	}

	for( unsigned int eIt = 0 ; eIt < fImpl2steps.minimal_edges.size() ; ++eIt ) {
		DCOctree::RegularOctree::Edge & min_edge = fImpl2steps.minimal_edges[eIt];

		Vec3 p1 = fImpl2steps.minimal_edges_extremities[2 * eIt] ,
				p2 = fImpl2steps.minimal_edges_extremities[2 * eIt + 1];

		Vec3 p1_proj = fImpl2steps.minimal_edges_extremities_projection[2 * eIt] ,
				p2_proj = fImpl2steps.minimal_edges_extremities_projection[2 * eIt + 1];

		Vec3 n1_proj = fImpl2steps.minimal_edges_extremities_projection_normal[2 * eIt] ,
				n2_proj = fImpl2steps.minimal_edges_extremities_projection_normal[2 * eIt + 1];

		// check if f changes sign
		double f1 = dot( p1 - p1_proj , n1_proj );
		double f2 = dot( p2 - p2_proj , n2_proj );

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


			Vec3 normal_to_fit = glm::normalize(f2 * n1_proj - f1 * n2_proj);

			DCOctree::QEM toAdd = DCOctree::QEM::fromPlan( normal_to_fit * length(p2 - p1) , plane_center_to_fit );
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




	for( int k = 0 ; k < 2 ; ++k ) {
		smoothMesh(mesh_vertices , mesh_triangles);
	}

	vertices.resize(mesh_vertices.size());
	normals.resize(mesh_vertices.size());
	for (unsigned i=0; i<mesh_vertices.size(); i++)
		vertices[i] = glm::vec3(mesh_vertices[i][0], mesh_vertices[i][1], mesh_vertices[i][2]);
	apss->erasePointsFromGPU();
	apss->copyPointsToGPU(mesh_vertices.size(), vertices.data());
	apss->project(mesh_vertices.size(), vertices.data(), normals.data(), 1, *ker);
	for (unsigned i=0; i<mesh_vertices.size(); i++)
		mesh_vertices[i] = Vec3(vertices[i].x, vertices[i].y, vertices[i].z);

	for( int k = 0 ; k < 2 ; ++k )
	{
		smoothMesh(mesh_vertices , mesh_triangles);
		subdivideMesh(mesh_vertices,mesh_triangles);
		vertices.resize(mesh_vertices.size());
		normals.resize(mesh_vertices.size());
		for (unsigned i=0; i<mesh_vertices.size(); i++)
			vertices[i] = glm::vec3(mesh_vertices[i][0], mesh_vertices[i][1], mesh_vertices[i][2]);
		apss->erasePointsFromGPU();
		apss->copyPointsToGPU(mesh_vertices.size(), vertices.data());
		apss->project(mesh_vertices.size(), vertices.data(), normals.data(), 1, *ker);
		for (unsigned i=0; i<mesh_vertices.size(); i++)
			mesh_vertices[i] = Vec3(vertices[i].x, vertices[i].y, vertices[i].z);
	}

}
