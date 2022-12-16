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

#ifndef SCENE_H
#define SCENE_H

#include <string>
#include <iostream>
#include <omp.h>
#include <cmath>
#include <chrono>
#include "../opengl/openglutils.h"
#include "../apss/Pn.h"
#include "../apss/Fastapss.h"
#include "../apss/GlobalAPSS.h"
#include "../apss/timer.h"

#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/component_wise.hpp>

#define COUT_DETAILED_STATISTICS

#ifdef COUT_DETAILED_STATISTICS
#endif

using namespace std;

class WidgetOpenGL;

class scene
{
public:
	scene();
	~scene();

	void reset(); //Free all memory for loading another model

	/**  Method called only once at the beginning (load off files ...) */
	bool load_scene(string filename);

	void compileShaders();
	void initBuffers();

	/**  Method called at every frame */
	void draw_scene();

	void changeMethod(bool fast);

	/** Set the pointer to the parent Widget */
	void set_widget(WidgetOpenGL* widget_param);

	void save_octree(const std::string& filename, unsigned int min_depth = 0, unsigned int max_depth = unsigned(-1), bool use_empty_triangles = true) const;

	void increasePtSize(){m_ptSize++;}
	void decreasePtSize(){if (m_ptSize > 1) m_ptSize--;}

	double projectPoints(bool updateBuffer = true);
	double projectPointsOnCPU(bool updateBuffer = true);

	void eraseFromGPU();

	bool renderPointSet = true;
	bool renderPointsToProject = false;
	bool renderProjectedPoints = true;
	bool renderInvalidPoints = false;
	bool renderSelection = false;
	bool renderKernel = false;
	bool normalMode = false;
	bool renderDualContour = false;
	bool renderComparison = false;
	bool wireframe = false;

	double drawTime(){return m_time.getTime();}
	void setNewKernel(const Kernel &newKernel);

	double resetPoints(){return resetPoints(m_useAdaptiveOctree);}
	double resetPoints(bool useAdaptiveOctree);
	double resetOctree(){return resetOctree(m_useAdaptiveOctree);}
	double resetOctree(bool useAdaptiveOctree);

	void setCutOff(float val){m_apss->setScalingProtectionSphere(val);}
	float getCutOff(){return m_apss->getScalingProtectionSphere();}

	void stop(){m_apss->stop();}

	void resamplePointSet(unsigned int nbOfPoints);
	unsigned int sizeOfPointSet(){return m_originalPn.size();}

	void clusterSkeleton();

	bool upsample();
	bool upsample(double &mortonCodeTime);
	void initUpsampling();

	unsigned getNbOfPtsToProject(){return m_nbOfPoints;}
	unsigned getmaxUpsampling(){return m_maxUpsamples;}
	unsigned getSizeOfGrid(){return m_sizeOfGrid;}
	unsigned getMaxDepth(){return m_maxDepth;}
	bool getCut(){return m_cut;}
	void setNbOfPtsToProject(unsigned nbOfPts){m_nbOfPoints = nbOfPts;}
	void setMaxUpsampling(unsigned maxUpsampling){m_maxUpsamples = maxUpsampling;}
	void setSizeOfGrid(unsigned sizeOfGrid){m_sizeOfGrid = sizeOfGrid;}
	void setMaxDepth(unsigned maxDepth){m_maxDepth = maxDepth;}
	void setCut(bool cut){m_cut = cut;}
	void setOnUpsample(bool upsampleOctree){m_onUpsampleOctree = upsampleOctree;}
	void updateMinDepthValue(unsigned newVal){m_minDepth = newVal; m_apss->setMinDepth(newVal);}
	void setOffset(std::vector<float> offset){for (unsigned i=0; i<6; i++)m_offset[i] = offset[i];}

	void selectPtsCloseToPlane(glm::vec3 normal, glm::vec3 pt, float threshold = 0.02f);
	void selectClosestPoint(glm::vec3 posCam, glm::vec3 farPosCam, float threshold = 0.1f);
	void saveSelection(string filename){m_selection.save(filename);}
	void saveInputFile(string filename){m_originalPn.save(filename);}
	void savePointsToProject(){m_pnProjectedPoints[m_currentProjection]->save("PointsToProject.ply");}
	void loadPointsToProject();

	void computeNoise();

	Kernel const getKernel() const{return *m_weightKernel;}

	void output_traversed_nodes_statistics();

	double computeDualContouring(string filename = "", bool recordOctree = false);
	void createTrianglesForMPU();

	void saveProjectedPoints(string filename);
	void comparePointSet(string filename);

	void checkForInvalidPoints(unsigned & nbOfInvalidPts);

	unsigned getMaxNbIterations()const{return maxNbIterations;}
	unsigned getFramesForSampling()const{return framesForSampling;}
	unsigned getMinDepth()const{return m_minDepth;}
	bool getComputedSigma()const{return computedSigma;}
	bool getResample()const{return resample;}
	const vector<float>& getOffset()const{return m_offset;}

	void mortonTest();

	void knnMode(bool setMode, unsigned knn = 20){m_apssOctree->knnMode(setMode, knn);}
	void ballMode(bool setMode, double ballRadius){m_apssOctree->ballMode(setMode, ballRadius);}

	bool getMode()const{return m_method;}
	size_t count_nodes_on_projected_points();

	double getAdaptiveOctreeTime() const {return m_adaptiveOctreeTime;}
	double getUsualOctreeTime() const {return m_usualOctreeTime;}
	double getUniformizeTime() const {return m_uniformizeTime;}
	unsigned getNbOfCrossedNodes() const {return m_nbOfCrossedNodes;}

private:
	void changeExtension(string & path, const string & newExtension);
	bool loadParameters(string filename);

	unsigned int m_nbOfPoints = 10000;

	bool renderPerPointColor = false;
	bool m_computeSkeleton = false;
	bool m_renderSkeleton = false;

	unsigned m_minDepth = 0;

	void createBuffers();
	void eraseBuffers();
	double updateBuffers(bool orderWithMorton = true);

	//Access to the parent object
	WidgetOpenGL* m_pwidget;
	//Values for UI
	unsigned maxNbIterations = 10;
	unsigned framesForSampling = 4;
	bool computedSigma = true;
	bool resample = false;

	//GL programs
	GLuint m_basicProgram;
	GLuint m_singleColorProgram;
	GLuint m_perPointColorProgram;
	GLuint m_triangleProgram;

	//PointSets
	Pn m_pn;
	Pn m_selection;
	Pn m_originalPn;
	Pn m_pnToProject;//Random points
	Pn m_pnToProject2;//Fitted to the geometry
	std::vector<Pn*> m_pnProjectedPoints;//  = vector<Pn*>(1);
	Pn m_pnForComparison;
	unsigned m_maxUpsamples = 1000;
	unsigned m_currentProjection = 0;
	unsigned m_sizeOfGrid = 64;
	vector<glm::vec3> m_perPointColor;//Color per point

	//Dual contouring
	unsigned m_nbOfIndices = 0;

	//Octree
	APSSOctree *m_apssOctree;
	APSS *m_apss;
	Kernel * m_weightKernel;

	//Octree for upsampling
	bool m_useAdaptiveOctree = true;
	unsigned m_maxDepth = 6;
	bool m_cut = false;
	bool m_onUpsampleOctree = false;
	std::vector<octree_id> m_adaptiveOctree;
	std::vector<octree_id> m_totalUpsamplingNodes;
	std::vector<octree_id> m_usualOctree;
	Pn generateFromUsualOctree();
	Pn generateFromAdaptiveOctree();
	Pn generateFromBoundingBox();
	Pn generateFromOctree(const std::vector<octree_id>& octree);
	void generate_adaptive_octree();
	void generate_usual_octree();

	double m_adaptiveOctreeTime = 0;
	double m_usualOctreeTime = 0;
	double m_uniformizeTime = 0;
	unsigned m_nbOfCrossedNodes = 0;

	//Storage for the VBOs
	GLuint m_vao;//Original pointSet loaded from file
	GLuint m_vao2;//Original pointSet that will be projected (either random or fitted to the geometry)
	GLuint m_vao3;//PointSet projected
	GLuint m_vao5;//Points with indi
	GLuint m_vaoSelection;
	GLuint m_vaoDualContouring;//Contains the triangles, result from the dual contouring
	GLuint m_vaoComparison;
	GLuint m_vboPointSet;
	GLuint m_vboPointSetToProject;
	GLuint m_vboPointSetProjected;
	GLuint m_vboWithColors;
	GLuint m_vboSelection;
	GLuint m_vboDualContouring;
	GLuint m_iboDualContouring;
	GLuint m_vboComparison;

	unsigned int m_ptSize = 1;

	unsigned m_nbOfInvalidPts = 0;

	Timer m_time;

	bool m_random = false;

	std::vector<float> m_offset = std::vector<float>(6, 0.f);

	bool m_method = true; //True = FastAPSS

#ifdef COUT_DETAILED_STATISTICS
	kdTreeNN::kdTree< glm::vec3 > debug_kdtree;
#endif
};

#endif // SCENE_H
