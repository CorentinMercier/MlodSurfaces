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

#ifndef CAMERA_H
#define CAMERA_H

#include <iostream>

#define _USE_MATH_DEFINES
#include <math.h>

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>

using namespace std;

class camera
{
public:
	camera();

	const glm::mat4 getModelview() {return m_modelview;}
	const glm::mat4 getProjection() {return m_projection;}
	const glm::mat4 getNormal() {return m_normal;}

	void setModelview(glm::mat4 modelview) {m_modelview=modelview;}
	void setProjection(glm::mat4 projection) {m_projection=projection;}
	void getNormal(glm::mat4 normal) {m_normal=normal;}

	void setupCamera();
	void computeProjection();
	void computeModelview();

	void setScreenSize(unsigned int width, unsigned int height){m_screenSizeX=width; m_screenSizeY=height;}
	void setDist(float dist){m_dist=dist;}
	void setQuaternion(glm::fquat q){m_quat=q;}
	void setTranslation(glm::vec3 t){m_translation=t;}

	unsigned int getScreenSizeX(){return m_screenSizeX;}
	unsigned int getScreenSizeY(){return m_screenSizeY;}
	float getDist(){return m_dist;}
	glm::fquat getQuat(){return m_quat;}
	glm::vec3 getTranslation(){return m_translation;}

	int &xPrevious(){return m_xPrevious;}
	int &yPrevious(){return m_yPrevious;}

	//Movements of camera
	void goUp(float dL);
	void goRight(float dL);
	void goForward(float dL);

	void rotation(const int x, const int y);
	void zoom(int const y);
	void zoomWheel(int const y);

	void alignX();
	void alignY();
	void alignZ();

	glm::vec3 getUpVector();
	glm::vec3 getLookAtVector();
	glm::vec3 getRightVector();

	glm::vec3 get3DPointFromClick(int x, int y, float z = 1.f);

private:

	float projectToDisc(float const x, float const y);

	//Matrices for 3D transformations
	glm::mat4 m_modelview;
	glm::mat4 m_projection;
	glm::mat4 m_normal;

	//Camera's parameters
	float m_fov=55.0f*M_PI/180.0f;
	float m_znear=1e-1f;
	float m_zfar=500.0f;
	float m_dist=-3.0f;
	unsigned int m_screenSizeX=800;
	unsigned int m_screenSizeY=800;
	glm::vec3 m_translation=glm::vec3(0.0f, 0.0f, 0.0f);

	//Rotations
	glm::fquat m_quat;
	float m_discRadius=0.8;

	//Mouse positions
	int m_xPrevious;
	int m_yPrevious;

	int m_zoomPrevious = -1325;
};

#endif // CAMERA_H
