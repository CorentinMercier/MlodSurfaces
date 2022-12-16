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

#include "camera.h"
#include <algorithm>

camera::camera()
{
	m_quat = glm::identity<glm::quat>();
}

void camera::setupCamera()
{
	computeProjection();
	computeModelview();
}

void camera::computeProjection()
{
	float ratio=static_cast<float>((float)m_screenSizeX/m_screenSizeY);
	const float f=1.0f/tan(m_fov/2.0f);
	const float fp=f/ratio;

	const float L=m_znear-m_zfar;
	//ASSERT(abs(L)>epsilon,"z_far-z_near too small");

	const float C=(m_zfar+m_znear)/L;
	const float D=(2.0f*m_zfar*m_znear)/L;

	m_projection = glm::transpose(glm::mat4 (fp, 0,  0, 0,
					 0, f,  0, 0,
					 0, 0,  C, D,
					 0, 0, -1, 0));
}

void camera::computeModelview()
{
	glm::mat4 worldMatrixZoom;
	worldMatrixZoom = glm::transpose(glm::mat4(1, 0, 0, 0,
					   0, 1, 0, 0,
					   0, 0, 1, m_dist,
					   0, 0, 0, 1));
	glm::mat4 worldMatrixTranslation;
	worldMatrixTranslation = glm::transpose(glm::mat4(1, 0, 0, m_translation[0],
							  0, 1, 0, m_translation[1],
							  0, 0, 1, m_translation[2],
							  0, 0, 0, 1));
	glm::mat4 worldMatrixRotation=glm::identity<glm::mat4>();
	worldMatrixRotation = glm::toMat4(m_quat);
	m_modelview=worldMatrixZoom*worldMatrixRotation*worldMatrixTranslation;
	m_normal = glm::transpose(glm::inverse(worldMatrixRotation));
}


void camera::goUp(float dL)
{
	glm::vec3 const y(0,-1,0);
	m_translation += dL*(glm::conjugate(m_quat)*y);
}

void camera::goRight(float dL)
{
	glm::vec3 const x(-1,0,0);
	m_translation += dL*(glm::conjugate(m_quat)*x);
}

void camera::goForward(float dL)
{
	glm::vec3 const z(0,0,1);
	m_translation += dL*(glm::conjugate(m_quat)*z);
}

float camera::projectToDisc(float const x, float const y)
{
	float const n=sqrt(x*x+y*y);
	if(n<m_discRadius*0.707107f)
		return sqrt(pow(m_discRadius,2)-n*n);
	else
		return pow(m_discRadius,2)/(2.0f*n);
}

void camera::rotation(int const x, int const y)
{
	float const xOld=m_xPrevious;
	float const yOld=m_yPrevious;

	float const x0=(2.0f*xOld-m_screenSizeX)/m_screenSizeX;
	float const y0=(m_screenSizeY-2.0f*yOld)/m_screenSizeY;
	float const x1=(2.0f*x-m_screenSizeX)/m_screenSizeX;
	float const y1=(m_screenSizeY-2.0f*y)/m_screenSizeY;
	if(sqrt(pow(x0-x1,2)+pow(y0-y1,2)>1e-6))
	{
		glm::vec3 const p1=glm::vec3(x0, y0, projectToDisc(x0, y0));
		glm::vec3 const p2=glm::vec3(x1, y1, projectToDisc(x1, y1));
		glm::vec3 const axis=glm::normalize(glm::cross(p1,p2));
		glm::vec3 const u=p1-p2;
		float t=glm::length(u)/(2.0f*m_discRadius);
		t=min(max(t,-1.0f),1.0f); //clamp
		float const phi = 2.0f*asin(t);
		//compute quaternion to apply
		m_quat = glm::angleAxis(phi, axis) * m_quat;
	}
	m_xPrevious=x;
	m_yPrevious=y;
}

void camera::alignX()
{
	m_quat = glm::angleAxis(float(M_PI_2), glm::vec3(0.f, -1.f, 0.f)) * glm::identity<glm::quat>();
}

void camera::alignY()
{
	m_quat = glm::angleAxis(float(M_PI_2), glm::vec3(1.f, 0.f, 0.f)) * glm::identity<glm::quat>();
}

void camera::alignZ()
{
	m_quat = glm::identity<glm::quat>();
}

void camera::zoom(int const y)
{
	m_dist += (fabs(m_dist)+1.0f)*(y-m_yPrevious)/500.0f;
	m_dist = min(m_dist,0.0f);
	m_yPrevious=y;
	m_zoomPrevious = y;
}

void camera::zoomWheel(int const y)
{
	m_dist += (fabs(m_dist)+1.0f)*(y-m_zoomPrevious)/500.0f;
	m_dist = min(m_dist,0.0f);
	m_zoomPrevious=y;
}

glm::vec3 camera::getUpVector()
{
	return glm::rotate(glm::inverse(m_quat), glm::vec3(0, 1, 0));
}

glm::vec3 camera::getLookAtVector()
{
	return glm::rotate(glm::inverse(m_quat), glm::vec3(0, 0, -1));
}

glm::vec3 camera::getRightVector()
{
	return glm::rotate(glm::inverse(m_quat), glm::vec3(1, 0, 0));
}

glm::vec3 camera::get3DPointFromClick(int x, int y, float z)
{
	glm::mat4 inverseProjection = glm::inverse(m_projection * m_modelview);
	glm::vec4 pt = glm::vec4(((2.f * x)/float(m_screenSizeX) - 1.f), -1.f * ((2.f * y)/float(m_screenSizeY) - 1.f), z, 1);
	pt = inverseProjection * pt;
	pt /= pt.w;
	return glm::vec3(pt);
}
