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

#ifndef KERNEL_H
#define KERNEL_H

#include "../opengl/openglutils.h"

#include <iostream>
#include <vector>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>

class Kernel
{
public:
	Kernel();
	Kernel(Kernel  const &newKernel);
	void operator =(Kernel const &k);
	float w( glm::vec3 const & p , glm::vec3 const & eta ) const;
	glm::vec4 gradW( glm::vec3 const & p , glm::vec3 const & eta ) const;

	void draw(GLuint program);

	void updateSigmas(bool computedSigmas = true, float radius = 1.f);

	enum KernelType { GAUSSIAN , SINGULAR , GAUSSIAN_MULTIPLE };
	float sigma = 0.0f;
	float exponent = 0.0f;
	float constant_shift = 0.0f;
	std::vector< float > sigmas;

	std::vector<float> fixedSigmas;//Vector to record fixed values
	float beginningSigma = 0.0f;
	float factor = 0.0f;
	unsigned nbOfComputedSigmas =1;

	KernelType type;

private:
	void init();

	glm::vec3 m_center = glm::vec3(0.f, 0.f, 0.f);
	float m_radius = 1.f;

	unsigned m_nbOfPtsForDrawing = 902;

	GLuint m_vao;
	GLuint m_vbo;
};

#endif // KERNEL_H
