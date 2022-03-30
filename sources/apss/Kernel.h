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
