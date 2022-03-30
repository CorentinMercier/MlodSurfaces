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

#include "GPUUtils.h"
#include "Kernel.h"

#include <vector>
#include <set>
#include <queue>

#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>

#define MAX_DEPTH 12
#define SHARED_MEM

//#define EXPFUNC

struct octreeNodeGPU2;
struct apssNodeStatsGPUSoA;

class Fastapss: public APSS
{
public:
	Fastapss(APSSOctree *apssOctree);
	Fastapss(APSSOctree *apssOctree, const Kernel &kernel, unsigned int nbOfVectors);
	virtual ~Fastapss(){eraseFromGPU();}

	void project(unsigned int nbOfVectors, glm::vec3 *outputPoints, glm::vec3 *outputNormals, unsigned int n_iterations, const Kernel &kernel, float stepMaxSize = -1.f) const final override;
	void copyToGPU(const Kernel &kernel, unsigned int nbOfVectors) final override;
	void eraseFromGPU() final override;

	void projectCPU(glm::vec3 const & qStart , glm::vec3 & outputPoint , glm::vec3 & outputNormal , unsigned int n_iterations ,
				  Kernel const & kernel = Kernel(), float stepMaxSize = -1.f) const final override;

private:
	octreeNodeGPU2 * m_nodesBis;
	apssNodeStatsGPUSoA * m_apssNodeStatsCudaSoA;
};
