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
	octreeNodeGPU2* m_nodesBisGPU = nullptr;
	octreeNodeGPU2* m_nodesBisCPU = nullptr;
	apssNodeStatsGPUSoA* m_apssNodeStatsGPU = nullptr;
	apssNodeStatsGPUSoA* m_apssNodeStatsCPU = nullptr;
};
