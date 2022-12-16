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

struct octreeNodeGPU;
struct apssNodeStatsGPU;
struct apssStatsGPU;

class GlobalAPSS final: public APSS {
public:
	GlobalAPSS(APSSOctree *apssOctree);
	GlobalAPSS(APSSOctree *apssOctree, const Kernel &kernel, unsigned int nbOfVectors);
	virtual ~GlobalAPSS(){eraseFromGPU();}

	void project(unsigned int nbOfVectors, glm::vec3 * outputPoints , glm::vec3 * outputNormals ,
							  unsigned int n_iterations, Kernel const & kernel, float stepMaxSize = -1.f) const final override;
	void copyToGPU(const Kernel &kernel, unsigned int nbOfVectors) final override;
	void eraseFromGPU() final override;

	void copyPointsToGPU(unsigned int nbOfVectors, glm::vec3 * outputPoints) final override;
	void erasePointsFromGPU() final override;

	void projectCPU(glm::vec3 const & qStart , glm::vec3 & outputPoint , glm::vec3 & outputNormal , unsigned int n_iterations ,
				  Kernel const & kernel = Kernel(), float stepMaxSize = -1.f) const final override;

	void cpuSphere(glm::vec3 const & qStart, float &u0, glm::vec3 &u123, float &u4, Kernel const & kernel = Kernel()) const;

private:

	octreeNodeGPU * m_octreeNodeCuda = NULL;
	apssNodeStatsGPU *m_apssNodeStatsCuda = NULL;
	apssStatsGPU *m_treeStats = NULL;
};
