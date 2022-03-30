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

#include "CPUUtils.h"
#include "Kernel.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <cudaGL.h>

struct KernelGPU;

class APSS {
public:
	virtual ~APSS(){erasePointsFromGPU();}

	virtual void project(unsigned int nbOfVectors, glm::vec3 * outputPoints , glm::vec3 * outputNormals ,
				 unsigned int n_iterations, Kernel const & kernel, float stepMaxSize = -1.f) const = 0;

	virtual void copyToGPU(const Kernel &kernel, unsigned int nbOfVectors) = 0;

	virtual void eraseFromGPU() = 0;

	virtual void projectCPU(glm::vec3 const & qStart , glm::vec3 & outputPoint , glm::vec3 & outputNormal , unsigned int n_iterations ,
				  Kernel const & kernel = Kernel(), float stepMaxSize = -1.f) const = 0;

	void printApss() const;

	void reorganizePoints(const std::vector<unsigned> *invalidPts, unsigned endValid);

	void setScalingProtectionSphere(float val){m_scalingProtectionSphere = val; m_apssOctree->setScalingProtectionSphere(val);}
	float getScalingProtectionSphere(){return m_scalingProtectionSphere;}
	void setMinDepth(unsigned val){m_minDepth = val; m_apssOctree->setMinDepth(val);}
	unsigned getMinDepth(){return m_minDepth;}

	void stop();
	void eraseAndRestart();

	virtual void copyPointsToGPU(unsigned int nbOfVectors, glm::vec3 * outputPoints);
	virtual void erasePointsFromGPU();

	void updateKernel(Kernel const & kernel);

	const APSSOctree* CPU(){return m_apssOctree;}

protected:
	float3 *m_outPts = NULL;
	float3 *m_outNmls = NULL;

	KernelGPU *m_kernelCuda = NULL;

	unsigned int m_numberOfNodes = 0;

	float m_scalingProtectionSphere = 2.f;
	unsigned m_minDepth = 0;

	APSSOctree* m_apssOctree;
};
