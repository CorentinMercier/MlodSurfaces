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

#include "GlobalAPSS.h"
#include "GPUUtils.cu"

/// Structures
///
///

struct apssNodeStatsGPU{
	float3 s_ai_pi , s_ai_ni;
	float s_ai_pi_ni , s_ai_pi_pi;
	float s_ai;

	void init()
	{
	s_ai_pi = make_float3(0.0, 0.0, 0.0);
	s_ai_ni = make_float3(0.0, 0.0, 0.0);
	s_ai_pi_ni = 0;
	s_ai_pi_pi = 0;
	s_ai = 0;
	}

	__host__ void copy(apssNodeStats const & a)
	{
	s_ai = a.s_ai;
	s_ai_ni = make_float3(a.s_ai_ni.x, a.s_ai_ni.y, a.s_ai_ni.z);
	s_ai_pi = make_float3(a.s_ai_pi.x, a.s_ai_pi.y, a.s_ai_pi.z);
	s_ai_pi_ni = a.s_ai_pi_ni;
	s_ai_pi_pi = a.s_ai_pi_pi;
	}
};

struct octreeNodeGPU {
	octreeNodeGPU * children[8];
	unsigned int numberOfChildren;
	unsigned int *indices;
	unsigned int indicesSize;
	unsigned int depth;
	BBOXGPU boundingBox;
	void setBoundingBox( const BBOXGPU & bbox ) { boundingBox = bbox; }
	BBOXGPU const & getBoundingBox() const { return boundingBox; }

	apssNodeStatsGPU nodeapssNodeStats;

	__host__ void erase()
	{
		if (indices != 0)
			gpuErrchk( cudaFree(indices));
		for (unsigned int i=0; i<8; i++)
		{
			if (children[i] != NULL)
			cudaFree(children[i]);
		}
		if(children) cudaFree(children);
	}
	__host__ void copy(octreeNode const * a)
	{
		numberOfChildren = a->numberOfChildren;
		depth = a->depth;
		boundingBox.copy(a->boundingBox);
		indicesSize = a->indices.size();
		if (indicesSize != 0)
			gpuErrchk( cudaMallocManaged(&indices, a->indices.size() * sizeof (unsigned int)));
		for (unsigned int i=0; i<a->indices.size(); i++)
			indices[i] = a->indices[i];
		nodeapssNodeStats.copy(a->nodeapssNodeStats);
		gpuErrchk( cudaMallocManaged((void**)(&children), 8 * sizeof (octreeNodeGPU*)));
		for (unsigned int i=0; i<8; i++)
		{
			if (a->children[i] != NULL)
			{
				gpuErrchk( cudaMallocManaged(&children[i], sizeof (octreeNodeGPU)));
				children[i]->copy(a->children[i]);
			}
			else
			{
				children[i] = NULL;
			}
		}
	}
};

GlobalAPSS::GlobalAPSS(APSSOctree *apssOctree)
{
	m_apssOctree = apssOctree;
	std::cout << "GlobalAPSS method used" << std::endl;
}

GlobalAPSS::GlobalAPSS(APSSOctree *apssOctree, const Kernel &kernel, unsigned int nbOfVectors): GlobalAPSS(apssOctree)
{
	copyToGPU(kernel, nbOfVectors);
}


//////////////////////////////////////////////////////////////////////////
///
/// GPU version 0 -> No octree, running on the GPU
///
//////////////////////////////////////////////////////////////////////////

__device__
void APSS( unsigned int nbOfPoints, float3 const & q , KernelGPU const * kernel,
		   apssNodeStatsGPU const * leavesapssNodeStats, apssStatsGPU & treeStats)
{
	for (unsigned int i=0; i<nbOfPoints; i++)
	{
		float w = kernel->w(leavesapssNodeStats[i].s_ai_pi / leavesapssNodeStats[i].s_ai, q);
		treeStats.s_ai_wi = __fmaf_rn(w, leavesapssNodeStats[i].s_ai, treeStats.s_ai_wi);
		treeStats.s_ai_wi_ni = treeStats.s_ai_wi_ni + w * leavesapssNodeStats[i].s_ai_ni;
		treeStats.s_ai_wi_pi = treeStats.s_ai_wi_pi + w * leavesapssNodeStats[i].s_ai_pi;
		treeStats.s_ai_wi_pi_ni = __fmaf_rn(w, leavesapssNodeStats[i].s_ai_pi_ni, treeStats.s_ai_wi_pi_ni);
		treeStats.s_ai_wi_pi_pi = __fmaf_rn(w, leavesapssNodeStats[i].s_ai_pi_pi, treeStats.s_ai_wi_pi_pi);
	}
}

__device__
void computeSphere(float3 const & q, float3 & outputPoint, float3 & outputNormal, apssStatsGPU const & treeStats, float stepMaxSize)
{
	// fit a sphere to the tree:

	// algebraic sphere: u4.||X||^2 + u123.X + u0 = 0
	// geometric sphere: ||X-C||^2 - r^2 = 0
	// geometric plane:  (X-C).n = 0
	float u4 = 0.5*( treeStats.s_ai_wi_pi_ni/*/treeStats.s_ai_wi*/ - dot((treeStats.s_ai_wi_pi/*/treeStats.s_ai_wi*/),(treeStats.s_ai_wi_ni/treeStats.s_ai_wi)) ) /
			( treeStats.s_ai_wi_pi_pi/*/treeStats.s_ai_wi*/ - dot((treeStats.s_ai_wi_pi/*/treeStats.s_ai_wi*/),(treeStats.s_ai_wi_pi/treeStats.s_ai_wi)) );
	float3 u123= (treeStats.s_ai_wi_ni - 2*u4*treeStats.s_ai_wi_pi)/treeStats.s_ai_wi;
	float u0= -(dot(treeStats.s_ai_wi_pi,u123) + u4*treeStats.s_ai_wi_pi_pi)/treeStats.s_ai_wi;

	outputPoint = q;

	if( fabs(u4) < 0.000000000001 ) {
	// then project on a plane (it's a degenerate sphere)
	float3 n = -u123;
	float lambda = ( u0 - dot(outputPoint,n) ) / squaredNorm(n);
	outputPoint = outputPoint + lambda * n;
	outputNormal = treeStats.s_ai_wi_ni;
	outputNormal = normalize(outputNormal);
	}
	else {
	float3 sphere_center = u123/(-2*u4);
	float val = squaredNorm(sphere_center) - u0/u4 > 0.0 ? squaredNorm(sphere_center) - u0/u4 : 0.0;
	float sphere_radius = sqrt( val );

	// projection of the inputpoint onto the sphere:
	float3 pc= outputPoint-sphere_center;
	pc = normalize(pc);
	outputPoint = sphere_center + sphere_radius*pc;

	// compute normal by looking at the gradient there:
	outputNormal = u123 + 2*u4*outputPoint;
	outputNormal = normalize(outputNormal);
	}
	float stepSize = norme(outputPoint - q);
	if (stepSize > stepMaxSize)
	outputPoint = q + (stepMaxSize / stepSize) * (outputPoint - q);
}

__global__
void launchComputeSphere(unsigned int n, float3 const * q, float3 * outputPoint, float3 * outputNormal, apssStatsGPU const * treeStats, float stepMaxSize)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int stride = blockDim.x * gridDim.x;
	for (int i = index; i < n; i+=stride)
	{
	computeSphere(outputPoint[i], outputPoint[i], outputNormal[i], treeStats[i], stepMaxSize);
	}
}

__global__
//__launch_bounds__(MAX_THREADS_PER_BLOCK, MIN_BLOCKS_PER_MP)
void projectCudaWithoutOctree(unsigned int n, unsigned int nbOfNodes, float3 const * outputPoint,
							  KernelGPU const * kernel, apssNodeStatsGPU const * leavesapssNodeStats, apssStatsGPU * treeStats)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int stride = blockDim.x * gridDim.x;
	for (int i = index; i < n; i+=stride)
	{
	APSS(nbOfNodes, outputPoint[i], kernel, leavesapssNodeStats, treeStats[i]);
	}
}

__global__
void initTreeStats(unsigned int n, apssStatsGPU * treeStats)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int stride = blockDim.x * gridDim.x;
	for (int i = index; i < n; i+=stride)
	treeStats[i].init();
}

void GlobalAPSS::project(unsigned int nbOfVectors, glm::vec3 * outputPoints , glm::vec3 * outputNormals , unsigned int n_iterations ,
						 Kernel const & kernel, float stepMaxSize) const
{
	if (nbOfVectors == 0)
	{
		std::cerr << "No points to project" << std::endl;
		return;
	}
	gpuErrchk( cudaProfilerStart() );
	int blockSize = 128;
	dim3 numBlocks = computeNbBlocks(nbOfVectors, blockSize);

	gpuErrchk( cudaPeekAtLastError() );

	//Init treeStats
	initTreeStats<<<numBlocks, blockSize>>>(nbOfVectors, m_treeStats);
	gpuErrchk( cudaDeviceSynchronize());
	//Main loop, treeStats filling
	gpuErrchk( cudaPeekAtLastError() );
	projectCudaWithoutOctree<<<numBlocks, blockSize>>>(nbOfVectors, m_apssOctree->getNbOfLeaves(), m_outPts, m_kernelCuda, m_apssNodeStatsCuda, m_treeStats);
	gpuErrchk( cudaDeviceSynchronize());
	gpuErrchk( cudaPeekAtLastError() );
	//Sphere computation
	if (stepMaxSize < 0)
		stepMaxSize =  m_apssOctree->getBoundingBox().radius();
	launchComputeSphere<<<numBlocks, blockSize>>>(nbOfVectors, m_outPts, m_outPts, m_outNmls, m_treeStats, stepMaxSize);
	gpuErrchk( cudaDeviceSynchronize());
	gpuErrchk( cudaPeekAtLastError() );

	gpuErrchk( cudaMemcpy(outputPoints, m_outPts, nbOfVectors * sizeof (float3), cudaMemcpyDeviceToHost));
	gpuErrchk( cudaMemcpy(outputNormals, m_outNmls, nbOfVectors * sizeof (float3), cudaMemcpyDeviceToHost));
	cudaProfilerStop();
}

void GlobalAPSS::copyToGPU(const Kernel &kernel, unsigned int nbOfVectors)
{
	//Allocate memory on GPU
	gpuErrchk( cudaPeekAtLastError() );
	std::cout << "Memory is being copied to GPU memory... " << std::flush;
	m_numberOfNodes = m_apssOctree->get_root()->nbOfNodes();
	std::cout << "Number of nodes : " << m_numberOfNodes << std::endl;
	unsigned pointsapssNodeStatsSize = m_apssOctree->getPointsApssNodeStats().size();
	gpuErrchk( cudaMallocManaged(&m_apssNodeStatsCuda, pointsapssNodeStatsSize * sizeof(apssNodeStatsGPU)));
	gpuErrchk( cudaMallocManaged(&m_kernelCuda, sizeof(KernelGPU)));

	gpuErrchk( cudaPeekAtLastError() );
	//Fill GPU memory
	m_kernelCuda->copy(kernel);
	gpuErrchk( cudaPeekAtLastError() );
	std::cout << pointsapssNodeStatsSize << std::endl;
	for (unsigned int i=0; i<pointsapssNodeStatsSize; i++)
	m_apssNodeStatsCuda[i].copy(m_apssOctree->getPointsApssNodeStats()[i]);
	std::cout << "    Done" << std::endl;

	gpuErrchk( cudaPeekAtLastError() );
}

void GlobalAPSS::eraseFromGPU()
{
	gpuErrchk(cudaFree(m_apssNodeStatsCuda));
	gpuErrchk(cudaFree(m_kernelCuda));
	gpuErrchk(cudaFree(m_treeStats));
}

void GlobalAPSS::copyPointsToGPU(unsigned int nbOfVectors, glm::vec3 *outputPoints)
{
	APSS::copyPointsToGPU(nbOfVectors, outputPoints);
	gpuErrchk(cudaMallocManaged(&m_treeStats, nbOfVectors * sizeof (apssStatsGPU)));
}

void GlobalAPSS::erasePointsFromGPU()
{
	APSS::erasePointsFromGPU();
	if (m_treeStats)
	{
	gpuErrchk(cudaFree(m_treeStats));
	m_treeStats = NULL;
	}
}

void GlobalAPSS::projectCPU(const glm::vec3 &qStart, glm::vec3 &outputPoint, glm::vec3 &outputNormal, unsigned int n_iterations, const Kernel &kernel, float stepMaxSize) const
{
	m_apssOctree->projectCPU(qStart, outputPoint, outputNormal, n_iterations, false, kernel, stepMaxSize);
}

void GlobalAPSS::cpuSphere(glm::vec3 const & qStart, float &u0, glm::vec3 &u123, float &u4, Kernel const & kernel) const
{
	m_apssOctree->cpuSphere(qStart, u0, u123, u4, kernel);
}
