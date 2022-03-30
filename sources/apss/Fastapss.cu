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

#include "Fastapss.h"
#include "GPUUtils.cu"
#include "GlobalAPSS.cu"

Fastapss::Fastapss(APSSOctree *apssOctree)
{
	m_apssOctree = apssOctree;
	std::cout << "FastAPSS method used" << std::endl;
}

Fastapss::Fastapss(APSSOctree *apssOctree, const Kernel &kernel, unsigned int nbOfVectors):Fastapss(apssOctree)
{
	copyToGPU(kernel, nbOfVectors);
}


//////////////////////////////////////////////////////////////////////////
///
/// GPU structs declarations and copies
///
//////////////////////////////////////////////////////////////////////////

struct apssNodeStatsGPUSoA{
	unsigned int m_nbOfLeaves;
	float3 * s_ai_pi = NULL;
	float3* s_ai_ni = NULL;
	float * s_ai_pi_ni = NULL;
	float * s_ai_pi_pi = NULL;
	float * s_ai = NULL;

	__host__ void init(unsigned int nbOfLeaves)
	{
		m_nbOfLeaves = nbOfLeaves;
		cudaMallocManaged(&s_ai_pi, m_nbOfLeaves * sizeof (float3));
		cudaMallocManaged(&s_ai_ni, m_nbOfLeaves * sizeof (float3));
		cudaMallocManaged(&s_ai_pi_ni, m_nbOfLeaves * sizeof (float));
		cudaMallocManaged(&s_ai_pi_pi, m_nbOfLeaves * sizeof (float));
		cudaMallocManaged(&s_ai, m_nbOfLeaves * sizeof (float));
	}
	__host__ void erase()
	{
		if (s_ai_pi){ gpuErrchk(cudaFree(s_ai_pi)); s_ai_pi = NULL;}
		if (s_ai_ni){ gpuErrchk(cudaFree(s_ai_ni)); s_ai_ni = NULL;}
		if (s_ai_pi_ni){ gpuErrchk(cudaFree(s_ai_pi_ni)); s_ai_pi_ni = NULL;}
		if (s_ai_pi_pi){ gpuErrchk(cudaFree(s_ai_pi_pi)); s_ai_pi_pi = NULL;}
		if (s_ai){ gpuErrchk(cudaFree(s_ai)); s_ai = NULL;}
	}
	__host__ void copy(std::vector<apssNodeStats> const & a)
	{
		init(a.size());
		for (unsigned int i=0; i<a.size(); i++)
		{
			s_ai_ni[i] = make_float3(a[i].s_ai_ni.x, a[i].s_ai_ni.y, a[i].s_ai_ni.z);
			s_ai_pi[i] = make_float3(a[i].s_ai_pi.x, a[i].s_ai_pi.y, a[i].s_ai_pi.z);
			s_ai_pi_ni[i] = a[i].s_ai_pi_ni;
			s_ai_pi_pi[i] = a[i].s_ai_pi_pi;
			s_ai[i] = a[i].s_ai;
		}
	}
};

__device__ void func(float x, float & result)
{
#ifdef EXPFUNC
	if (x>=1){
		result = 1;
	}
	else
		result = __expf(__fdiv_rn(-(__expf(__fdiv_rn(1, x-1))), x*x));
#else
	result = (1 - cos(3.14159265 * x))/2.0;
#endif
}

void copy(octreeNode const * a, octreeNodeGPU2 * b, unsigned int numberOfNodes)
{
	std::vector<const octreeNode*> nodes;
	std::vector<relations> related;
	nodes.resize(numberOfNodes, nullptr);
	related.resize(numberOfNodes);
	a->getRelations(nodes, related);
	for(unsigned int i=0; i<numberOfNodes; i++)
	{
		b[i].depth = nodes[i]->depth;
		b[i].radius = float(nodes[i]->getBoundingBox().radius());
		b[i].center = glmToFloat3(nodes[i]->getBoundingBox().center());
		b[i].s_ai_pi = glmToFloat3(nodes[i]->nodeapssNodeStats.s_ai_pi);
		b[i].s_ai_ni = glmToFloat3(nodes[i]->nodeapssNodeStats.s_ai_ni);
		b[i].s_ai_pi_ni = nodes[i]->nodeapssNodeStats.s_ai_pi_ni;
		b[i].s_ai_pi_pi = nodes[i]->nodeapssNodeStats.s_ai_pi_pi;
		b[i].s_ai = nodes[i]->nodeapssNodeStats.s_ai;
		b[i].indicesSize = nodes[i]->indices.size();
		if (b[i].indicesSize != 0)
		{
			cudaMallocManaged(&b[i].indices, b[i].indicesSize * sizeof (unsigned int));
			for (unsigned int j=0; j<b[i].indicesSize; j++)
			{
				b[i].indices[j] = nodes[i]->indices[j];
			}
		}
		b[i].numberOfChildren = nodes[i]->numberOfChildren;
		b[i].father = related[i].father;
		b[i].nextBrother = related[i].nextBrother;
		b[i].firstChild = related[i].firstChild;
	}
}


//////////////////////////////////////////////////////////////////////////
///
/// GPU Version of FastAPSS -> Non recursive version (use of while) with Array of Structs
///
//////////////////////////////////////////////////////////////////////////


__device__
#ifdef SHARED_MEM
void approximateAPSSNonRecursive( float3 const & q , KernelGPU const * kernel , unsigned int const & minimal_depth , float const & scalingProtectionSphere ,
								  apssNodeStatsGPUSoA const & leavesapssNodeStats,
							 const octreeNodeGPU2 *__restrict  node,
							 apssStatsGPU & treeStats, const octreeNodeGPU2 * nodeInSharedMemory, const unsigned sharedMemorySize)
#else
void approximateAPSSNonRecursive( float3 const & q , KernelGPU const * kernel , unsigned int const & minimal_depth , float const & scalingProtectionSphere ,
								  apssNodeStatsGPUSoA const & leavesapssNodeStats,
							 const octreeNodeGPU2 *__restrict  node,
							 apssStatsGPU & treeStats )
#endif
{
	apssStatsGPU depthStats[MAX_DEPTH];
	int currentChild[MAX_DEPTH - 1];
	for (unsigned int i=0; i<MAX_DEPTH - 1; i++)
	{
		depthStats[i].init();
		currentChild[i] = -1;
	}
	depthStats[MAX_DEPTH - 1].init();
	bool not_finished = true;
#ifdef SHARED_MEM
	octreeNodeGPU2 currentNode = nodeInSharedMemory[0];
#else
	octreeNodeGPU2 currentNode = node[0];
#endif
	while(not_finished)
	{
		if( currentNode.indicesSize > 0 ) { //Leaf
#pragma unroll 4
			for( int lIt = 0 ; lIt < currentNode.indicesSize ; ++lIt ) {
				unsigned int leafIdx = currentNode.indices[lIt];
				float w = kernel->w(leavesapssNodeStats.s_ai_pi[leafIdx] / leavesapssNodeStats.s_ai[leafIdx], q);

				depthStats[currentNode.depth].s_ai_wi = __fmaf_rn(w, leavesapssNodeStats.s_ai[leafIdx], depthStats[currentNode.depth].s_ai_wi);
				depthStats[currentNode.depth].s_ai_wi_ni = depthStats[currentNode.depth].s_ai_wi_ni + w * leavesapssNodeStats.s_ai_ni[leafIdx];
				depthStats[currentNode.depth].s_ai_wi_pi = depthStats[currentNode.depth].s_ai_wi_pi + w * leavesapssNodeStats.s_ai_pi[leafIdx];
				depthStats[currentNode.depth].s_ai_wi_pi_ni = __fmaf_rn(w, leavesapssNodeStats.s_ai_pi_ni[leafIdx], depthStats[currentNode.depth].s_ai_wi_pi_ni);
				depthStats[currentNode.depth].s_ai_wi_pi_pi = __fmaf_rn(w, leavesapssNodeStats.s_ai_pi_pi[leafIdx], depthStats[currentNode.depth].s_ai_wi_pi_pi);
			}
			int nextNode = currentNode.father;
			if (nextNode != -1)
			{
#ifdef SHARED_MEM
				if (nextNode < sharedMemorySize)
					currentNode = nodeInSharedMemory[nextNode];
				else
#endif
					currentNode = node[nextNode];
			}
			else
				not_finished = false;
		}
		else if (currentNode.depth >= minimal_depth) {
			float3 p = currentNode.s_ai_pi / currentNode.s_ai;
			if( norme(q - currentNode.center)  > scalingProtectionSphere * currentNode.radius )
			{
				float w = kernel->w( p , q );
				depthStats[currentNode.depth].s_ai_wi = __fmaf_rn(w, currentNode.s_ai, depthStats[currentNode.depth].s_ai_wi);
				depthStats[currentNode.depth].s_ai_wi_ni = depthStats[currentNode.depth].s_ai_wi_ni + w * currentNode.s_ai_ni;
				depthStats[currentNode.depth].s_ai_wi_pi = depthStats[currentNode.depth].s_ai_wi_pi + w * currentNode.s_ai_pi;
				depthStats[currentNode.depth].s_ai_wi_pi_ni = __fmaf_rn(w, currentNode.s_ai_pi_ni, depthStats[currentNode.depth].s_ai_wi_pi_ni);
				depthStats[currentNode.depth].s_ai_wi_pi_pi = __fmaf_rn(w, currentNode.s_ai_pi_pi, depthStats[currentNode.depth].s_ai_wi_pi_pi);

				int nextNode = currentNode.father;
				if (nextNode != -1)
				{
#ifdef SHARED_MEM
					if (nextNode < sharedMemorySize)
						currentNode = nodeInSharedMemory[nextNode];
					else
#endif
						currentNode = node[nextNode];
				}
				else
					not_finished = false;
			}
			else {
				// do blending between levels

				if (currentChild[currentNode.depth] != -1) //Computation done for the child currentChild
				{
					apssStatsGPU parentStats;
					float w = kernel->w( p , q ) * node[currentChild[currentNode.depth]].s_ai / currentNode.s_ai ;
					parentStats.init(w, currentNode);

					float distToParent = fabs(norme(q - currentNode.center) - (currentNode.radius * scalingProtectionSphere));

					float distChild = max(norme(q - node[currentChild[currentNode.depth]].center) - (node[currentChild[currentNode.depth]].radius * scalingProtectionSphere),0.f);
					float wChild = distToParent / (distToParent + distChild);
					float valChild;
					func(wChild, valChild);
					float valParent = 1 - valChild;

					depthStats[currentNode.depth+1].times(depthStats[currentNode.depth+1], valChild);
					depthStats[currentNode.depth].add(depthStats[currentNode.depth+1]);
					parentStats.times(parentStats, valParent);
					depthStats[currentNode.depth].add(parentStats);

					//Reinitialize this level
					depthStats[currentNode.depth+1].init();

					int nextNode = node[currentChild[currentNode.depth]].nextBrother;
					if (nextNode != -1)
					{
						currentChild[currentNode.depth] = nextNode;
#ifdef SHARED_MEM
						if (nextNode < sharedMemorySize)
							currentNode = nodeInSharedMemory[nextNode];
						else
#endif
							currentNode = node[nextNode];
					}
					else {
						currentChild[currentNode.depth] = -1;
						nextNode = currentNode.father;
						if (nextNode != -1)
						{
#ifdef SHARED_MEM
							if (nextNode < sharedMemorySize)
								currentNode = nodeInSharedMemory[nextNode];
							else
#endif
								currentNode = node[nextNode];
						}
						else
							not_finished = false;
					}
				}
				else {
					currentChild[currentNode.depth] = currentNode.firstChild;
#ifdef SHARED_MEM
					if (currentNode.firstChild < sharedMemorySize)
						currentNode = nodeInSharedMemory[currentNode.firstChild];
					else
#endif
						currentNode = node[currentNode.firstChild];
				}
			}
		}

		// if not, then just accumulate the wn of the children:
		else
		{
			if (currentChild[currentNode.depth] != -1) //Computation done for the child currentChild
			{
				depthStats[currentNode.depth].add(depthStats[currentNode.depth+1]);
				//Reinitialize this level
				depthStats[currentNode.depth+1].init();
				int nextNode = node[currentChild[currentNode.depth]].nextBrother;
				if (nextNode != -1)
				{
					currentChild[currentNode.depth] = nextNode;
					currentNode = node[nextNode];
				}
				else {
					currentChild[currentNode.depth] = -1;
					nextNode = currentNode.father;
					if (nextNode != -1)
						currentNode = node[nextNode];
					else
						not_finished = false;
				}
			}
			else {
				currentChild[currentNode.depth] = currentNode.firstChild;
				currentNode = node[currentNode.firstChild];
			}
		}

	}
	treeStats.add(depthStats[0]);
}

__device__
#ifdef SHARED_MEM
void oneStepProjection( float3 const & q , float3 & outputPoint , float3 & outputNormal ,
						  const KernelGPU * kernel, unsigned int minimal_depth , float scalingProtectionSphere,
						  apssNodeStatsGPUSoA const & leavesapssNodeStats, const octreeNodeGPU2 *__restrict node, float stepMaxSize, const octreeNodeGPU2 * nodeInSharedMemory, const unsigned sharedMemorySize)
#else
void oneStepProjection( float3 const & q , float3 & outputPoint , float3 & outputNormal ,
								  KernelGPU const * kernel, unsigned int minimal_depth , float scalingProtectionSphere,
						apssNodeStatsGPUSoA const & leavesapssNodeStats, const octreeNodeGPU2 *__restrict node, float stepMaxSize)
#endif
{
	apssStatsGPU treeStats;
	treeStats.init();
#ifdef SHARED_MEM
	approximateAPSSNonRecursive( q , kernel , minimal_depth , scalingProtectionSphere , leavesapssNodeStats, node, treeStats, nodeInSharedMemory, sharedMemorySize);
#else
	approximateAPSSNonRecursive( q , kernel , minimal_depth , scalingProtectionSphere , leavesapssNodeStats, node, treeStats);
#endif

	float u4 = 0.5*( treeStats.s_ai_wi_pi_ni/*/treeStats.s_ai_wi*/ - dot((treeStats.s_ai_wi_pi/*/treeStats.s_ai_wi*/),(treeStats.s_ai_wi_ni/treeStats.s_ai_wi)) ) /
			( treeStats.s_ai_wi_pi_pi/*/treeStats.s_ai_wi*/ - dot((treeStats.s_ai_wi_pi/*/treeStats.s_ai_wi*/),(treeStats.s_ai_wi_pi/treeStats.s_ai_wi)) );
	float3 u123= (treeStats.s_ai_wi_ni - 2*u4*treeStats.s_ai_wi_pi)/treeStats.s_ai_wi;
	float u0= -(dot(treeStats.s_ai_wi_pi,u123) + u4*treeStats.s_ai_wi_pi_pi)/treeStats.s_ai_wi;

	float3 inputPt = q;
	outputPoint = q;

	if( fabs(u4) < 0.000001 ) {
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

	//Invalid points have a normal equals to (0, 0, 0) and their position is kept for debug
	if (isnan(outputPoint.x) || isnan(outputPoint.y) || isnan(outputPoint.z))
	{
		outputPoint = inputPt;
		outputNormal = make_float3(0, 0, 0);
		return;
	}

	float stepSize = norme(outputPoint - inputPt);
	if (stepSize > stepMaxSize)
		outputPoint = inputPt + (stepMaxSize / stepSize) * (outputPoint - inputPt);
}

__global__
void projectCuda(unsigned int n, float3 const * qStart , float3 * outputPoint , float3 * outputNormal , unsigned int n_iterations ,
				 KernelGPU const * kernel, unsigned int minimal_depth, float scalingProtectionSphere,
				 apssNodeStatsGPUSoA const & leavesapssNodeStats, const octreeNodeGPU2 *__restrict node, float stepMaxSize, unsigned blockSize)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int stride = blockDim.x * gridDim.x;

#ifdef SHARED_MEM
	const unsigned sharedMemorySize = 73;
	__shared__ octreeNodeGPU2 nodeInSharedMemory[sharedMemorySize];
	for (int i = threadIdx.x; i < sharedMemorySize; i+=blockSize)
	{
		nodeInSharedMemory[i] = node[i];
	}
	__syncthreads();
	__shared__ KernelGPU k;
	if (threadIdx.x == 0)
		k = *kernel;
	__syncthreads();
	for (int i = threadIdx.x; i < kernel->sigmasSize; i+=blockSize)
	{
		k.sigmas[i] = kernel->sigmas[i];
	}
	__syncthreads();
#endif

	for (int i = index; i < n; i+=stride)
	{
		outputPoint[i] = qStart[i];
		for (unsigned int j=0; j<n_iterations; j++)
		{
#ifdef SHARED_MEM
			oneStepProjection(outputPoint[i], outputPoint[i], outputNormal[i], &k, minimal_depth, scalingProtectionSphere, leavesapssNodeStats, node, stepMaxSize, nodeInSharedMemory, sharedMemorySize);
#else
			oneStepProjection(outputPoint[i], outputPoint[i], outputNormal[i], kernel, minimal_depth, scalingProtectionSphere, leavesapssNodeStats, node, stepMaxSize);
#endif
		}
	}
}

void Fastapss::project(unsigned int nbOfVectors, glm::vec3 *outputPoints, glm::vec3 *outputNormals, unsigned int n_iterations, const Kernel &kernel, float stepMaxSize) const
{
	if (nbOfVectors == 0)
	{
		std::cerr << "No points to project" << std::endl;
		return;
	}
	gpuErrchk( cudaPeekAtLastError() );
	int blockSize = 128;
	dim3 numBlocks = computeNbBlocks(nbOfVectors, blockSize);
	if (stepMaxSize < 0)
		stepMaxSize =  m_apssOctree->getBoundingBox().radius();
	projectCuda<<<numBlocks, blockSize>>>(nbOfVectors, m_outPts, m_outPts, m_outNmls, n_iterations, m_kernelCuda, m_minDepth, m_scalingProtectionSphere, *m_apssNodeStatsCudaSoA, m_nodesBis, stepMaxSize, blockSize);
	gpuErrchk( cudaPeekAtLastError() );

	gpuErrchk( cudaMemcpy(outputPoints, m_outPts, nbOfVectors * sizeof (float3), cudaMemcpyDeviceToHost));
	gpuErrchk( cudaMemcpy(outputNormals, m_outNmls, nbOfVectors * sizeof (float3), cudaMemcpyDeviceToHost));
}

void Fastapss::copyToGPU(const Kernel &kernel, unsigned int nbOfVectors)
{
	//Allocate memory on GPU
	gpuErrchk( cudaPeekAtLastError() );
	std::cout << "Tree is being copied to GPU memory... " << std::flush;
	m_numberOfNodes = m_apssOctree->get_root()->nbOfNodes();
	std::cout << "Number of nodes : " << m_numberOfNodes << std::endl;
	gpuErrchk( cudaMallocManaged(&m_nodesBis, m_numberOfNodes * sizeof (octreeNodeGPU2)));
	gpuErrchk( cudaMallocManaged(&m_apssNodeStatsCudaSoA, sizeof(apssNodeStatsGPUSoA)));
	gpuErrchk( cudaMallocManaged(&m_kernelCuda, sizeof(KernelGPU)));

	gpuErrchk( cudaPeekAtLastError() );
	//Fill GPU memory
	copy(m_apssOctree->get_root(), m_nodesBis, m_numberOfNodes);
	gpuErrchk( cudaPeekAtLastError() );
	m_kernelCuda->copy(kernel);
	gpuErrchk( cudaPeekAtLastError() );
	m_apssNodeStatsCudaSoA->copy(m_apssOctree->getPointsApssNodeStats());
	std::cout << "    Done" << std::endl;

	gpuErrchk( cudaPeekAtLastError() );
	std::cout << "Tree depth: " << m_apssOctree->get_root()->computeDepth() << std::endl;
}

void Fastapss::eraseFromGPU()
{
	gpuErrchk( cudaPeekAtLastError() );
	cudaFree(m_kernelCuda);
	m_apssNodeStatsCudaSoA->erase();
	cudaFree(m_apssNodeStatsCudaSoA);
	gpuErrchk( cudaPeekAtLastError() );
	for (unsigned int i=0; i<m_numberOfNodes; i++)
		m_nodesBis[i].erase();
	cudaFree(m_nodesBis);
	gpuErrchk( cudaPeekAtLastError() );
}

void Fastapss::projectCPU(const glm::vec3 &qStart, glm::vec3 &outputPoint, glm::vec3 &outputNormal, unsigned int n_iterations, const Kernel &kernel, float stepMaxSize) const
{
	m_apssOctree->projectCPU(qStart, outputPoint, outputNormal, n_iterations, true, kernel, stepMaxSize);
}
