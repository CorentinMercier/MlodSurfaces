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

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/scan.h>
#include <thrust/reduce.h>
#include <thrust/fill.h>
#include <thrust/system_error.h>
#include <cuda_profiler_api.h>
#include <cuda_gl_interop.h>

//////////////////////////////////////////////////////////////////////////
///
/// General GPU functions
///
//////////////////////////////////////////////////////////////////////////

#define gpuErrchk(ans) { gpuAssert((ans), #ans, __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char* command, const char *file, int line, bool abort=true)
{
	//printf("%s = %i\n", command, (int)code);
	if (code != cudaSuccess)
	{
		fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}

//basic maths device functions to use in cuda kernels //////////////////////////////////////////////

typedef struct{
	float q[9];
}mat33;

inline __host__ __device__ void setZero(mat33 &m){
	for (int i = 0; i < 9; i++){
		m.q[i] = 0.0;
	}
}

inline __host__ __device__ mat33 setMat33(float m00, float m01, float m02, float m10, float m11, float m12, float m20, float m21, float m22){
	mat33 res;
	res.q[0] = m00;
	res.q[1] = m01;
	res.q[2] = m02;
	res.q[3] = m10;
	res.q[4] = m11;
	res.q[5] = m12;
	res.q[6] = m20;
	res.q[7] = m21;
	res.q[8] = m22;
	return res;
}

inline __host__ __device__ mat33 identity33(){
	mat33 id = setMat33(1.0, 0.0, 0.0,
						0.0, 1.0, 0.0,
						0.0, 0.0, 1.0);
	return id;
}

inline __host__ __device__ mat33 getTranspose(mat33 m){
	mat33 res;
	res = setMat33(m.q[0], m.q[3], m.q[6], m.q[1], m.q[4], m.q[7], m.q[2], m.q[5], m.q[8]);
	return res;
}

inline __host__ __device__ float3 operator*(mat33 m, float3 p){
	float3 res = make_float3(m.q[0] * p.x + m.q[1] * p.y + m.q[2] * p.z,
			m.q[3] * p.x + m.q[4] * p.y + m.q[5] * p.z,
			m.q[6] * p.x + m.q[7] * p.y + m.q[8] * p.z);
	return res;
}

inline __host__ __device__ mat33 operator*(float f, mat33 m){
	mat33 res = setMat33(f * m.q[0], f * m.q[1], f * m.q[2],
			f * m.q[3], f * m.q[4], f * m.q[5],
			f * m.q[6], f * m.q[7], f * m.q[8]);
	return res;
}

inline __host__ __device__ mat33 operator+(mat33 m1, mat33 m2){
	mat33 res = setMat33(m1.q[0] + m2.q[0], m1.q[1] + m2.q[1], m1.q[2] + m2.q[2],
			m1.q[3] + m2.q[3], m1.q[4] + m2.q[4], m1.q[5] + m2.q[5],
			m1.q[6] + m2.q[6], m1.q[7] + m2.q[7], m1.q[8] + m2.q[8]);
	return res;
}

inline __host__ __device__ mat33 operator*(mat33 m1, mat33 m2){
	mat33 res = setMat33(m1.q[0] * m2.q[0] + m1.q[1] * m2.q[3] + m1.q[2] * m2.q[6],
			m1.q[0] * m2.q[1] + m1.q[1] * m2.q[4] + m1.q[2] * m2.q[7],
			m1.q[0] * m2.q[2] + m1.q[1] * m2.q[5] + m1.q[2] * m2.q[8],
			m1.q[3] * m2.q[0] + m1.q[4] * m2.q[3] + m1.q[5] * m2.q[6],
			m1.q[3] * m2.q[1] + m1.q[4] * m2.q[4] + m1.q[5] * m2.q[7],
			m1.q[3] * m2.q[2] + m1.q[4] * m2.q[5] + m1.q[5] * m2.q[8],
			m1.q[6] * m2.q[0] + m1.q[7] * m2.q[3] + m1.q[8] * m2.q[6],
			m1.q[6] * m2.q[1] + m1.q[7] * m2.q[4] + m1.q[8] * m2.q[7],
			m1.q[6] * m2.q[2] + m1.q[7] * m2.q[5] + m1.q[8] * m2.q[8]);
	return res;
}

inline __host__ __device__ float3 operator-( float3 a){
	return make_float3(-a.x, -a.y, -a.z);
}

inline __host__ __device__ float3 operator-( float3 a, float3 b){
	return make_float3(a.x-b.x, a.y-b.y, a.z-b.z);
}

inline __host__ __device__ double3 operator-( double3 a, double3 b){
	return make_double3(a.x-b.x, a.y-b.y, a.z-b.z);
}

inline __host__ __device__ float3 operator+( float3 a, float3 b){
	return make_float3(a.x+b.x, a.y+b.y, a.z+b.z);
}

inline __host__ __device__ double3 operator+( double3 a, double3 b){
	return make_double3(a.x+b.x, a.y+b.y, a.z+b.z);
}

inline __host__ __device__ float3 operator*( float3 a, float3 b){
	return make_float3(a.x*b.x, a.y*b.y, a.z*b.z);
}

inline __host__ __device__ float3 operator/( float3 a, float b){
	return make_float3(a.x/b, a.y/b, a.z/b);
}

inline __host__ __device__ double3 operator/( double3 a, double b){
	return make_double3(a.x/b, a.y/b, a.z/b);
}

inline __host__ __device__ float3 operator*( float a, float3 b){
	return make_float3(a*b.x, a*b.y, a*b.z);
}


inline __device__ double3 cross(double3 a, double3 b){
	return make_double3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

inline __device__ double dot(double3 a, double3 b){
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline __device__ double norme(double3 x){
	return sqrt(x.x * x.x + x.y * x.y + x.z * x.z);
}

inline __device__ double squaredNorm(double3 x){
	return (x.x * x.x + x.y * x.y + x.z * x.z);
}

inline __device__ float3 cross(float3 a, float3 b){
	return make_float3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

inline __device__ float dot(float3 a, float3 b){
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline __device__ float norme(float3 const & a){
	return __fsqrt_rn(a.x * a.x + a.y * a.y + a.z * a.z);
}

inline __device__ float squaredNorm(float3 x){
	return (x.x * x.x + x.y * x.y + x.z * x.z);
}

inline __device__ float squaredDist(float3 x1, float3 x2){
	return ((x1.x - x2.x) * (x1.x - x2.x) + (x1.y - x2.y) * (x1.y - x2.y) + (x1.z - x2.z) * (x1.z - x2.z));
}

inline __device__ float dist(float3 x1, float3 x2){
	return __fsqrt_rn(squaredDist(x1, x2));
}

inline __device__ float3 normalize(float3 x){
	float norm = __fsqrt_rn(x.x * x.x + x.y * x.y + x.z * x.z);
	if(norm!=0){
		return make_float3(x.x / norm, x.y / norm, x.z / norm);
	}else{
		return make_float3(0.0, 0.0, 0.0);
	}
}

inline __device__ double3 normalize(double3 x){
	double norm = sqrt(x.x * x.x + x.y * x.y + x.z * x.z);
	if(norm!=0){
		return make_double3(x.x / norm, x.y / norm, x.z / norm);
	}else{
		return make_double3(0.0, 0.0, 0.0);
	}
}

inline float3 glmToFloat3(glm::vec3 a)
{
	return make_float3(a.x, a.y, a.z);
}

//End of basic math

inline dim3 computeNbBlocks(unsigned int nbThreads, unsigned int nbThreadsPerBlock){
	dim3 nbBlocks = ceil((float)(nbThreads)/(float)nbThreadsPerBlock);
	if (nbBlocks.x > 65535){
		nbBlocks.y = ceil((float)nbBlocks.x / (float)65535);
		nbBlocks.x = 65535;
	}
	return nbBlocks;
}

//////////////////////////////////////////////////////////////////////////
///
/// GPU Kernel
///
//////////////////////////////////////////////////////////////////////////

struct KernelGPU{
	enum KernelType { GAUSSIAN , WENDLAND , SINGULAR , GAUSSIAN_MULTIPLE };
	float sigma;
	float exponent;
	float constant_shift;
	float *sigmas;
	unsigned int sigmasSize;

	KernelType type;

	__device__ float w( float3 const & p , float3 const & eta) const {
		if( type == KernelGPU::GAUSSIAN_MULTIPLE ) {
			float sum_g = 0.f;//constant_shift;
			for( unsigned int i = 0 ; i < sigmasSize ; ++i )
				sum_g = sum_g + __fdiv_rn( __expf( __fdiv_rn( - squaredNorm(p-eta) , (2.f * sigmas[i]*sigmas[i]) ) ) , (sigmas[i]*sigmas[i]*sigmas[i]) )  ;
			return sum_g;
		}
		if( type == KernelGPU::SINGULAR )
			return __fdiv_rn(__powf(max(constant_shift, 0.0001), exponent), __powf(__powf(norme(p - eta), 2) + constant_shift, exponent));
		if( type == KernelGPU::GAUSSIAN )
		{
			return __expf( __fdiv_rn(- squaredNorm(p-eta) , (2.0 * sigma*sigma)) ) + constant_shift;
		}
		return 1.f;
	}
	__host__ void copy(Kernel const & a)
	{
		if(a.type == Kernel::GAUSSIAN)
			type = KernelGPU::GAUSSIAN;
		if(a.type == Kernel::GAUSSIAN_MULTIPLE)
			type = KernelGPU::GAUSSIAN_MULTIPLE;
		if(a.type == Kernel::SINGULAR)
			type = KernelGPU::SINGULAR;
		sigma = a.sigma;
		exponent = a.exponent;
		constant_shift = a.constant_shift;
		sigmasSize = a.sigmas.size();
		if (sigmasSize !=0)
		{
			gpuErrchk( cudaMalloc(&sigmas, a.sigmas.size() * sizeof(float)) );
			gpuErrchk( cudaMemcpy(sigmas, a.sigmas.data(), a.sigmas.size() * sizeof(float), cudaMemcpyHostToDevice) );
		}
	}
	__host__ void print() const
	{
		std::cout << "sigma: " << sigma << std::endl;
		std::cout << "exponent: " << exponent << std::endl;
		std::cout << "constant_shift: " << constant_shift << std::endl;
		std::cout << "sigmasSize: " << sigmasSize << std::endl;
		for (unsigned i=0; i<sigmasSize; i++)
			std::cout << i << " " << sigmas[i] << std::endl;
	}
};

struct BBOXGPU
{
	float3 bb,BB;

	__device__ inline void squareDiagonal(float & sd) const { sd = squaredNorm(BB - bb); }
	__device__ inline void diagonal(float & d) const { d = sqrt( (float)squaredNorm(BB - bb) ); }
	__device__ inline void radius(float & r) const {float d; diagonal(d); r = d / 2.0; }
	__device__ inline void center(float3 & c) const { c = (bb + BB) / 2.0; }
	__device__ inline void squareRadius(float & sr) const { float sd; squareDiagonal(sd); sr = sd / 4.0; }
	__host__ void copy(BBOX const & a)
	{
		bb = make_float3(a.bb.x, a.bb.y, a.bb.z);
		BB = make_float3(a.BB.x, a.BB.y, a.BB.z);
	}
};

struct octreeNodeGPU2 {
	unsigned int depth;
	float3 center;
	float radius;
	float3 s_ai_pi , s_ai_ni;
	float s_ai_pi_ni , s_ai_pi_pi;
	float s_ai;
	unsigned int indicesSize;
	unsigned int* indices = nullptr;
	unsigned int numberOfChildren;
	int firstChild;
	int nextBrother;
	int father;
};

struct apssStatsGPU{
	float3 s_ai_wi_pi , s_ai_wi_ni;
	float s_ai_wi, s_ai_wi_pi_ni , s_ai_wi_pi_pi;

	__device__ void init()
	{
		s_ai_wi = 0;
		s_ai_wi_ni = make_float3(0.0, 0.0, 0.0);
		s_ai_wi_pi = make_float3(0.0, 0.0, 0.0);
		s_ai_wi_pi_ni = 0;
		s_ai_wi_pi_pi = 0;
	}

	__device__ void init(float w, octreeNodeGPU2 * const node)
	{
		s_ai_wi = w * node->s_ai;
		s_ai_wi_ni = w * node->s_ai_ni;
		s_ai_wi_pi = w * node->s_ai_pi;
		s_ai_wi_pi_ni = w * node->s_ai_pi_ni;
		s_ai_wi_pi_pi = w * node->s_ai_pi_pi;
	}

	__device__ void init(float w, octreeNodeGPU2 const node)
	{
		s_ai_wi = w * node.s_ai;
		s_ai_wi_ni = w * node.s_ai_ni;
		s_ai_wi_pi = w * node.s_ai_pi;
		s_ai_wi_pi_ni = w * node.s_ai_pi_ni;
		s_ai_wi_pi_pi = w * node.s_ai_pi_pi;
	}

	__device__ void operator += (apssStatsGPU const * o) {
		s_ai_wi_pi = s_ai_wi_pi + o->s_ai_wi_pi;
		s_ai_wi_ni = s_ai_wi_ni + o->s_ai_wi_ni;
		s_ai_wi += o->s_ai_wi;
		s_ai_wi_pi_ni += o->s_ai_wi_pi_ni;
		s_ai_wi_pi_pi += o->s_ai_wi_pi_pi;
	}
	__device__ void add (apssStatsGPU const *i1){
		this->s_ai_wi_pi = i1->s_ai_wi_pi + s_ai_wi_pi;
		this->s_ai_wi_ni = i1->s_ai_wi_ni + s_ai_wi_ni;
		this->s_ai_wi = i1->s_ai_wi + s_ai_wi;
		this->s_ai_wi_pi_ni = i1->s_ai_wi_pi_ni + s_ai_wi_pi_ni;
		this->s_ai_wi_pi_pi = i1->s_ai_wi_pi_pi + s_ai_wi_pi_pi;
	}
	__device__ void add (apssStatsGPU const &i1){
		s_ai_wi_pi = i1.s_ai_wi_pi + s_ai_wi_pi;
		s_ai_wi_ni = i1.s_ai_wi_ni + s_ai_wi_ni;
		s_ai_wi = i1.s_ai_wi + s_ai_wi;
		s_ai_wi_pi_ni = i1.s_ai_wi_pi_ni + s_ai_wi_pi_ni;
		s_ai_wi_pi_pi = i1.s_ai_wi_pi_pi + s_ai_wi_pi_pi;
	}
	__device__ void times (apssStatsGPU * const i, float const f){
		this->s_ai_wi_pi = f * i->s_ai_wi_pi;
		this->s_ai_wi_ni = f * i->s_ai_wi_ni;
		this->s_ai_wi = i->s_ai_wi * f;
		this->s_ai_wi_pi_ni = i->s_ai_wi_pi_ni * f;
		this->s_ai_wi_pi_pi = i->s_ai_wi_pi_pi * f;
	}
	__device__ void times (apssStatsGPU const & i, float const f){
		s_ai_wi_pi = f * i.s_ai_wi_pi;
		s_ai_wi_ni = f * i.s_ai_wi_ni;
		s_ai_wi = i.s_ai_wi * f;
		s_ai_wi_pi_ni = i.s_ai_wi_pi_ni * f;
		s_ai_wi_pi_pi = i.s_ai_wi_pi_pi * f;
	}
	__host__ void copy(apssStats const & a)
	{
		s_ai_wi = a.s_ai_wi;
		s_ai_wi_ni = make_float3(a.s_ai_wi_ni.x, a.s_ai_wi_ni.y, a.s_ai_wi_ni.z);
		s_ai_wi_pi = make_float3(a.s_ai_wi_pi.x, a.s_ai_wi_pi.y, a.s_ai_wi_pi.z);
		s_ai_wi_pi_ni = a.s_ai_wi_pi_ni;
		s_ai_wi_pi_pi = a.s_ai_wi_pi_pi;
	}
};

void APSS::stop()
{
	gpuErrchk( cudaProfilerStop() );
}

void APSS::eraseAndRestart()
{
	eraseFromGPU();
	delete(m_apssOctree);
	m_apssOctree = new APSSOctree();
}

void APSS::copyPointsToGPU(unsigned int nbOfVectors, glm::vec3 * outputPoints)
{
	gpuErrchk( cudaMalloc(&m_outPts,  nbOfVectors * sizeof(float3)));
	gpuErrchk( cudaMalloc(&m_outNmls, nbOfVectors * sizeof(float3)));
	gpuErrchk( cudaMemcpy(m_outPts, outputPoints, nbOfVectors * sizeof (float3), cudaMemcpyHostToDevice));
}

__global__
void reorganizePts(unsigned n, unsigned const * invalidPts, unsigned endValid, float3 * outputPoint)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int stride = blockDim.x * gridDim.x;
	for (int i = index; i < n; i+=stride)
	{
		float3 temp = outputPoint[invalidPts[i]];
		outputPoint[invalidPts[i]] = outputPoint[endValid - 1 - i];
		outputPoint[endValid - 1 - i] = temp;
	}
}

void APSS::printApss() const
{
	std::cout << "m_numberOfNodes: " << m_numberOfNodes << std::endl;
	std::cout << "m_scalingProtectionSphere: " << m_scalingProtectionSphere << std::endl;
	std::cout << "m_minDepth: " << m_minDepth << std::endl;
	std::cout << "KernelCuda:" << std::endl;
	m_kernelCuda->print();
	m_apssOctree->printState();
}

void APSS::reorganizePoints(const std::vector<unsigned>* invalidPts, unsigned endValid)
{
	gpuErrchk( cudaPeekAtLastError() );
	unsigned nbOfElmts = invalidPts->size();
	unsigned * invalid = NULL;
	gpuErrchk( cudaMalloc(&invalid, nbOfElmts * sizeof(unsigned)));
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaMemcpy(invalid, invalidPts->data(), nbOfElmts * sizeof (int), cudaMemcpyHostToDevice));
	gpuErrchk( cudaPeekAtLastError() );
	int blockSize = 128;
	dim3 numBlocks = computeNbBlocks(nbOfElmts, blockSize);
	reorganizePts<<<numBlocks, blockSize>>>(nbOfElmts, invalid, endValid, m_outPts);
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk(cudaFree(invalid));
}

void APSS::erasePointsFromGPU()
{
	if (m_outPts)
	{
		gpuErrchk(cudaFree(m_outPts));
		m_outPts = NULL;
	}
	if (m_outNmls)
	{
		gpuErrchk(cudaFree(m_outNmls));
		m_outNmls = NULL;
	}
}

void APSS::updateKernel(Kernel const & kernel)
{
	KernelGPU kernelOnCPU;
	kernelOnCPU.copy(kernel);
	gpuErrchk(cudaMalloc(&m_kernelCuda, sizeof(KernelGPU)));
	gpuErrchk(cudaMemcpy(m_kernelCuda, &kernelOnCPU, sizeof(KernelGPU), cudaMemcpyHostToDevice));
}
