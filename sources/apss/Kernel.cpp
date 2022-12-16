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

#include "Kernel.h"

Kernel::Kernel() : sigma(2.0) , exponent(1.0) , constant_shift(0) , type(GAUSSIAN)
{
	init();
}

Kernel::Kernel(Kernel  const &newKernel)
{
	sigma = newKernel.sigma;
	type = newKernel.type;
	exponent = newKernel.exponent;
	constant_shift = newKernel.constant_shift;
	if (newKernel.sigmas.size() != 0)
	{
		sigmas.resize(newKernel.sigmas.size());
		for (unsigned int i=0; i<newKernel.sigmas.size(); i++)
			sigmas[i] = newKernel.sigmas[i];
	}
	fixedSigmas.resize(newKernel.fixedSigmas.size());
	for (unsigned i=0; i<fixedSigmas.size(); i++)
		fixedSigmas[i] = newKernel.fixedSigmas[i];
	beginningSigma = newKernel.beginningSigma;
	factor = newKernel.factor;
	nbOfComputedSigmas = newKernel.nbOfComputedSigmas;
	init();
}

void Kernel::init()
{
}

void Kernel::operator =(Kernel const &k)
{
	sigma = k.sigma;
	type = k.type;
	exponent = k.exponent;
	constant_shift = k.constant_shift;
	if (k.sigmas.size() != 0)
	{
		sigmas.resize(k.sigmas.size());
		for (unsigned int i=0; i<k.sigmas.size(); i++)
			sigmas[i] = k.sigmas[i];
	}
	fixedSigmas.resize(k.fixedSigmas.size());
	for (unsigned i=0; i<fixedSigmas.size(); i++)
		fixedSigmas[i] = k.fixedSigmas[i];
	beginningSigma = k.beginningSigma;
	factor = k.factor;
	nbOfComputedSigmas = k.nbOfComputedSigmas;
}

float Kernel::w( glm::vec3 const & p , glm::vec3 const & eta ) const {
	if( type == GAUSSIAN_MULTIPLE ) {
		double sum_g = 0.0;//constant_shift;
		for( unsigned int i = 0 ; i < sigmas.size() ; ++i )
			sum_g += exp( - glm::distance2(p,eta) / (2.0 * sigmas[i]*sigmas[i]) ) / (sigmas[i]*sigmas[i]*sigmas[i])  ;
		if (sum_g == 0)
		{
			std::cerr << "sum_g is 0" << std::endl;
			sum_g = constant_shift;
		}
		return sum_g;
	}
	if( type == GAUSSIAN )
		return exp( - glm::distance2(p,eta) / (2.0 * sigma*sigma) ) + constant_shift;
	if( type == SINGULAR ) {
		const float val = std::max(constant_shift, 0.0001f);
		return std::pow(val, exponent)/(std::pow(glm::distance2(p,eta) + constant_shift, exponent));
	}
	return 1.0;
}

glm::vec4 Kernel::gradW( glm::vec3 const & p , glm::vec3 const & eta ) const {
	float w = constant_shift;
	glm::vec3 gradW(0.f);
	if( type == GAUSSIAN_MULTIPLE ) {
		for( unsigned int i = 0 ; i < sigmas.size() ; ++i )
		{
			w += exp( - glm::distance2(p,eta) / (2.0 * sigmas[i]*sigmas[i]) ) / (sigmas[i]*sigmas[i]*sigmas[i])  ;
			gradW -= (eta - p) / (sigmas[i]*sigmas[i]);
		}
		if (w == 0)
			std::cerr << "w is 0" << std::endl;
		gradW *= w;
	}
	else if( type == GAUSSIAN )
		w = exp( - glm::distance2(p,eta) / (2.0 * sigma*sigma) ) + constant_shift;
	else if( type == SINGULAR ) {
		w = 1.0/pow( std::max<double>( glm::distance(p,eta), 0.000001) / sigma , exponent) + constant_shift;
	}
	return glm::vec4(gradW, w);
}

void Kernel::draw(GLuint program)
{
}

void Kernel::updateSigmas(bool computedSigmas, float radius)
{
	if (computedSigmas)
	{
		sigmas.resize(nbOfComputedSigmas);
		sigmas[0] = beginningSigma;
		for (unsigned i=1; i<nbOfComputedSigmas; i++)
			sigmas[i] = sigmas[i-1] * factor;
	}
	else
	{
		sigmas.resize(fixedSigmas.size());
		for (unsigned i=0; i<sigmas.size(); i++)
			sigmas[i] = fixedSigmas[i];
	}
	sigma *= radius;
	for (unsigned i=0; i<sigmas.size(); i++)
		sigmas[i] *= radius;
}
