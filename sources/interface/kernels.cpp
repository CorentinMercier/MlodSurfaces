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

#include "kernels.h"
#include "ui_kernels.h"

Kernels::Kernels(QWidget *parent) :
	QWidget(parent),
	ui(new Ui::Kernels)
{
	ui->setupUi(this);

	connect(ui->addGaussian, SIGNAL(clicked()), this, SLOT(add_gaussian()));
	connect(ui->removeGaussian, SIGNAL(clicked()), this, SLOT(remove_gaussian()));

	connect(ui->tabWidget, SIGNAL(currentChanged(int)), this, SLOT(changed_tab()));
}

Kernels::~Kernels()
{
	delete ui;
}

void Kernels::changed_tab()
{
	switch(ui->tabWidget->currentIndex())
	{
	case 0:
		singular_tab();
		break;
	case 1:
		gaussian_tab();
		break;
	case 2:
		multiple_gaussian_tab();
		break;
	default:
		std::cerr << "Non existent tab" << std::endl;
	}
}

void Kernels::multiple_gaussian_tab()
{
	std::cout << "Multiple gaussian tab selected" << std::endl;
}

void Kernels::gaussian_tab()
{
	std::cout << "Gaussian tab selected" << std::endl;
}

void Kernels::singular_tab()
{
	std::cout << "Singular tab selected" << std::endl;
}

void Kernels::initKernels()
{
	gaussianKernel = new Kernel();
	gaussianKernel->type = Kernel::GAUSSIAN;
	gaussianKernel->constant_shift = 0.000001f;
	gaussianKernel->sigma = 0.01f;

	singularKernel = new Kernel();
	singularKernel->type = Kernel::SINGULAR;
	singularKernel->constant_shift = 0.000001f;
	singularKernel->exponent = 3.f;
	singularKernel->sigma = 0.01f;

	multipleGaussianKernel = new Kernel();
	multipleGaussianKernel->type = Kernel::GAUSSIAN_MULTIPLE;
	multipleGaussianKernel->constant_shift = 0.000001f;
	multipleGaussianKernel->sigmas.resize(5);
	multipleGaussianKernel->sigmas[0] = 0.02f;
	multipleGaussianKernel->sigmas[1] = 0.05f;
	multipleGaussianKernel->sigmas[2] = 0.1f;
	multipleGaussianKernel->sigmas[3] = 0.3f;
	multipleGaussianKernel->sigmas[4] = 1.5f;

	currentKernel = new Kernel(*multipleGaussianKernel);
	m_firstUse = false;
}

void Kernels::kernelChanged()
{
	std::cout << "emission of kernel changed" << std::endl;
	emit kernel_changed();
}

void Kernels::add_gaussian()
{
}

void Kernels::remove_gaussian()
{
}
