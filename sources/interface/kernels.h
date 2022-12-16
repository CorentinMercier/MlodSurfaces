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

#ifndef KERNELS_H
#define KERNELS_H

#include "../apss/Fastapss.h"

#include <QWidget>
#include <iostream>


namespace Ui {
class Kernels;
}

class Kernels : public QWidget
{
	Q_OBJECT

public:
	explicit Kernels(QWidget *parent = nullptr);
	~Kernels();
	void initKernels();

	bool firstUse(){return m_firstUse;}

signals:
	void kernel_changed();

private slots:
	void kernelChanged();
	void add_gaussian();
	void remove_gaussian();
	void changed_tab();

private:
	void multiple_gaussian_tab();
	void gaussian_tab();
	void singular_tab();

	bool m_firstUse = true;

	Ui::Kernels *ui;

	Kernel *currentKernel;
	Kernel *gaussianKernel;
	Kernel *multipleGaussianKernel;
	Kernel *singularKernel;
};

#endif // KERNELS_H
