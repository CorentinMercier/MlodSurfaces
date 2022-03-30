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
