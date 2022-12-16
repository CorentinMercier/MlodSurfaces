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

#include "sliderWindow.h"

sliderWindow::sliderWindow(QWidget *parent) :
	QWidget(parent), ui(new Ui::sliderWindow)
{
	ui->setupUi(this);

	connect(ui->xPoint1, SIGNAL(sliderMoved(int)), this, SLOT(sliderMoved()));
	connect(ui->xPoint2, SIGNAL(sliderMoved(int)), this, SLOT(sliderMoved()));
	connect(ui->yPoint1, SIGNAL(sliderMoved(int)), this, SLOT(sliderMoved()));
	connect(ui->yPoint2, SIGNAL(sliderMoved(int)), this, SLOT(sliderMoved()));
	connect(ui->zPoint1, SIGNAL(sliderMoved(int)), this, SLOT(sliderMoved()));
	connect(ui->zPoint2, SIGNAL(sliderMoved(int)), this, SLOT(sliderMoved()));
	connect(ui->reset, SIGNAL(clicked()), this, SLOT(reset()));
}

sliderWindow::~sliderWindow()
{
	delete ui;
}


std::vector<float> sliderWindow::readSliders()
{
	std::vector<float> sliders(6);
	sliders[0] = (ui->xPoint1->value() - 50) / 50.f;
	sliders[1] = (ui->yPoint1->value() - 50) / 50.f;
	sliders[2] = (ui->zPoint1->value() - 50) / 50.f;
	sliders[3] = (ui->xPoint2->value() - 50) / 50.f;
	sliders[4] = (ui->yPoint2->value() - 50) / 50.f;
	sliders[5] = (ui->zPoint2->value() - 50) / 50.f;
	return sliders;
}

void sliderWindow::setSliders(std::vector<float> vals)
{
	ui->xPoint1->setValue(vals[0] * 50.f + 50);
	ui->yPoint1->setValue(vals[1] * 50.f + 50);
	ui->zPoint1->setValue(vals[2] * 50.f + 50);
	ui->xPoint2->setValue(vals[3] * 50.f + 50);
	ui->yPoint2->setValue(vals[4] * 50.f + 50);
	ui->zPoint2->setValue(vals[5] * 50.f + 50);
}

void sliderWindow::sliderMoved()
{
	emit slider_moved();
}

void sliderWindow::reset()
{
	ui->xPoint1->setValue(50);
	ui->xPoint2->setValue(50);
	ui->yPoint1->setValue(50);
	ui->yPoint2->setValue(50);
	ui->zPoint1->setValue(50);
	ui->zPoint2->setValue(50);
	emit slider_moved();
}
