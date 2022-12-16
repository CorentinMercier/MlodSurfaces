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

#include "parameters.h"
#include "ui_parameters.h"

parameters::parameters(QWidget *parent) :
	QWidget(parent),
	ui(new Ui::parameters)
{
	ui->setupUi(this);

	connect(ui->finishButtons, SIGNAL(rejected()), this, SLOT(hide()));
	connect(ui->finishButtons, SIGNAL(accepted()), this, SLOT(recordValues()));
}

parameters::~parameters()
{
	delete ui;
}

void parameters::setValues(unsigned nbOfPtsToProject, unsigned maxUpsampling, unsigned sizeOfGrid)
{
	ui->nbOfPtsToProject->setValue(nbOfPtsToProject);
	ui->maxUpsampling->setValue(maxUpsampling);
	ui->sizeOfMortonCodeGrid->setValue(sizeOfGrid);
}

void parameters::recordValues()
{
	emit record_values();
	hide();
}

void parameters::getValues(unsigned & nbOfPts, unsigned & maxUpsampling, unsigned & sizeOfGrid)
{
	nbOfPts = unsigned(ui->nbOfPtsToProject->value());
	maxUpsampling = unsigned(ui->maxUpsampling->value());
	sizeOfGrid = unsigned(ui->sizeOfMortonCodeGrid->value());
}
