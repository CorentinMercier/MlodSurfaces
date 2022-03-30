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
