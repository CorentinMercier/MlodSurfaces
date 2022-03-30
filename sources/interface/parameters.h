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

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <QWidget>
#include <iostream>

namespace Ui {
class parameters;
}

class parameters : public QWidget
{
	Q_OBJECT

public:
	explicit parameters(QWidget *parent = 0);
	~parameters();

	void setValues(unsigned nbOfPtsToProject, unsigned maxUpsampling, unsigned sizeOfGrid);
	void getValues(unsigned &nbOfPts, unsigned &maxUpsampling, unsigned &sizeOfGrid);

signals:
	void record_values();

private slots:
	void recordValues();

private:
	Ui::parameters *ui;
};

#endif // PARAMETERS_H
