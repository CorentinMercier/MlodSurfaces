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
