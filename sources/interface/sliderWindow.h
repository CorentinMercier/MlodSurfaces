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

#ifndef SLIDERWINDOW_H
#define SLIDERWINDOW_H

#include <QWidget>
#include <QSlider>
#include "ui_sliderWindow.h"

namespace Ui {
class sliderWindow;
}

class sliderWindow : public QWidget
{
	Q_OBJECT

public:
	explicit sliderWindow(QWidget *parent = nullptr);
	~sliderWindow();

	std::vector<float> readSliders();
	void setSliders(std::vector<float>);

signals:
	void slider_moved();

private slots:
	void sliderMoved();
	void reset();

private:
	Ui::sliderWindow *ui;
};

#endif // SLIDERWINDOW_H
