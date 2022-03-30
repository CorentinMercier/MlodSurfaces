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
