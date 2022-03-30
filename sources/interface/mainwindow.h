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

#pragma once

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "kernels.h"
#include "widgetopengl.h"
#include <QTimer>
#include <QMainWindow>
#include <QThread>
#include <sstream>
#include <fstream>
#include "sliderWindow.h"
#include "parameters.h"

inline void changeExtension(string & path, const string & newExtension)
{
	std::size_t endPath = path.find_last_of('.');
	path = path.substr(0, endPath) + newExtension;
}

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	explicit MainWindow(const CommandArgs &cmd, QWidget *parent = nullptr);
	~MainWindow();

private slots:

	/** Quit the application */
	void action_quit();

	void open();
	void save_selection();
	void save_input_file();
	void save_points_to_project();
	void load_points_to_project();
	void save_projected_points();
	void save_triangles();
	void save_timings(string f = "");
	void editInitBox();
	void editParameters();
	void editKernels();
	void planeMode();
	void computeNoise();
	void enableSaveSelection();
	void disableSaveSelection();

	void auto_iteration();
	void action_reset_usual_octree();
	void action_reset_adaptive_octree();
	void change_upsample_octree();

	void projectionMode();
	void changeMethod();
	void changeRadius(int val);

	void upsampleChanged();
	void display_condition();
	void display_normal();
	void background();

	void compare_point_set();

	void autoUpdateKernel();
	void updateKernel();
	void updateMode();
	void updateUpsamplingOctree();
	void updateCut();
	void connectUIWithKernel();

	void compute_stats();
	void compute_dual_contouring();

	/** Function called in a timer loop */
	void timerEvent(QTimerEvent *event);
	/**Loop for checking user modifications - invoqued less often */
	void secondTimerEvent();

	void recordValues();
	void sliderMoved();

	void saveParameters();

	void updateValues();

private:
	void recordAutomatic();

	/** Layout for the Window */
	Ui::MainWindow *ui;
	/** The OpenGL Widget */
	WidgetOpenGL *glWidget;

	//Menu bar
	QMenu *fileMenu;
	QAction *openFile;
	QAction *saveSelection;
	QAction *saveInputFile;
	QAction *loadPointsToProject;
	QAction *savePointsToProject;
	QAction *saveProjectedPoints;
	QAction *saveTriangles;
	QAction *saveTimings;
	QAction *quit;

	QMenu *editMenu;
	QAction *editBox;
	QAction *parametersMenu;
	QAction *kernelsMenu;
	QAction *planesMenu;
	QAction *noiseMenu;

	QMenu *display;
	QAction *originalPointSet;
	QAction *projectedPoints;
	QAction *centers;

	QMenu *viewerMenu;
	QAction *onlyWhenFinishedIteration;
	QAction *normalDisplay;
	QAction *whiteBackground;

	QMenu *compareMenu;
	QAction *comparePointSet;

	Kernel::KernelType typeChecked();

	bool m_GPUMode = true;
	bool m_method = true; //True = FastAPSS
	bool m_knnMode = false;

	//Values from UI
	unsigned m_framesForSampling = 4;
	unsigned m_nbPtsForSampling = 20000;

	unsigned int m_maxNbIterations = 10;
	unsigned int m_currentIteration = 0;
	unsigned int m_minDepth = 1;

	bool m_autoIteration = false;
	bool m_displayAll = true;

	QTimer * secondTimer;

	enum lastReset{RANDOM, FITTED, OCTREE};
	lastReset m_lastReset = lastReset::OCTREE;
	void reset();
	bool m_computing = false;

	//Timings
	vector<double> m_timings;
	vector<double> m_mortonCodeTimings;

	//Windows
	sliderWindow *sw;
	parameters *param;
	Kernels *k;
};

#endif // MAINWINDOW_H
