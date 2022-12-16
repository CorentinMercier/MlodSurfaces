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

#include "mainwindow.h"
#include "ui_mainwindow.h"
#ifdef _WIN32
#include <io.h>
#define F_OK 0
#else
#include <sys/stat.h>
#include <unistd.h>
#endif

MainWindow::MainWindow(const CommandArgs &cmd, QWidget *parent) :
	QMainWindow(parent), ui(new Ui::MainWindow)
{
	//Setup window layout
	ui->setupUi(this);

	this->resize(1920*2/3,1280*2/3);
	this->setWindowTitle("FastAPSS");

	//Create openGL context
	QGLFormat qglFormat;
	//qglFormat.setVersion(1,2);
	qglFormat.setVersion(4,5);
	qglFormat.setDoubleBuffer(false);
	qglFormat.setSampleBuffers(false);
	//Create OpenGL Widget renderer
	glWidget=new WidgetOpenGL(qglFormat, cmd);

	//Add the OpenGL Widget into the layout
	ui->layout_scene->addWidget(glWidget);

	//loadParameters();

	//Connect slot and signals
	connect(ui->quit,SIGNAL(clicked()),this,SLOT(action_quit()));
	connect(ui->Reset, SIGNAL(clicked()),this,SLOT(action_reset_usual_octree()));
	connect(ui->Reset3, SIGNAL(clicked()),this,SLOT(action_reset_adaptive_octree()));
	connect(ui->ComputeStats, SIGNAL(clicked()), this, SLOT(compute_stats()));
	connect(ui->autoIteration, SIGNAL(clicked()), this, SLOT(auto_iteration()));
	connect(ui->dualContour, SIGNAL(clicked()), this, SLOT(compute_dual_contouring()));

	connect(ui->onUpsampled, SIGNAL(clicked()), this, SLOT(change_upsample_octree()));

	connect(ui->upsampling, SIGNAL(toggled(bool)), this, SLOT(upsampleChanged()));

	connect(ui->CPUMode, SIGNAL(toggled(bool)), this, SLOT(projectionMode()));
	connect(ui->GPUMode, SIGNAL(toggled(bool)), this, SLOT(projectionMode()));
	connect(ui->Fast, SIGNAL(toggled(bool)), this, SLOT(changeMethod()));
	connect(ui->Global, SIGNAL(toggled(bool)), this, SLOT(changeMethod()));
	connect(ui->Apss, SIGNAL(toggled(bool)), this, SLOT(changeMethod()));

	startTimer(0); //start timer
	secondTimer = new QTimer(this);
	connect(secondTimer, SIGNAL(timeout()), this, SLOT(secondTimerEvent()));
	secondTimer->start(100);//Start timer every 100ms

	fileMenu = menuBar()->addMenu(tr("&File"));
	openFile = new QAction(tr("&Open"), this);
	fileMenu->addAction(openFile);
	connect(openFile, SIGNAL(triggered()), this, SLOT(open()));
	saveSelection = new QAction(tr("&Save selection"), this);
	saveSelection->setEnabled(false);
	fileMenu->addAction(saveSelection);
	connect(saveSelection, SIGNAL(triggered()), this, SLOT(save_selection()));
	connect(glWidget, SIGNAL(thereIsASelection()), this, SLOT(enableSaveSelection()));
	connect(glWidget, SIGNAL(thereIsNoSelection()), this, SLOT(disableSaveSelection()));
	saveInputFile = new QAction(tr("&Save input file"), this);
	fileMenu->addAction(saveInputFile);
	connect(saveInputFile, SIGNAL(triggered()), this, SLOT(save_input_file()));
	loadPointsToProject = new QAction(tr("&Load points to project"), this);
	fileMenu->addAction(loadPointsToProject);
	connect(loadPointsToProject, SIGNAL(triggered()), this, SLOT(load_points_to_project()));
	savePointsToProject = new QAction(tr("&Save points to project"), this);
	fileMenu->addAction(savePointsToProject);
	connect(savePointsToProject, SIGNAL(triggered()), this, SLOT(save_points_to_project()));
	saveProjectedPoints = new QAction(tr("&Save projected points"), this);
	fileMenu->addAction(saveProjectedPoints);
	connect(saveProjectedPoints, SIGNAL(triggered()), this, SLOT(save_projected_points()));
	saveTimings = new QAction(tr("&Save timings"), this);
	fileMenu->addAction(saveTimings);
	connect(saveTimings, SIGNAL(triggered()), this, SLOT(save_timings()));
	quit = new QAction(tr("&Quit"), this);
	fileMenu->addAction(quit);
	connect(quit, SIGNAL(triggered()), this, SLOT(action_quit()));

	editMenu = menuBar()->addMenu(tr("&Edit"));
	parametersMenu = new QAction(tr("&Parameters"), this);
	editMenu->addAction(parametersMenu);
	connect(parametersMenu, SIGNAL(triggered()), this, SLOT(editParameters()));

	viewerMenu = menuBar()->addMenu(tr("&Viewer"));
	onlyWhenFinishedIteration = new QAction(tr("Update only after iterations"), this);
	onlyWhenFinishedIteration->setCheckable(true);
	onlyWhenFinishedIteration->setChecked(false);
	viewerMenu->addAction(onlyWhenFinishedIteration);
	connect(onlyWhenFinishedIteration, SIGNAL(triggered()), this, SLOT(display_condition()));
	normalDisplay = new QAction(tr("Display normals"), this);
	normalDisplay->setCheckable(true);
	normalDisplay->setChecked(false);
	viewerMenu->addAction(normalDisplay);
	connect(normalDisplay, SIGNAL(triggered()), this, SLOT(display_normal()));
	whiteBackground = new QAction(tr("White background"), this);
	whiteBackground->setCheckable(true);
	whiteBackground->setChecked(true);
	viewerMenu->addAction(whiteBackground);
	connect(whiteBackground, SIGNAL(triggered()), this, SLOT(background()));

	param = new parameters(this);
	param->setValues(glWidget->get_scene().getNbOfPtsToProject(), glWidget->get_scene().getmaxUpsampling(), glWidget->get_scene().getSizeOfGrid());
	param->setWindowFlags(Qt::Dialog);
	param->setFixedSize(param->geometry().width(), param->geometry().height());

	sw = new sliderWindow(this);
	sw->setWindowFlags(Qt::Dialog);
	sw->setFixedSize(sw->geometry().width(), sw->geometry().height());

	k = new Kernels(this);
	k->setWindowFlags(Qt::Dialog);

	connect(param, SIGNAL(record_values()), this, SLOT(recordValues()));
	connect(sw, SIGNAL(slider_moved()), this, SLOT(sliderMoved()));
	connect(k, SIGNAL(kernel_changed()), this, SLOT(updateKernel()));
	connect(glWidget, SIGNAL(update_UI()), this, SLOT(updateValues()));
	connect(glWidget, SIGNAL(connect_kernel()), this, SLOT(connectUIWithKernel()));

	if (glWidget->get_cmd().automatic)
	{
		m_autoIteration = true;
		if (glWidget->get_cmd().cpu || glWidget->get_cmd().knn || glWidget->get_cmd().ball)
			m_GPUMode = false;
	}
}

void MainWindow::connectUIWithKernel()
{
	connect(ui->sigma1, SIGNAL(valueChanged(double)), this, SLOT(updateKernel()));
	connect(ui->sigma2, SIGNAL(valueChanged(double)), this, SLOT(updateKernel()));
	connect(ui->sigma3, SIGNAL(valueChanged(double)), this, SLOT(updateKernel()));
	connect(ui->sigma4, SIGNAL(valueChanged(double)), this, SLOT(updateKernel()));
	connect(ui->sigma5, SIGNAL(valueChanged(double)), this, SLOT(updateKernel()));
	connect(ui->ExponentValue, SIGNAL(valueChanged(double)), this, SLOT(updateKernel()));
	connect(ui->SigmaValue, SIGNAL(valueChanged(double)), this, SLOT(updateKernel()));
	connect(ui->ConstantShiftValue, SIGNAL(valueChanged(double)), this, SLOT(updateKernel()));
	connect(ui->cutoffValue, SIGNAL(valueChanged(double)), this, SLOT(updateKernel()));
	connect(ui->computedSigmas1, SIGNAL(valueChanged(double)), this, SLOT(updateKernel()));
	connect(ui->computedSigmasFactor, SIGNAL(valueChanged(double)), this, SLOT(updateKernel()));
	connect(ui->computedSigmasNber, SIGNAL(valueChanged(int)), this, SLOT(updateKernel()));
	connect(ui->GaussianButton, SIGNAL(toggled(bool)), this, SLOT(updateKernel()));
	connect(ui->SingularButton, SIGNAL(toggled(bool)), this, SLOT(updateKernel()));
	connect(ui->MultipleGaussiansButton, SIGNAL(toggled(bool)), this, SLOT(updateKernel()));
	connect(ui->computedSigmas, SIGNAL(toggled(bool)), this, SLOT(updateKernel()));
	connect(ui->individualSigmas, SIGNAL(toggled(bool)), this, SLOT(updateKernel()));
	connect(ui->AutoUpdate, SIGNAL(toggled(bool)), this, SLOT(autoUpdateKernel()));

	connect(ui->maxDepth, SIGNAL(valueChanged(int)), this, SLOT(updateUpsamplingOctree()));
	connect(ui->cut, SIGNAL(toggled(bool)), this, SLOT(updateCut()));
	connect(ui->iterationValue, SIGNAL(valueChanged(int)), this, SLOT(updateMode()));
	connect(ui->minDepth, SIGNAL(valueChanged(int)), this, SLOT(updateMode()));
}

MainWindow::~MainWindow()
{
	if (glWidget)
	{
		secondTimer->stop();
		glWidget->close();
		delete glWidget;
		glWidget = nullptr;
	}
	delete ui;
	delete secondTimer;
	delete param;
	delete sw;
}

void MainWindow::recordAutomatic()
{
	cout << "Recording files ..." << endl;
	std::size_t endPath = glWidget->get_cmd().input_file.find_last_of('.');
	string nameOfFile, fileForProjectedPoints, fileForTimings, fileForDualContour;
	if (glWidget->get_cmd().output_file != "")
	{
		fileForProjectedPoints = glWidget->get_cmd().output_file;
		fileForTimings = glWidget->get_cmd().output_file;
		changeExtension(fileForTimings, ".timings");
		fileForDualContour = glWidget->get_cmd().output_file;
		changeExtension(fileForDualContour, "dc.ply");
	}
	else
	{
		if (glWidget->get_cmd().output_folder != "")
		{
			std::size_t beginPath = glWidget->get_cmd().input_file.find_last_of('/')+1;
			nameOfFile = glWidget->get_cmd().input_file.substr(beginPath, endPath - beginPath);
		}
		else
			nameOfFile = glWidget->get_cmd().input_file.substr(0, endPath);
		fileForProjectedPoints = glWidget->get_cmd().output_folder + nameOfFile + "Points.ply";
		fileForTimings = glWidget->get_cmd().output_folder + nameOfFile + ".timings";
		fileForDualContour = glWidget->get_cmd().output_folder + nameOfFile + ".ply";
	}
	glWidget->get_scene().saveProjectedPoints(fileForProjectedPoints);
	if (!glWidget->get_cmd().ball && !glWidget->get_cmd().knn)
		glWidget->get_scene().computeDualContouring(fileForDualContour);
	save_timings(fileForTimings);
	cout << "Computation finished, exiting" << endl;
	exit(EXIT_SUCCESS);
}

void MainWindow::timerEvent(QTimerEvent *event)
{
	event->accept();
	if (m_autoIteration)// && m_currentIteration < m_maxNbIterations)
	{
		if (m_currentIteration == m_maxNbIterations)
		{
			if (ui->upsampling->isChecked())
			{
				double mortonCodeTime;
				if (!glWidget->get_scene().upsample(mortonCodeTime))
				{
					m_autoIteration = false;
					if (glWidget->get_cmd().automatic)
						recordAutomatic();
					return;
				}
				else
					m_mortonCodeTimings.push_back(mortonCodeTime);
				m_currentIteration = 0;
			}
			else
			{
				m_autoIteration = false;
				if (glWidget->get_cmd().automatic)
					recordAutomatic();
				return;
			}
		}
		if (m_GPUMode)
			m_timings.push_back(glWidget->get_scene().projectPoints(m_displayAll || (m_currentIteration + 1) == m_maxNbIterations));
		else
			m_timings.push_back(glWidget->get_scene().projectPointsOnCPU(m_displayAll || (m_currentIteration + 1) == m_maxNbIterations));
		m_currentIteration++;
	}
}

void MainWindow::secondTimerEvent()
{
	double time = glWidget->drawTime();
	std::ostringstream out, out2;
	out.precision(2);
	out << std::fixed << time;
	double fps = 1000.0 / time;
	out2.precision(2);
	out2 << std::fixed << fps;
	string textToPrint = out.str() + "ms     " + out2.str() + "fps";
	ui->informations->setText(QString::fromStdString(textToPrint));
	if (!m_computing)
	{
		std::ostringstream it;
		it << "Iteration " << m_currentIteration << "/" << m_maxNbIterations;
		string sndTextToPrint = it.str();
		ui->informations2->setText(QString::fromStdString(sndTextToPrint));
	}
}

void MainWindow::updateMode()
{
	cout << "updating mode..." << endl;
	if (m_maxNbIterations != unsigned(ui->iterationValue->value()))
	{
		m_maxNbIterations = unsigned(ui->iterationValue->value());
		glWidget->get_scene().initUpsampling();
		reset();
		m_currentIteration = 0;
	}
	if (m_minDepth != unsigned(ui->minDepth->value()))
	{
		m_minDepth = unsigned(ui->minDepth->value());
		glWidget->get_scene().updateMinDepthValue(m_minDepth);
		glWidget->get_scene().initUpsampling();
		reset();
		m_currentIteration = 0;
	}
	if (ui->autoIteration->isChecked() && !m_autoIteration)
	{
		reset();
		m_autoIteration = true;
		m_currentIteration = 0;
	}
}

void MainWindow::updateUpsamplingOctree()
{
	glWidget->get_scene().setCut(ui->cut->isChecked());
	glWidget->get_scene().setMaxDepth(unsigned(ui->maxDepth->value()));
	m_mortonCodeTimings.clear();
	m_mortonCodeTimings.push_back(glWidget->get_scene().resetOctree());
	m_currentIteration = 0;
	if (ui->autoIteration->isChecked() && !m_autoIteration)
		m_autoIteration = true;
}

void MainWindow::updateCut()
{
	glWidget->get_scene().setCut(ui->cut->isChecked());
}

void MainWindow::autoUpdateKernel()
{
	if (ui->AutoUpdate->isChecked())
	{
		reset();
		if (ui->autoIteration->isChecked() && !m_autoIteration)
			m_autoIteration = true;
	}
}

void MainWindow::updateKernel()
{
	cout << "updating kernel..." << endl;
	glWidget->get_scene().initUpsampling();
	m_currentIteration = 0;
	Kernel newKernel;
	newKernel.type = typeChecked();
	newKernel.sigma = float(ui->SigmaValue->value());
	if (ui->individualSigmas->isChecked())
	{
		newKernel.sigmas.resize(5);
		newKernel.sigmas[0] = float(ui->sigma1->value());
		newKernel.sigmas[1] = float(ui->sigma2->value());
		newKernel.sigmas[2] = float(ui->sigma3->value());
		newKernel.sigmas[3] = float(ui->sigma4->value());
		newKernel.sigmas[4] = float(ui->sigma5->value());
	}
	else
	{
		newKernel.sigmas.resize(ui->computedSigmasNber->value());
		newKernel.sigmas[0] = ui->computedSigmas1->value();
		for (unsigned int i=1; i<newKernel.sigmas.size(); i++)
			newKernel.sigmas[i] = newKernel.sigmas[i - 1] * ui->computedSigmasFactor->value();
	}
	newKernel.beginningSigma = ui->computedSigmas1->value();
	newKernel.factor = ui->computedSigmasFactor->value();
	newKernel.exponent = float(ui->ExponentValue->value());
	newKernel.constant_shift = float(ui->ConstantShiftValue->value());
	glWidget->get_scene().setCutOff(float(ui->cutoffValue->value()));
	glWidget->get_scene().setNewKernel(newKernel);

	if (ui->autoIteration->isChecked() && !m_autoIteration)
		m_autoIteration = true;
	if (ui->AutoUpdate->isChecked())
		reset();
	else
		m_autoIteration = false;
}

void MainWindow::reset()
{
	m_mortonCodeTimings.clear();
	m_mortonCodeTimings.push_back(glWidget->get_scene().resetPoints(m_lastReset == lastReset::OCTREE));
}

void MainWindow::auto_iteration()
{
	m_autoIteration = ui->autoIteration->isChecked();
	if (m_autoIteration && m_currentIteration == m_maxNbIterations)
	{
		reset();
		m_currentIteration = 0;
	}
}

void MainWindow::action_reset_usual_octree()
{
	if (!m_computing)
	{
		m_timings.clear();
		m_mortonCodeTimings.clear();
		m_lastReset = lastReset::RANDOM;
		m_mortonCodeTimings.push_back(glWidget->get_scene().resetPoints(false));
		m_currentIteration = 0;
		if (ui->autoIteration->isChecked() && !m_autoIteration)
			m_autoIteration = true;
	}
}

void MainWindow::action_reset_adaptive_octree()
{
	if (!m_computing)
	{
		m_timings.clear();
		m_mortonCodeTimings.clear();
		m_lastReset = lastReset::OCTREE;
		m_mortonCodeTimings.push_back(glWidget->get_scene().resetPoints(true));
		m_currentIteration = 0;
		if (ui->autoIteration->isChecked() && !m_autoIteration)
			m_autoIteration = true;
	}
}

void MainWindow::change_upsample_octree()
{
	glWidget->get_scene().setOnUpsample(ui->onUpsampled->isChecked());
}

void MainWindow::projectionMode()
{
	m_GPUMode = ui->GPUMode->isChecked();
}

void MainWindow::changeMethod()
{
	m_knnMode = false;
	if (ui->Apss->isChecked())
	{
		ui->CPUMode->setChecked(true);
		ui->GPUMode->setChecked(false);
		ui->GPUMode->setCheckable(false);
		glWidget->get_scene().knnMode(true);
		m_knnMode = true;
	}
	else
	{
		glWidget->get_scene().knnMode(false);
		ui->GPUMode->setCheckable(true);
		if (ui->Global->isChecked())//GlobalAPSS
			updateMode();
		if (m_method != ui->Fast->isChecked())
		{
			m_method = ui->Fast->isChecked();
			glWidget->get_scene().changeMethod(m_method);
		}
	}
}

void MainWindow::changeRadius(int val)
{
	glWidget->changeRadiusFactor(val / 1000.f);
}

void MainWindow::upsampleChanged()
{
	if (ui->upsampling->isChecked() && ui->autoIteration->isChecked())
		m_autoIteration = true;
}

void MainWindow::display_condition()
{
	m_displayAll = !m_displayAll;
}

void MainWindow::display_normal()
{
	glWidget->get_scene().normalMode = normalDisplay->isChecked();
}

void MainWindow::background()
{
	glWidget->changeBackground(whiteBackground->isChecked());
}

void MainWindow::compare_point_set()
{
	glWidget->stopTimer();
	string filename = QFileDialog::getOpenFileName(this, "Open a file containing a point set", ".", "pointset files (*.pn *.ply)").toStdString();
	glWidget->restartTimer();
	if (filename == "")
	{
		cerr << "Error, no file selected" << endl;
		return;
	}
	glWidget->get_scene().comparePointSet(filename);
}

void MainWindow::compute_stats()
{
	m_computing = true;
	glWidget->stopTimer();
	string filename = QFileDialog::getSaveFileName(this , "Choose a file in which to record stats", ".", "txt files (*.txt)", nullptr, QFileDialog::DontUseNativeDialog).toStdString();
	glWidget->restartTimer();

	glWidget->updateGL();
	update();
	QCoreApplication::processEvents();
	QCoreApplication::sendPostedEvents();

	fstream file;
	file.open(filename, ios_base::out);
	if (!file.is_open())
	{
		std::cerr << "Impossible to write file " << filename << std::endl;
		return;
	}
	unsigned int nbOfProjections = m_maxNbIterations;
	for (unsigned int i=1000; i<glWidget->get_scene().sizeOfPointSet(); i+=1000)
	{
		std::ostringstream it;
		it << "Stats computation " << i/(float)glWidget->get_scene().sizeOfPointSet() * 100 << "%";
		string textToPrint = it.str();
		ui->informations2->setText(QString::fromStdString(textToPrint));
		//Resample pointSet
		glWidget->get_scene().resamplePointSet(i);
		//Reset projected points
		reset();
		//Project 10 times
		double time = 0.0;
		for (unsigned int j=0; j<nbOfProjections; j++)
		{
			time += glWidget->get_scene().projectPoints();
			glWidget->updateGL();
			update();
			QCoreApplication::processEvents();
		}
		//Compute average
		time /= double(nbOfProjections);
		//Write value in file
		file << i << " " << time << std::endl;
	}
	file.close();
	glWidget->get_scene().resamplePointSet(glWidget->get_scene().sizeOfPointSet());
	//Reset projected points
	reset();
	m_computing = false;
}

void MainWindow::compute_dual_contouring()
{
	glWidget->stopTimer();
	string filename = QFileDialog::getSaveFileName(this, "Choose a file to save the dual contour", "../data/", "ply files (*.ply)", nullptr, QFileDialog::DontUseNativeDialog).toStdString();
	glWidget->restartTimer();
	if (filename == "")
	{
		cerr << "Error, no file selected" << endl;
		return;
	}
	glWidget->get_scene().computeDualContouring(filename, false);
}

void MainWindow::action_quit()
{
	saveParameters();
	glWidget->get_scene().stop();
	cudaDeviceReset();
	close();
}

Kernel::KernelType MainWindow::typeChecked()
{
	if (ui->SingularButton->isChecked())
		return Kernel::SINGULAR;
	if (ui->GaussianButton->isChecked())
		return Kernel::GAUSSIAN;
	if(ui->MultipleGaussiansButton->isChecked())
		return Kernel::GAUSSIAN_MULTIPLE;
	return Kernel::SINGULAR;
}

void MainWindow::open()
{
	bool autoIteration = m_autoIteration;
	m_autoIteration = false;
	saveParameters();
	glWidget->openNewScene();
	m_autoIteration = autoIteration;
}

void MainWindow::save_selection()
{
	glWidget->stopTimer();
	string filename = QFileDialog::getSaveFileName(this , "Choose a file in which to record selection", "../data/", "ply files (*.ply)", nullptr, QFileDialog::DontUseNativeDialog).toStdString();
	glWidget->restartTimer();
	if (filename == "")
	{
		cerr << "Error, no file selected" << endl;
		return;
	}
	glWidget->get_scene().saveSelection(filename);
}

void MainWindow::save_input_file()
{
	glWidget->stopTimer();
	string filename = QFileDialog::getSaveFileName(this , "Choose a file in which to record input file", "../data/", "ply files (*.ply)", nullptr, QFileDialog::DontUseNativeDialog).toStdString();
	glWidget->restartTimer();
	if (filename == "")
	{
		cerr << "Error, no file selected" << endl;
		return;
	}
	glWidget->get_scene().saveInputFile(filename);
}

void MainWindow::load_points_to_project()
{
	glWidget->get_scene().loadPointsToProject();
	glWidget->update();
}

void MainWindow::save_points_to_project()
{
	glWidget->get_scene().savePointsToProject();
}

void MainWindow::save_projected_points()
{
	glWidget->stopTimer();
	string filename = QFileDialog::getSaveFileName(this , "Choose a file in which to record projected points", "../data/", "ply files (*.ply)", nullptr, QFileDialog::DontUseNativeDialog).toStdString();
	glWidget->restartTimer();
	if (filename == "")
	{
		cerr << "Error, no file selected" << endl;
		return;
	}
	glWidget->get_scene().saveProjectedPoints(filename);
}

void MainWindow::save_triangles()
{
	glWidget->get_scene().createTrianglesForMPU();
}

void MainWindow::save_timings(string f)
{
	if (m_timings.size() == 0)
	{
		cerr << "No timings to record" << endl;
		return;
	}
	glWidget->stopTimer();
	string filename;
	if (f == "")
	{
		filename = glWidget->getFilename();
		changeExtension(filename, ".timings");
	}
	else
		filename = f;
	unsigned nb = 0;
	while(access( filename.c_str(), F_OK ) != -1) /// TODO: make cross platform
	{
		nb++;
		changeExtension(filename, ".timings" + to_string(nb));
	}
	glWidget->restartTimer();
	fstream file;
	file.open(filename, ios_base::out);
	if (!file.is_open())
	{
		cerr << "Impossible to create file : " << filename << endl;
		cerr << "No timings were recorded" << endl;
		return;
	}
	file << "File: " << glWidget->getFilename() << endl;
	file << "Nb of input points: " << glWidget->get_scene().sizeOfPointSet() << endl;
	file << "Nb of projected points: " << glWidget->get_scene().getNbOfPtsToProject() << endl;
	if (m_GPUMode)
		file << "GPU - Mode: ";
	else
		file << "CPU - Mode: ";
	if (m_knnMode)
		file << "knn APSS" << endl;
	else if (m_method)//True == fast
		file << "FastAPSS" << endl;
	else
		file << "Global APSS" << endl;
	double averagePerPoint = 0;
	double standardDeviation = 0;
	for (unsigned i=0; i<m_timings.size(); i++)
	{
		averagePerPoint += m_timings[i];
		standardDeviation += m_timings[i] * m_timings[i];
	}
	averagePerPoint /= double(m_timings.size() * glWidget->get_scene().getNbOfPtsToProject());
	standardDeviation = sqrt(standardDeviation / double(m_timings.size() * glWidget->get_scene().getNbOfPtsToProject() * glWidget->get_scene().getNbOfPtsToProject()) - averagePerPoint * averagePerPoint);
	file << "Timings for projection per point: " << averagePerPoint << "+-" << standardDeviation << "ms" <<  endl;
	double averageMortonCode = 0;
	double standardDeviationMorton = 0;
	for (unsigned i=0; i<m_mortonCodeTimings.size(); i++)
	{
		averageMortonCode += m_mortonCodeTimings[i];
		standardDeviationMorton += m_mortonCodeTimings[i] * m_mortonCodeTimings[i];
	}
	averageMortonCode /= double(m_mortonCodeTimings.size());
	standardDeviationMorton = sqrt(standardDeviationMorton / double(m_mortonCodeTimings.size()) - averageMortonCode * averageMortonCode);
	file << "Timings for Morton Code: " << averageMortonCode << "+-" << standardDeviationMorton << "ms" <<  endl;

	file << "Octrees max depth: " << glWidget->get_scene().getMaxDepth() << endl;
	file << "Adaptive octree: " << glWidget->get_scene().getAdaptiveOctreeTime() << endl;
	file << "Usual octree: " << glWidget->get_scene().getUsualOctreeTime() << endl;
	size_t nodesOnProjectedPoints =  glWidget->get_scene().count_nodes_on_projected_points();//This computes the time to uniformize
	file << "Uniformize octree: " << glWidget->get_scene().getUniformizeTime() << endl;
	file << "Crossed nodes: " << glWidget->get_scene().getNbOfCrossedNodes() << endl;
	file << "Nodes created on projected points: " << nodesOnProjectedPoints << endl;


	file.close();
	cout << "File " << filename << " saved !" << endl;
}

void MainWindow::enableSaveSelection()
{
	saveSelection->setEnabled(true);
}

void MainWindow::disableSaveSelection()
{
	saveSelection->setEnabled(false);
}

void MainWindow::editInitBox()
{
	sw->show();
}

void MainWindow::editParameters()
{
	bool autoIteration = m_autoIteration;
	m_autoIteration = false;
	param->show();
	m_autoIteration = autoIteration;
}

void MainWindow::editKernels()
{
	if (k->firstUse())
		k->initKernels();
	k->show();
}

void MainWindow::planeMode()
{
	glWidget->changePlaneMode();
}

void MainWindow::computeNoise()
{
	glWidget->get_scene().computeNoise();
}

void MainWindow::recordValues()
{
	unsigned nbOfPts, maxUpsampling, sizeOfGrid;
	param->getValues(nbOfPts, maxUpsampling, sizeOfGrid);
	glWidget->get_scene().setSizeOfGrid(sizeOfGrid);
	if (glWidget->get_scene().getNbOfPtsToProject() == nbOfPts && glWidget->get_scene().getmaxUpsampling() == maxUpsampling)
		return;
	glWidget->get_scene().setNbOfPtsToProject(nbOfPts);
	glWidget->get_scene().setMaxUpsampling(maxUpsampling);
	glWidget->get_scene().setOffset(sw->readSliders());
	glWidget->get_scene().eraseFromGPU();
	glWidget->get_scene().initBuffers();
}

void MainWindow::sliderMoved()
{
}

void MainWindow::saveParameters()
{
	string filename = glWidget->getFilename();
	changeExtension(filename, ".params");
	fstream file;
	file.open(filename, ios_base::out | ios_base::binary);
	if (!file.is_open())
	{
		cerr << "Impossible to create file : " << filename << endl;
		cerr << "No parameters were recorded" << endl;
		return;
	}
	glm::fquat q = glWidget->cam.getQuat();
	file.write(reinterpret_cast<const char *>(&q), sizeof(glm::fquat));
	glm::vec3 t=glWidget->cam.getTranslation();
	file.write(reinterpret_cast<const char *>(&t), sizeof(glm::vec3));
	float d=glWidget->cam.getDist();
	file.write(reinterpret_cast<const char *>(&d), sizeof(float));

	float val = float(ui->SigmaValue->value());
	file.write(reinterpret_cast<const char *>(&val), sizeof(float));
	val = float(ui->sigma1->value());
	file.write(reinterpret_cast<const char *>(&val), sizeof(float));
	val = float(ui->sigma2->value());
	file.write(reinterpret_cast<const char *>(&val), sizeof(float));
	val = float(ui->sigma3->value());
	file.write(reinterpret_cast<const char *>(&val), sizeof(float));
	val = float(ui->sigma4->value());
	file.write(reinterpret_cast<const char *>(&val), sizeof(float));
	val = float(ui->sigma5->value());
	file.write(reinterpret_cast<const char *>(&val), sizeof(float));
	val = float(ui->computedSigmas1->value());
	file.write(reinterpret_cast<const char *>(&val), sizeof(float));
	val = float(ui->computedSigmasFactor->value());
	file.write(reinterpret_cast<const char *>(&val), sizeof(float));
	val = float(ui->ConstantShiftValue->value());
	file.write(reinterpret_cast<const char *>(&val), sizeof(float));
	val = float(ui->ExponentValue->value());
	file.write(reinterpret_cast<const char *>(&val), sizeof(float));
	val = float(ui->cutoffValue->value());
	file.write(reinterpret_cast<const char *>(&val), sizeof(float));

	unsigned v = unsigned(ui->computedSigmasNber->value());
	file.write(reinterpret_cast<const char *>(&v), sizeof(unsigned));
	v = unsigned(ui->iterationValue->value());
	file.write(reinterpret_cast<const char *>(&v), sizeof(unsigned));
	v = unsigned(ui->maxDepth->value());
	file.write(reinterpret_cast<const char *>(&v), sizeof(unsigned));
	v = unsigned(ui->minDepth->value());
	file.write(reinterpret_cast<const char *>(&v), sizeof(unsigned));
	v = unsigned(glWidget->get_scene().getMode());
	file.write(reinterpret_cast<const char *>(&v), sizeof(unsigned));

	Kernel::KernelType type = typeChecked();
	file.write(reinterpret_cast<const char *>(&type), sizeof(Kernel::KernelType));
	bool computed = ui->computedSigmas->isChecked();
	file.write(reinterpret_cast<const char *>(&computed), sizeof(bool));
	bool cut = ui->cut->isChecked();
	file.write(reinterpret_cast<const char *>(&cut), sizeof(bool));

	unsigned nbOfPtsToProject = glWidget->get_scene().getNbOfPtsToProject();
	file.write(reinterpret_cast<const char *>(&nbOfPtsToProject), sizeof(unsigned));
	unsigned nbOfUpsampling = glWidget->get_scene().getmaxUpsampling();
	file.write(reinterpret_cast<const char *>(&nbOfUpsampling), sizeof(unsigned));
	std::vector<float> sliders = sw->readSliders();
	file.write(reinterpret_cast<const char *>(sliders.data()), 6 * sizeof(float));

	file.close();
	cout << "Parameters saved in : " << filename << endl;
}

void MainWindow::updateValues()
{
	ui->maxDepth->setValue(int(glWidget->get_scene().getMaxDepth()));
	ui->cut->setChecked(glWidget->get_scene().getCut());
	m_maxNbIterations = glWidget->get_scene().getMaxNbIterations();
	ui->iterationValue->setValue(int(m_maxNbIterations));
	m_currentIteration = 0;
	m_minDepth = glWidget->get_scene().getMinDepth();
	ui->minDepth->setValue(int(m_minDepth));

	Kernel newKernel = glWidget->get_scene().getKernel();
	ui->SigmaValue->setValue(double(newKernel.sigma));
	ui->sigma1->setValue(double(newKernel.fixedSigmas[0]));
	ui->sigma2->setValue(double(newKernel.fixedSigmas[1]));
	ui->sigma3->setValue(double(newKernel.fixedSigmas[2]));
	ui->sigma4->setValue(double(newKernel.fixedSigmas[3]));
	ui->sigma5->setValue(double(newKernel.fixedSigmas[4]));
	ui->computedSigmas1->setValue(double(newKernel.beginningSigma));
	ui->computedSigmasFactor->setValue(double(newKernel.factor));
	ui->ConstantShiftValue->setValue(double(newKernel.constant_shift));
	ui->ExponentValue->setValue(double(newKernel.exponent));
	ui->cutoffValue->setValue(double(glWidget->get_scene().getCutOff()));
	ui->computedSigmasNber->setValue(int(newKernel.nbOfComputedSigmas));

	Kernel::KernelType type = newKernel.type;
	ui->GaussianButton->setChecked(type == Kernel::KernelType::GAUSSIAN);
	ui->MultipleGaussiansButton->setChecked(type == Kernel::KernelType::GAUSSIAN_MULTIPLE);
	ui->SingularButton->setChecked(type == Kernel::KernelType::SINGULAR);

	ui->computedSigmas->setChecked(glWidget->get_scene().getComputedSigma());
	ui->individualSigmas->setChecked(!glWidget->get_scene().getComputedSigma());

	m_method = glWidget->get_scene().getMode();
	ui->Fast->setChecked(m_method);
	ui->Global->setChecked(!m_method);

	sw->setSliders(glWidget->get_scene().getOffset());
	param->setValues(glWidget->get_scene().getNbOfPtsToProject(), glWidget->get_scene().getmaxUpsampling(), glWidget->get_scene().getSizeOfGrid());
}

void MainWindow::closeEvent(QCloseEvent* event)
{
	secondTimer->stop();
	glWidget->close();
	delete glWidget;
	glWidget = nullptr;
}
