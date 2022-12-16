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

#include "widgetopengl.h"

WidgetOpenGL::WidgetOpenGL(const QGLFormat& format, const CommandArgs& cmd, QGLWidget *parent) :
	QGLWidget(format, parent), m_cmd(cmd)
{
	if (!m_cmd.automatic)
		m_filename = QFileDialog::getOpenFileName(this, "Open a file containing a point set", "../data/", "pointset files (*.pn *.ply)").toStdString();
	else
		m_filename = m_cmd.input_file;
	QWidget::setFocusPolicy(Qt::WheelFocus);
	m_mainTimer = startTimer(0);
}

WidgetOpenGL::~WidgetOpenGL()
{

}

void WidgetOpenGL::initializeGL()
{
	//Init OpenGL
	setup_opengl();

	//Init Camera
	cam.setupCamera();

	//Init Scene 3D
	m_scene.set_widget(this);
	//emit(scene_loaded());
	if (!m_scene.load_scene(m_filename))
	{
		this->window()->close();
		exit(0);
	}
	if (m_cmd.automatic && m_cmd.customKernel)// for kernel parameters in command line
	{
		m_cmd.kernel.updateSigmas(true);
		m_scene.setNewKernel(m_cmd.kernel);
		m_scene.setCutOff(m_cmd.cutoff);
	}
	//Activate depth buffer
	glEnable(GL_DEPTH_TEST); PRINT_OPENGL_ERROR();

	if (m_cmd.automatic)
	{
		m_scene.changeMethod(m_cmd.fast);
		if (m_scene.getMaxDepth() != m_cmd.depth)
			m_scene.setMaxDepth(m_cmd.depth);
		m_scene.knnMode(m_cmd.knn, m_cmd.nb_knn);
		m_scene.ballMode(m_cmd.ball, m_cmd.ball_radius);
		m_scene.resetOctree(!m_cmd.inputOctree);
	}

	if (m_cmd.dc)
	{
		std::size_t endPath = m_filename.find_last_of('.');
		string nameOfFile, fileForDualContour;
		if (m_cmd.output_file != "")
			fileForDualContour = m_cmd.output_file;
		else
		{
			if (m_cmd.output_folder != "")
			{
				std::size_t beginPath = m_filename.find_last_of('/')+1;
				nameOfFile = m_filename.substr(beginPath, endPath - beginPath);
			}
			else
				nameOfFile = m_filename.substr(0, endPath);
			fileForDualContour = m_cmd.output_folder + nameOfFile + ".ply";
		}
		m_scene.computeDualContouring(fileForDualContour);
		exit(EXIT_SUCCESS);
	}
}

bool WidgetOpenGL::openNewScene()
{
	killTimer(m_mainTimer);
	string filename = QFileDialog::getOpenFileName(this, "Open a file containing a point set", "../data/", "pointset files (*.pn *.ply)").toStdString();
	m_mainTimer = startTimer(0);
	if (filename == "")
	{
		cerr << "Error, no file selected" << endl;
		return false;
	}
	emit thereIsNoSelection();
	m_filename = filename;
	m_scene.reset();
	if (!m_scene.load_scene(m_filename))
	{
		this->window()->close();
		exit(0);
		return false;
	}
	return true;
}

void WidgetOpenGL::paintGL()
{
	//Compute current cameras
	cam.setupCamera();

	//clear screen
	glViewport (0, 0, cam.getScreenSizeX(), cam.getScreenSizeY()); PRINT_OPENGL_ERROR();
	if (m_whiteBackground)
		glClearColor (1.0f, 1.0f, 1.0f, 1.0f);
	else
		glClearColor (0.2f, 0.2f, 0.3f, 1.0f);
	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);        PRINT_OPENGL_ERROR();

	m_scene.draw_scene();
	m_scene.drawTime();
}


void WidgetOpenGL::resizeGL(int const width, int const height)
{
	cam.setScreenSize(width, height);
	glViewport(0,0, width, height); PRINT_OPENGL_ERROR();
}

void WidgetOpenGL::setup_opengl()
{
	setup_glad();
	print_current_opengl_context();
}

void WidgetOpenGL::setup_glad()
{
	if(!gladLoadGL())
	{
		std::cerr<<"Error initializing GLAD\n";
		exit(EXIT_FAILURE);
	}
	std::cout << "GLAD initialized\n";
}

void WidgetOpenGL::print_current_opengl_context() const
{
	std::cout << "OpenGl informations: VENDOR:       " << glGetString(GL_VENDOR)<<std::endl;
	std::cout << "                     RENDERER:     " << glGetString(GL_RENDERER)<<std::endl;
	std::cout << "                     VERSION:      " << glGetString(GL_VERSION)<<std::endl;
	std::cout << "                     GLSL VERSION: " << glGetString(GL_SHADING_LANGUAGE_VERSION)<<std::endl;
	std::cout<<"Current OpenGL context: "<< context()->format().majorVersion() << "." << context()->format().minorVersion()<<std::endl;
}

scene& WidgetOpenGL::get_scene()
{
	return m_scene;
}

void WidgetOpenGL::keyPressEvent(QKeyEvent *event)
{
	int current=event->key();
	Qt::KeyboardModifiers mod=event->modifiers();

	// We can quit the scene with 'Q'
	if( (mod&Qt::ShiftModifier)!=0 && (current==Qt::Key_Q) )
	{
		std::cout<<"\n[EXIT OK]\n\n"<<std::endl;
		this->window()->close();
	}
	if( (mod&Qt::ShiftModifier)!=0 && (current==Qt::Key_H) )
	{
		std::cout<<"\nDisplay traversed nodes\n"<<std::endl;
		m_scene.output_traversed_nodes_statistics();
	}
	if ((mod&Qt::ShiftModifier)!=0 && current==Qt::Key_S)
	{
		fstream file;
		file.open("Camera.txt" , ios_base::out);
		glm::fquat q = cam.getQuat();
		file << q.x << " " << q.y << " " << q.z << " " << q.w << endl;
		glm::vec3 t=cam.getTranslation();
		file << t[0] << " " << t[1] << " " << t[2] << endl;
		file << cam.getDist() << endl;
		file.close();
		cout << "Position de la caméra enregistrée" << endl;
	}
	if ((mod&Qt::ShiftModifier)!=0 && current==Qt::Key_C)
	{
		fstream file;
		file.open("Camera.txt", ios_base::in);
		float x, y, z, w;
		file >> x; file >> y; file >> z; file >> w;
		glm::fquat q(w, x, y, z);
		cam.setQuaternion(q);
		file >> x; file >> y; file >> z;
		cam.setTranslation(glm::vec3(x, y, z));
		file >> x;
		cam.setDist(x);
		m_scene.draw_scene();
		updateGL(); PRINT_OPENGL_ERROR();
		cout << "Position de la caméra chargée" << endl;
	}
	if (current==Qt::Key_X)
	{
		cam.alignX();
		updateGL(); PRINT_OPENGL_ERROR();
	}
	if (current==Qt::Key_Y)
	{
		cam.alignY();
		updateGL(); PRINT_OPENGL_ERROR();
	}
	if (current==Qt::Key_Z)
	{
		cam.alignZ();
		updateGL(); PRINT_OPENGL_ERROR();
	}
	if (current==Qt::Key_Up)
	{
		m_scene.increasePtSize();
		updateGL(); PRINT_OPENGL_ERROR();
	}
	if (current==Qt::Key_Down)
	{
		m_scene.decreasePtSize();
		updateGL(); PRINT_OPENGL_ERROR();
	}
	if (current==Qt::Key_P)
	{
		m_scene.projectPoints();
		updateGL(); PRINT_OPENGL_ERROR();
	}
	if (current==Qt::Key_1)
		m_scene.renderPointSet = !m_scene.renderPointSet;
	if (current==Qt::Key_2)
		m_scene.renderPointsToProject = !m_scene.renderPointsToProject;
	if (current==Qt::Key_3)
		m_scene.renderProjectedPoints = !m_scene.renderProjectedPoints;
	if (current==Qt::Key_4)
		m_scene.renderDualContour = !m_scene.renderDualContour;
	if (current==Qt::Key_5)
		m_scene.renderInvalidPoints = !m_scene.renderInvalidPoints;
	if (current==Qt::Key_6)
		m_scene.renderComparison = !m_scene.renderComparison;

	if ((mod&Qt::ShiftModifier)==0 && current==Qt::Key_D)// Without shift, do not record octree
	{
		killTimer(m_mainTimer);
		string filename = QFileDialog::getSaveFileName(this, "Choose a file to save the dual contour", "../data/", "ply files (*.ply)", nullptr, QFileDialog::DontUseNativeDialog).toStdString();
		m_mainTimer = startTimer(0);
		if (filename == "")
		{
			cerr << "Error, no file selected" << endl;
			return;
		}
		m_scene.computeDualContouring(filename, false);
	}

	if ((mod&Qt::ShiftModifier)!=0 && current==Qt::Key_D)// With shift, record octree
	{
		killTimer(m_mainTimer);
		string filename = QFileDialog::getSaveFileName(this, "Choose a file to save the dual contour", "../data/", "ply files (*.ply)", nullptr, QFileDialog::DontUseNativeDialog).toStdString();
		m_mainTimer = startTimer(0);
		if (filename == "")
		{
			cerr << "Error, no file selected" << endl;
			return;
		}
		m_scene.computeDualContouring(filename, true);
	}

	if (current==Qt::Key_W)
		m_scene.wireframe = !m_scene.wireframe;

	if (current==Qt::Key_M)
	{
		stopTimer();
		m_scene.mortonTest();
		restartTimer();
	}

	if(current == Qt::Key_O)
	{
		m_scene.save_octree("CompleteOctree.ply");
		for(unsigned int d = 0; d <= 12; d++)
			m_scene.save_octree("OctreeDepth-" + std::to_string(d) + ".ply", d, d);
	}

	if (current == Qt::Key_T)
	{
		m_scene.resetOctree();
	}

	if (current == Qt::Key_I)
	{
		m_scene.count_nodes_on_projected_points();
	}

	if (current == Qt::Key_R)
	{
		m_scene.compileShaders();
		cout << "Shaders compiled" << endl;
	}

	QGLWidget::keyPressEvent(event);
}

void WidgetOpenGL::timerEvent(QTimerEvent *event)
{
	event->accept();
	update(); PRINT_OPENGL_ERROR();
}


void WidgetOpenGL::mousePressEvent(QMouseEvent *event)
{
	int const x = event->x();
	int const y = event->y();
	cam.xPrevious()= x;
	cam.yPrevious()= y;
	Qt::KeyboardModifiers mod=event->modifiers();
	if (m_planeMode && (mod&Qt::AltModifier))//Plane selection
	{
		glm::vec3 pt1 = cam.get3DPointFromClick(x, y, -1.f);
		glm::vec3 pt2 = cam.get3DPointFromClick(x, y);
		m_scene.selectPtsCloseToPlane(glm::normalize(glm::cross(glm::normalize(cam.getUpVector()), glm::normalize(pt2 - pt1))), -pt1, 0.002f);
	}
	else if (mod&Qt::ControlModifier)//Selection of closest point to click
	{
		glm::vec3 pt1 = cam.get3DPointFromClick(x, y, -1.f);
		glm::vec3 pt2 = cam.get3DPointFromClick(x, y);
		m_scene.selectClosestPoint(pt1, pt2, 0.002f);
	}
}

void WidgetOpenGL::mouseMoveEvent(QMouseEvent *event)
{
	int const x=event->x();
	int const y=event->y();

	int const ctrl_pressed  = (event->modifiers() & Qt::ControlModifier);
	int const shift_pressed = (event->modifiers() & Qt::ShiftModifier);

	if(!ctrl_pressed && !shift_pressed && (event->buttons() & Qt::LeftButton))
		cam.rotation(x, y);
	if(!ctrl_pressed && !shift_pressed && (event->buttons() & Qt::RightButton))
	{
		cam.zoom(y);
		m_steps = y;
	}

	// Shift+Left button controls the window translation (left/right, bottom/up)
	float const dL=0.0001f*(1+10*fabs(cam.getDist()));
	if( (!ctrl_pressed && shift_pressed && (event->buttons() & Qt::LeftButton)) || (event->buttons() & Qt::MiddleButton))
	{
		cam.goUp(dL*(y-cam.yPrevious()));
		cam.goRight(-dL*(x-cam.xPrevious()));
	}

	// Shift+Right button enables to translate forward/backward
	if( !ctrl_pressed && shift_pressed && (event->buttons() & Qt::RightButton) )
		cam.goForward(5.0f*dL*(y-cam.yPrevious()));

	cam.xPrevious()=x;
	cam.yPrevious()=y;
}

void WidgetOpenGL::wheelEvent(QWheelEvent * event)
{
	int numSteps = event->angleDelta().y() / 8;

	m_steps += numSteps;
	cam.zoomWheel(m_steps);
	event->accept();
}

double WidgetOpenGL::drawTime()
{
	return m_scene.drawTime();
}

void WidgetOpenGL::updateUI()
{
	emit update_UI();
}

void WidgetOpenGL::connectKernel()
{
	emit connect_kernel();
}
