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

#pragma once

#ifndef WIDGETOPENGL_H
#define WIDGETOPENGL_H

#if _WIN32
#define NOMINMAX
#include <Windows.h>
#endif
#include <glad/glad.h>
#include <GL/gl.h>

#include <QtOpenGL/QGLWidget>
#include <QMouseEvent>
#include <QFileDialog>
#include <QDialogButtonBox>

#include <cstring>
#include <iostream>

#include "../scene/scene.h"
#include "../opengl/openglutils.h"
#include "../opengl/camera.h"


struct CommandArgs
{
	string input_file = "";
	string output_folder = "";
	string output_file = "";
	Kernel kernel = Kernel();
	unsigned depth = 6;
	unsigned nb_knn = 20;
	float ball_radius = 0.1;
	float cutoff = 2.0;
	bool cpu = false;
	bool knn = false;
	bool ball = false;
	bool dc = false;
	bool automatic = false;
	bool fast = true;
	bool inputOctree = false;
	bool customKernel = false;
};

using namespace std;

/** Qt Widget to render OpenGL scene */
class WidgetOpenGL : public QGLWidget
{
	Q_OBJECT

public:

	WidgetOpenGL(const QGLFormat& format, const CommandArgs& cmd, QGLWidget *parent = 0);
	~WidgetOpenGL();

	bool openNewScene();

	scene& get_scene();
	//Camera
	camera cam;

	double drawTime();

	void changePlaneMode(){m_planeMode = !m_planeMode; m_scene.renderSelection = m_planeMode;}
	void selectionNotEmpty(){emit thereIsASelection();}
	void selectionEmpty(){emit thereIsNoSelection();}

	void stopTimer(){killTimer(m_mainTimer);}
	void restartTimer(){m_mainTimer = startTimer(0);}

	string getFilename(){return m_filename;}

	void changeBackground(bool whiteBackground){m_whiteBackground = whiteBackground;}

	void changeRadiusFactor(float factor){m_radiusFactor = factor;}
	float getRadiusFactor(){return m_radiusFactor;}

	void updateUI();
	void connectKernel();

	CommandArgs& get_cmd(){return m_cmd;}

signals:
	void thereIsASelection();
	void thereIsNoSelection();
	void scene_loaded();
	void update_UI();
	void connect_kernel();

protected:

	/** Setup the OpenGL rendering mode */
	void initializeGL();
	/** The actual rendering function */
	void paintGL();
	/** Function called when the window is resized */
	void resizeGL(const int width, const int height);

	/** Function called a button of the mouse is pressed */
	void mousePressEvent(QMouseEvent *event);
	/** Function called when the wheel of the mouse is moved **/
	void wheelEvent(QWheelEvent *event);
	/** Function called when the mouse is moved */
	void mouseMoveEvent(QMouseEvent *event);
	/** Function called in a timer loop */
	void timerEvent(QTimerEvent *event);
	/** Function called when keyboard is pressed */
	void keyPressEvent(QKeyEvent *event);

private:
	/** Init the OpenGL rendering mode once at the beginning */
	void setup_opengl();
	/** Init GLAD once at the beginning */
	void setup_glad();
	/** Print on the command line the actual version of the OpenGL Context */
	void print_current_opengl_context() const;
	/** All the content of the 3D scene */
	scene m_scene;

	string m_filename;
	int m_steps = -1325;

	int m_mainTimer;

	bool m_planeMode = false;
	bool m_whiteBackground = true;

	float m_radiusFactor = 1.f;

	CommandArgs m_cmd;
};

#endif // WIDGETOPENGL_H
