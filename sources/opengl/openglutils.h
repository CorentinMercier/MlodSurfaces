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
#ifndef OPENGLUTILS_H
#define OPENGLUTILS_H

#if _WIN32
#define NOMINMAX
#include <Windows.h>
#endif
#include <glad/glad.h>
#include <string>
#include <iostream>
#include <fstream>
#include <streambuf>
#include <algorithm>

/** Macro indicating OpenGL errors (file and line) */
#define PRINT_OPENGL_ERROR() print_opengl_error(__FILE__, __LINE__)

typedef struct
{
	GLenum type;
	const char* filename;
	GLuint shader;
}ShaderInfo;

//Function returning the ID of the program using the shaders included in ShaderInfo
GLuint loadShaders(ShaderInfo* shaders);

/** Print OpenGL Error given the information of the file and the line
 *  Function called by the macro PRINT_OPENGL_ERROR */
bool print_opengl_error(const char *file, int line);

GLint get_uni_loc(GLuint program, const GLchar *name);

class Timer {

public:
	void start();
	void stop();
	double getTime(){return time;} //time in ms

private:
	bool queryGenerated = false;
	GLuint queryId[2] = {0, 0};
	double time = 0.;

};

#endif // OPENGLUTILS_H
