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

#include "openglutils.h"
#if _WIN32
#include <windows.h>
#endif
#include <GL/glu.h>
#include <GL/gl.h>


static const GLchar* readShader(const char* filename)
{
	FILE* file=fopen(filename, "rb");
	if (!file)
		return nullptr;
	fseek( file, 0, SEEK_END );
	unsigned len = unsigned(ftell( file ));
	fseek( file, 0, SEEK_SET );
	GLchar* content=new GLchar[len+1];
	if (fread(content, 1, len, file)==0)
	{
		std::cout << "Error while reading shaders" << std::endl;
	}
	fclose(file);
	content[len]=0;
	return const_cast<const GLchar*>(content);
}


GLuint loadShaders(ShaderInfo *shaders)
{
	if (shaders == nullptr)
	{
		std::cout << "No shaders found" << std::endl;
		return 0;
	}
	GLuint program = glCreateProgram();

	ShaderInfo* entry = shaders;
	while(entry->type!=GL_NONE)
	{
		GLuint shader=glCreateShader(entry->type);
		const GLchar* source=readShader(entry->filename);
		if (source==nullptr)
		{
			std::cout << "No content in shader " << entry->filename << std::endl;
			for (entry=shaders; entry->type!=GL_NONE; ++entry)
			{
				glDeleteShader(entry->shader);
				entry->shader=0;
			}
			return 0;
		}
		glShaderSource(shader, 1, &source, nullptr);
		delete[] source;
		glCompileShader(shader);
		GLint compiled;
		glGetShaderiv(shader, GL_COMPILE_STATUS, &compiled);
		if(!compiled)
		{
			GLsizei len;
			glGetShaderiv( shader, GL_INFO_LOG_LENGTH, &len );

			GLchar* log = new GLchar[unsigned(len)+1];
			glGetShaderInfoLog( shader, len, &len, log );
			std::cerr << "Shader compilation failed: " << log << std::endl;
			delete [] log;
			return 0;
		}
		glAttachShader(program, shader);
		++entry;
	}
	glLinkProgram(program);
	GLint linked;
	glGetProgramiv(program, GL_LINK_STATUS, &linked);
	if(!linked)
	{
		GLsizei len;
		glGetProgramiv(program, GL_INFO_LOG_LENGTH, &len);
		GLchar* log = new GLchar[unsigned(len+1)];
		glGetProgramInfoLog( program, len, &len, log );
		std::cerr << "Program linkage failed: " << log << std::endl;
		delete [] log;

		for ( entry = shaders; entry->type != GL_NONE; ++entry ) {
			glDeleteShader( entry->shader );
			entry->shader = 0;
		}
		std::cout << "Shader program linkage error" << std::endl;
		return 0;
	}
	return program;
}

/*****************************************************************************\
 * print_opengl_error                                                        *
\*****************************************************************************/
bool print_opengl_error(const char *file, int line)
{
	// Returns true if an OpenGL error occurred, false otherwise.
	GLenum error;
	bool   ret_code = false;

	error = glGetError();
	while (error != GL_NO_ERROR)
	{
		std::cout << "glError in file " << file << ", line " << line << ": " << gluErrorString(error) << std::endl;
		ret_code = true;
		error = glGetError();
	}
	return ret_code;
}

GLint get_uni_loc(GLuint program, const GLchar *name)
{
  GLint loc;

  loc = glGetUniformLocation(program, name); PRINT_OPENGL_ERROR();
  if (loc == -1)
	std::cerr << "No such uniform named \"" << name << "\"" << std::endl;

  return loc;
}

void Timer::start(){
	if(!queryGenerated){
		glGenQueries(2, queryId);
		queryGenerated = true;
	}
	glQueryCounter(queryId[0], GL_TIMESTAMP);
}

void Timer::stop(){
	glQueryCounter(queryId[1], GL_TIMESTAMP);
	GLuint64 startTime, stopTime;
	glGetQueryObjectui64v(queryId[0], GL_QUERY_RESULT, &startTime);
	glGetQueryObjectui64v(queryId[1], GL_QUERY_RESULT, &stopTime);
	time = (stopTime - startTime) / 1000000.;
}
