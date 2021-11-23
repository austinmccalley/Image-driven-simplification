#ifndef INC_GLERROR_H
#define INC_GLERROR_H

#include "gl/glew.h"
#include <stdio.h>

static int checkGLError(char *file, int line){
	GLenum glError;
	int returnCode = 0;

	glError = glGetError();
	while (glError != GL_NO_ERROR) 
	{
		printf("GL Error #%d(%s) in File %s at line: %d\n", glError, gluErrorString(glError), file, line);
		returnCode = 1;
		glError = glGetError();
	}
	return returnCode;
}
#define CHECK_GL_ERROR() checkGLError(__FILE__, __LINE__)


#endif

