#pragma once
#include "glError.h"
#include "gl/glew.h"
#include "gl/freeglut.h"
#include "polyline.h"

// Draws a dot at the specified location
// x, y, z are the coordinates of the dot
// radius: radius of the dot
// R: red channel for the dot color [0,1]
// B: blue channel for the dot color [0,1]
// G: green channel of the dot color [0,1]
void drawDot(double x, double y, double z, double radius = 0.15, float R = 0.0, float G = 0.0, float B = 0.0)
{
	glDisable(GL_POLYGON_OFFSET_FILL);
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glMatrixMode(GL_MODELVIEW);

	CHECK_GL_ERROR();

	GLfloat mat_diffuse[4];
	mat_diffuse[0] = R;
	mat_diffuse[1] = G;
	mat_diffuse[2] = B;
	mat_diffuse[3] = 1.0;

	CHECK_GL_ERROR();

	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);

	CHECK_GL_ERROR();

	GLUquadric* quadric = gluNewQuadric();
	glPushMatrix();
	glTranslated(x, y, z);
	glColor3f(R, G, B);
	gluSphere(quadric, radius, 16, 16);
	glPopMatrix();
	gluDeleteQuadric(quadric);
}

// Draws a single line segment (LineSegment defined in polyline.h file)
// width: width of the line segment
// R: red channel for the line color [0,1]
// B: blue channel for the line color [0,1]
// G: green channel of the line color [0,1]
void drawLineSegment(LineSegment ls, double width = 1.0, float R = 0.0, float G = 0.0, float B = 0.0)
{
	glDisable(GL_LIGHTING);
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glLineWidth(width);

	glBegin(GL_LINES);
	glColor3f(R, G, B);
	glVertex3f(ls.start.x, ls.start.y, ls.start.z);
	glVertex3f(ls.end.x, ls.end.y, ls.end.z);
	glEnd();
}

// Draws a polyline (PolyLine defined in polyline.h file)
// width: width of the polyline
// R: red channel for the polyline color [0,1]
// B: blue channel for the polyline color [0,1]
// G: green channel of the polyline color [0,1]
void drawPolyLine(PolyLine pl, double width = 1.0, float R = 0.0, float G = 0.0, float B = 0.0)
{
	glDisable(GL_LIGHTING);
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glLineWidth(width);

	glBegin(GL_LINES);
	glColor3f(R, G, B);

	for (int i = 0; i < pl.size(); i++)
	{
		glVertex3f(pl[i].start.x, pl[i].start.y, pl[i].start.z);
		glVertex3f(pl[i].end.x, pl[i].end.y, pl[i].end.z);
	}
	
	glEnd();
}

// example function for using dots and polylines
void dots_and_lines_example(std::vector<icVector3>* points, std::vector<PolyLine>* lines)
{

	// make polylines for linear, quadratic, and cubic functions
	PolyLine linear, quadratic, cubic;
	for (int x = -10; x < 10; x++)
	{
		double y_linear = (double)x;
		double y_quadratic = (double)x * (double)x / 10.0;
		double y_cubic = (double)x * (double)x * (double)x / 100.0;

		double x1 = x + 1;
		double y1_linear = (double)x1;
		double y1_quadratic = (double)x1 * (double)x1 / 10.0;
		double y1_cubic = (double)x1 * (double)x1 * (double)x1 / 100.0;

		LineSegment linear_seg = LineSegment(x, y_linear, 0, x1, y1_linear, 0);
		LineSegment quadratic_seg = LineSegment(x, y_quadratic, 0, x1, y1_quadratic, 0);
		LineSegment cubic_seg = LineSegment(x, y_cubic, 0, x1, y1_cubic, 0);

		linear.push_back(linear_seg);
		quadratic.push_back(quadratic_seg);
		cubic.push_back(cubic_seg);
	}
	
	lines->push_back(linear);
	lines->push_back(quadratic);
	lines->push_back(cubic);

	// make dots along x and y axes
	for (int i = -10; i <= 10; i++)
	{
		icVector3 x_ax = icVector3(i, 0, 0);
		icVector3 y_ax = icVector3(0, i, 0);
		points->push_back(x_ax);
		points->push_back(y_ax);
	}
}
