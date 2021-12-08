#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <iostream>

#include "glError.h"
#include "gl/glew.h"
#include "gl/freeglut.h"
#include "ply.h"
#include "icVector.H"
#include "icMatrix.H"
#include "polyhedron.h"
#include "trackball.h"
#include "tmatrix.h"

#include "drawUtil.h"

Polyhedron* poly;
std::vector<PolyLine> lines;
std::vector<icVector3> points;

/*scene related variables*/
const float zoomspeed = 0.9;
int win_width = 1024;
int win_height = 1024;
float aspectRatio = win_width / win_height;
const int view_mode = 0;		// 0 = othogonal, 1=perspective
const double radius_factor = 0.9;

/*
Use keys 1 to 0 to switch among different display modes.
Each display mode can be designed to show one type
visualization result.

Predefined ones:
display mode 1: solid rendering
display mode 2: show wire frames
display mode 3: render each quad with colors of vertices
*/
int display_mode = 1;

/*User Interaction related variables*/
float s_old, t_old;
float rotmat[4][4];
double zoom = 1.0;
double translation[2] = { 0, 0 };
int mouse_mode = -2;	// -1 = no action, 1 = translate y, 2 = rotate

int N = 50;
const double STEP = 1e-2;
const int STEP_MAX = 10000;

std::vector<PolyLine> streamlines;
std::vector<icVector3> sources;
std::vector<icVector3> sinks;
std::vector<icVector3> saddles;
std::vector<icVector3> higher_order;

int SCALE_MIN = 0;
int SCALE_MAX = 25;

// IBFV related variables (Van Wijk 2002)
//https://www.win.tue.nl/~vanwijk/ibfv/
#define NPN		64
#define SCALE	4.0
#define ALPHA	8
float tmax = win_width / (SCALE * NPN);
float dmax = SCALE / win_width;
unsigned char* pixels;

/******************************************************************************
Forward declaration of functions
******************************************************************************/

void init(void);
void initIBFV();

/*glut attaching functions*/
void keyboard(unsigned char key, int x, int y);
void motion(int x, int y);
void displayIBFV();
void display(void);
void mouse(int button, int state, int x, int y);
void mousewheel(int wheel, int direction, int x, int y);
void reshape(int width, int height);

/*functions for element picking*/
void display_vertices(GLenum mode, Polyhedron* poly);
void display_quads(GLenum mode, Polyhedron* poly);
void display_selected_vertex(Polyhedron* poly);
void display_selected_quad(Polyhedron* poly);

/*display vis results*/
void display_polyhedron(Polyhedron* poly);

/******************************************************************************
Main program.
******************************************************************************/

int main(int argc, char* argv[])
{
	/*load mesh from ply file*/
	FILE* this_file = fopen("../data/vector_data/v6.ply", "r");
	poly = new Polyhedron(this_file);
	fclose(this_file);

	/*initialize the mesh*/
	poly->initialize(); // initialize the mesh
	poly->write_info();


	/*init glut and create window*/
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowPosition(20, 20);
	glutInitWindowSize(win_width, win_height);
	glutCreateWindow("Scientific Visualization");


	/*initialize openGL*/
	init();

	/*the render function and callback registration*/
	glutKeyboardFunc(keyboard);
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutMotionFunc(motion);
	glutMouseFunc(mouse);
	glutMouseWheelFunc(mousewheel);

	/*event processing loop*/
	glutMainLoop();

	/*clear memory before exit*/
	poly->finalize();	// finalize everything
	free(pixels);
	return 0;
}

/******************************************************************************
Set projection mode
******************************************************************************/

void set_view(GLenum mode)
{
	GLfloat light_ambient0[] = { 0.3, 0.3, 0.3, 1.0 };
	GLfloat light_diffuse0[] = { 0.7, 0.7, 0.7, 1.0 };
	GLfloat light_specular0[] = { 0.0, 0.0, 0.0, 1.0 };

	GLfloat light_ambient1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_diffuse1[] = { 0.5, 0.5, 0.5, 1.0 };
	GLfloat light_specular1[] = { 0.0, 0.0, 0.0, 1.0 };

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);

	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);


	glMatrixMode(GL_PROJECTION);
	if (mode == GL_RENDER)
		glLoadIdentity();

	if (aspectRatio >= 1.0) {
		if (view_mode == 0)
			glOrtho(-radius_factor * zoom * aspectRatio, radius_factor * zoom * aspectRatio, -radius_factor * zoom, radius_factor * zoom, -1000, 1000);
		else
			glFrustum(-radius_factor * zoom * aspectRatio, radius_factor * zoom * aspectRatio, -radius_factor * zoom, radius_factor * zoom, 0.1, 1000);
	}
	else {
		if (view_mode == 0)
			glOrtho(-radius_factor * zoom, radius_factor * zoom, -radius_factor * zoom / aspectRatio, radius_factor * zoom / aspectRatio, -1000, 1000);
		else
			glFrustum(-radius_factor * zoom, radius_factor * zoom, -radius_factor * zoom / aspectRatio, radius_factor * zoom / aspectRatio, 0.1, 1000);
	}

	GLfloat light_position[3];
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	light_position[0] = 5.5;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	light_position[0] = -0.1;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT2, GL_POSITION, light_position);
}

/******************************************************************************
Update the scene
******************************************************************************/

void set_scene(GLenum mode, Polyhedron* poly)
{
	glTranslatef(translation[0], translation[1], -3.0);

	/*multiply rotmat to current mat*/
	{
		int i, j, index = 0;

		GLfloat mat[16];

		for (i = 0; i < 4; i++)
			for (j = 0; j < 4; j++)
				mat[index++] = rotmat[i][j];

		glMultMatrixf(mat);
	}

	glScalef(0.9 / poly->radius, 0.9 / poly->radius, 0.9 / poly->radius);
	glTranslatef(-poly->center.entry[0], -poly->center.entry[1], -poly->center.entry[2]);
}

/******************************************************************************
Init scene
******************************************************************************/

void init(void) {

	mat_ident(rotmat);

	/* select clearing color */
	glClearColor(0.0, 0.0, 0.0, 0.0);  // background
	glShadeModel(GL_FLAT);
	glPolygonMode(GL_FRONT, GL_FILL);

	glDisable(GL_DITHER);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	//set pixel storage modes
	glPixelStorei(GL_PACK_ALIGNMENT, 1);

	glEnable(GL_NORMALIZE);
	if (poly->orientation == 0)
		glFrontFace(GL_CW);
	else
		glFrontFace(GL_CCW);
}

/******************************************************************************
Initialize IBFV patterns
******************************************************************************/

void initIBFV()
{
	pixels = (unsigned char*)malloc(sizeof(unsigned char) * win_width * win_height * 3);
	memset(pixels, 255, sizeof(unsigned char) * win_width * win_height * 3);

	tmax = win_width / (SCALE * NPN);
	dmax = SCALE / win_width;

	int lut[256];
	int phase[NPN][NPN];
	GLubyte pat[NPN][NPN][4];
	int i, j, k;

	for (i = 0; i < 256; i++) lut[i] = i < 127 ? 0 : 255;
	for (i = 0; i < NPN; i++)
		for (j = 0; j < NPN; j++) phase[i][j] = rand() % 256;

	for (i = 0; i < NPN; i++)
	{
		for (j = 0; j < NPN; j++)
		{
			pat[i][j][0] =
				pat[i][j][1] =
				pat[i][j][2] = lut[(phase[i][j]) % 255];
			pat[i][j][3] = ALPHA;
		}
	}

	glNewList(1, GL_COMPILE);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0, GL_RGBA, GL_UNSIGNED_BYTE, pat);
	glEndList();
}

/******************************************************************************
Pick objects from the scene
******************************************************************************/

int processHits(GLint hits, GLuint buffer[])
{
	unsigned int i, j;
	GLuint names, * ptr;
	double smallest_depth = 1.0e+20, current_depth;
	int seed_id = -1;
	unsigned char need_to_update;

	ptr = (GLuint*)buffer;
	for (i = 0; i < hits; i++) {  /* for each hit  */
		need_to_update = 0;
		names = *ptr;
		ptr++;

		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		for (j = 0; j < names; j++) {  /* for each name */
			if (need_to_update == 1)
				seed_id = *ptr - 1;
			ptr++;
		}
	}
	return seed_id;
}

/******************************************************************************
Diaplay all quads for selection
******************************************************************************/

void display_quads(GLenum mode, Polyhedron* this_poly)
{
	unsigned int i, j;
	GLfloat mat_diffuse[4];

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glDisable(GL_LIGHTING);

	for (i = 0; i < this_poly->nquads; i++) {
		if (mode == GL_SELECT)
			glLoadName(i + 1);

		Quad* temp_q = this_poly->qlist[i];

		glBegin(GL_POLYGON);
		for (j = 0; j < 4; j++) {
			Vertex* temp_v = temp_q->verts[j];
			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
		}
		glEnd();
	}
}

/******************************************************************************
Diaplay all vertices for selection
******************************************************************************/

void display_vertices(GLenum mode, Polyhedron* this_poly)
{
	for (int i = 0; i < this_poly->nverts; i++) {
		if (mode == GL_SELECT)
			glLoadName(i + 1);

		CHECK_GL_ERROR();

		Vertex* temp_v = this_poly->vlist[i];
		drawDot(temp_v->x, temp_v->y, temp_v->z, 0.15);
	}
	CHECK_GL_ERROR();
}

/******************************************************************************
Diaplay selected quad
******************************************************************************/

void display_selected_quad(Polyhedron* this_poly)
{
	if (this_poly->selected_quad == -1)
	{
		return;
	}

	unsigned int i, j;

	glDisable(GL_POLYGON_OFFSET_FILL);
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glDisable(GL_LIGHTING);

	Quad* temp_q = this_poly->qlist[this_poly->selected_quad];

	glBegin(GL_POLYGON);
	for (j = 0; j < 4; j++) {
		Vertex* temp_v = temp_q->verts[j];
		glColor3f(1.0, 0.0, 0.0);
		glVertex3d(temp_v->x, temp_v->y, temp_v->z);
	}
	glEnd();
}

/******************************************************************************
Diaplay selected vertex
******************************************************************************/

void display_selected_vertex(Polyhedron* this_poly)
{
	if (this_poly->selected_vertex == -1)
	{
		return;
	}

	Vertex* temp_v = this_poly->vlist[this_poly->selected_vertex];
	drawDot(temp_v->x, temp_v->y, temp_v->z, 0.15, 1.0, 0.0, 0.0);

	CHECK_GL_ERROR();
}



void classifyVector(double s) {
	for (int i = 0; i < poly->nverts; i++) {
		Vertex* vTemp = poly->vlist[i];

		if (vTemp->scalar >= s) {
			vTemp->level = 1;
		}
		else {
			vTemp->level = 0;
		}
	}
}

void calcCrossingPoints(double s, int z) {

	/* Calculate crossing points */
	for (int i = 0; i < poly->nedges; i++) {
		Edge* eTemp = poly->elist[i];
		Vertex* v1 = eTemp->verts[0];
		Vertex* v2 = eTemp->verts[1];

		int a = SCALE_MIN;
		int b = SCALE_MAX;
		int dba = b - a;
		int scalarDiff = poly->maxScalar - poly->minScalar;

		if (z == 1) {
			v1->z = dba * ((v1->scalar - poly->minScalar) / scalarDiff) + a;
			v2->z = dba * ((v2->scalar - poly->minScalar) / scalarDiff) + a;
		}

		/* Delete any old crossing points */
		if (eTemp->crossing) {
			delete eTemp->crossing;
			eTemp->crossing = NULL;
		}

		if (v1->level != v2->level) {
			/* Find cross point using linear interpolation */

			double s1 = v1->scalar;
			double s2 = v2->scalar;

			/*
			The s is the scalar value of the contour you are trying to find, and scalar value at the crossing point should be equal to s.
			Alpha is just the interpolation factor that tells you how far along the edge the crossing is from v1
			*/
			double alpha = (s - s1) / (s2 - s1);

			double vx = alpha * (v2->x - v1->x) + v1->x;
			double vy = alpha * (v2->y - v1->y) + v1->y;
			double vz = alpha * (v2->z - v1->z) + v1->z;


			icVector3* cp = new icVector3(vx, vy, vz);
			eTemp->crossing = cp;
		}
	}
}


void connectCP(double s) {
	for (int i = 0; i < poly->nquads; i++) {
		PolyLine plCountour;
		Quad* qTemp = poly->qlist[i];

		std::vector<icVector3*> crossPoints;

		/* Find saved crossing points on each edge */
		for (int j = 0; j < 4; j++) {
			Edge* eTemp = qTemp->edges[j];

			if (eTemp->crossing != NULL)
				crossPoints.push_back(eTemp->crossing);
		}


		switch (crossPoints.size())
		{

		case 2: {

			/* If there are only two crossing points we only compare crossPoints[0] and crossPoints[1] */

			icVector3* cp1 = crossPoints[0];
			icVector3* cp2 = crossPoints[1];

			LineSegment line(cp1->x, cp1->y, cp1->z, cp2->x, cp2->y, cp2->z);
			plCountour.push_back(line);

			break;
		}
		case 4: {

			/* If there are 4 crossing points we compare all crossPoints[0...3]*/

			double closestDelta = INT_MAX;
			icVector3* p1 = NULL;
			icVector3* p2 = NULL;
			int p1I = -1;
			int p2I = -1;

			for (int j = 0; j < 4; j++) {
				for (int k = 0; k < 4; k++) {
					if (j != k) {
						icVector3* cp1 = crossPoints[j];
						icVector3* cp2 = crossPoints[k];

						if (dot(*cp1, *cp2) == 0) {
							std::cout << "HIT INTERSECTION" << std::endl;
							continue;
						}

						LineSegment line(cp1->x, cp1->y, cp1->z, cp2->x, cp2->y, cp2->z);

						icVector3 mp = line.midpoint();

						double scalar = length(mp);

						double delta = s - scalar;

						if (delta < closestDelta) {
							closestDelta = delta;
							p1 = cp1;
							p2 = cp2;
						}
						else {
							p1I = j;
							p2I = k;
						}
					}
				}
			}

			if (p1 && p2) {
				std::cout << "Case hit for 4 crossing points" << std::endl;
				std::cout << p1I << " " << p2I << std::endl;

				icVector3* cp3 = crossPoints[p1I];
				icVector3* cp4 = crossPoints[p2I];

				LineSegment line2(cp3->x, cp3->y, cp3->z, cp4->x, cp4->y, cp4->z);



				LineSegment line(p1->x, p1->y, p1->z, p2->x, p2->y, p2->z);
				plCountour.push_back(line);
				plCountour.push_back(line2);
			}
			break;
		}
		default:
			break;
		}

		lines.push_back(plCountour);
	}
}

bool between10(double x) {
	return x >= 0 && x <= 1;
}

double sing_prox(icVector2 pos) {
	double prox = DBL_MAX;

	for (int i = 0; i < sources.size(); i++) {
		icVector3 sing = sources[i];
		icVector2 spos = icVector2(sing.x, sing.y);

		double dist = length(pos - spos);
		if (dist < prox)
			prox = dist;
	}

	for (int i = 0; i < sinks.size(); i++)
	{
		icVector3 sing = sinks[i];
		icVector2 spos = icVector2(sing.x, sing.y);

		double dist = length(pos - spos);
		if (dist < prox)
			prox = dist;
	}

	for (int i = 0; i < higher_order.size(); i++) {
		icVector3 sing = higher_order[i];
		icVector2 spos = icVector2(sing.x, sing.y);

		double dist = length(pos - spos);
		if (dist < prox)
			prox = dist;
	}

	for (int i = 0; i < saddles.size(); i++) {
		icVector3 sing = saddles[i];
		icVector2 spos = icVector2(sing.x, sing.y);

		double dist = length(pos - spos);
		if (dist < prox)
			prox = dist;
	}

	return prox;
}

void find_singularities() {
	for (int i = 0; i < poly->nquads; i++) {
		Quad* qTemp = poly->qlist[i];

		double x1, x2, y1, y2, fx1y1, fx2y1, fx1y2, fx2y2, gx1y1, gx2y1, gx1y2, gx2y2;

		for (int j = 0; j < 4; j++) {
			Vertex* vTemp = qTemp->verts[j];

			x1 = INT_MAX;
			x2 = INT_MIN;
			y1 = INT_MAX;
			y2 = INT_MIN;

			for (int j = 0; j < 4; j++)
			{
				Vertex* v = qTemp->verts[j];
				// Current x is smaller than x1
				if (v->x < x1) {
					x1 = v->x;
				}

				// Current y is smaller than y1
				if (v->y < y1) {
					y1 = v->y;
				}

				// Current x is bigger than x2
				if (v->x > x2) {
					x2 = v->x;
				}

				// Current y is bigger than y2 
				if (v->y > y2) {
					y2 = v->y;
				}
			}

			for (int j = 0; j < 4; j++) {
				Vertex* v = qTemp->verts[j];
				if (v->x == x1 && v->y == y1) {
					fx1y1 = v->vx;
					gx1y1 = v->vy;
				}
				if (v->x == x2 && v->y == y2) {
					fx2y2 = v->vx;
					gx2y2 = v->vy;
				}
				if (v->x == x1 && v->y == y2) {
					fx1y2 = v->vx;
					gx1y2 = v->vy;
				}
				if (v->x == x2 && v->y == y1) {
					fx2y1 = v->vx;
					gx2y1 = v->vy;
				}
			}
		}

		double a00 = fx1y1;
		double a10 = fx2y1 - fx1y1;
		double a01 = fx1y2 - fx1y1;
		double a11 = fx1y1 - fx2y1 - fx1y2 + fx2y2;
		double b00 = gx1y1;
		double b10 = gx2y1 - gx1y1;
		double b01 = gx1y2 - gx1y1;
		double b11 = gx1y1 - gx2y1 - gx1y2 + gx2y2;
		double c00 = a11 * b00 - a00 * b11;
		double c10 = a11 * b10 - a10 * b11;
		double c01 = a11 * b01 - a01 * b11;

		double a = -a11 * c10;
		double b = -a11 * c00 - a01 * c10 + a10 * c01;
		double c = a00 * c01 - a01 * c00;


		double innerSqrt = (pow(b, 2)) - (4 * a * c);
		double bottomS = 2 * a;

		if (innerSqrt < 0) {
			continue;
		}
		else if (bottomS == 0) {
			continue;
		}

		double sNeg = (-b - sqrt(innerSqrt)) / bottomS;
		double sPos = (-b + sqrt(innerSqrt)) / bottomS;


		double tNeg = -(c00 / c01) - (c10 / c01) * sNeg;
		double tPos = -(c00 / c01) - (c10 / c01) * sPos;

		if (between10(sNeg) && between10(tNeg)) {
			icVector2* sing = new icVector2(x1 + sNeg * (x2 - x1), y1 + tNeg * (y2 - y1));
			qTemp->singularity = sing;
		}
		else if (between10(sPos) && between10(tPos)) {
			icVector2* sing = new icVector2(x1 + sPos * (x2 - x1), y1 + tPos * (y2 - y1));
			qTemp->singularity = sing;
		}
		else {
			continue;
		}

	}
}

void classify_singularities() {
	sources.clear();
	saddles.clear();
	higher_order.clear();
	for (int i = 0; i < poly->nquads; i++)
	{
		Quad* qTemp = poly->qlist[i];

		if (qTemp->singularity != NULL) {
			double x1, x2, y1, y2, fx1y1, fx2y1, fx1y2, fx2y2, gx1y1, gx2y1, gx1y2, gx2y2;

			for (int j = 0; j < 4; j++) {
				Vertex* vTemp = qTemp->verts[j];

				x1 = INT_MAX;
				x2 = INT_MIN;
				y1 = INT_MAX;
				y2 = INT_MIN;

				for (int j = 0; j < 4; j++)
				{
					Vertex* v = qTemp->verts[j];
					// Current x is smaller than x1
					if (v->x < x1) {
						x1 = v->x;
					}

					// Current y is smaller than y1
					if (v->y < y1) {
						y1 = v->y;
					}

					// Current x is bigger than x2
					if (v->x > x2) {
						x2 = v->x;
					}

					// Current y is bigger than y2 
					if (v->y > y2) {
						y2 = v->y;
					}
				}

				for (int j = 0; j < 4; j++) {
					Vertex* v = qTemp->verts[j];
					if (v->x == x1 && v->y == y1) {
						fx1y1 = v->vx;
						gx1y1 = v->vy;
					}
					if (v->x == x2 && v->y == y2) {
						fx2y2 = v->vx;
						gx2y2 = v->vy;
					}
					if (v->x == x1 && v->y == y2) {
						fx1y2 = v->vx;
						gx1y2 = v->vy;
					}
					if (v->x == x2 && v->y == y1) {
						fx2y1 = v->vx;
						gx2y1 = v->vy;
					}
				}
			}

			double x0 = qTemp->singularity->x;
			double y0 = qTemp->singularity->y;

			double dfdx = (-(y2 - y0) * fx1y1 + (y2 - y0) * fx2y1 - (y0 - y1) * fx1y2 + (y0 - y1) * fx2y2) / (x2 - x1) * (y2 - y1);
			double dfdy = (-(x2 - x0) * fx1y1 - (x0 - x1) * fx2y1 + (x2 - x0) * fx1y2 + (x0 - x1) * fx2y2) / (x2 - x1) * (y2 - y1);
			double dgdx = (-(y2 - y0) * gx1y1 + (y2 - y0) * gx2y1 - (y0 - y1) * gx1y2 + (y0 - y1) * gx2y2) / (x2 - x1) * (y2 - y1);
			double dgdy = (-(x2 - x0) * gx1y1 - (x0 - x1) * gx2y1 + (x2 - x0) * gx1y2 + (x0 - x1) * gx2y2) / (x2 - x1) * (y2 - y1);

			icMatrix2x2 matrix = icMatrix2x2(dfdx, dfdy, dgdx, dgdy);

			double det = determinant(matrix);

			if (det > 0) {
				std::cout << x0 << " " << y0 << "\t" << det << std::endl;
				sources.push_back(icVector3(x0, y0, 0));
			}
			else if (det == 0) {
				std::cout << x0 << " " << y0 << "\t" << det << std::endl;
				higher_order.push_back(icVector3(x0, y0, 0));
			}
			else if (det < 0) {
				std::cout << x0 << " " << y0 << "\t" << det << std::endl;
				saddles.push_back(icVector3(x0, y0, 0));
			}
		}

	}
}

Quad* streamline_step(icVector2& cpos, icVector2& npos, Quad* cquad, bool forward)
{
	double x1, y1, x2, y2, f11, f12, f21, f22, g11, g21, g12, g22;
	Vertex* v11, * v12, * v21, * v22;
	// 1. Find x1, x2, y1, y2, f11, f21, f12, f22, g11, g21, g12, g22, v11, v21, v12, v22
	// x1 is the smallest x coordinate of the 4 vertices
	// x2 is the largest x coordinate of the 4 vertices
	// y1 is the smallest y coordinate of the 4 vertices
	// y2 is the largest y coordinate of the 4 vertices
	// v11 is the vertex with coordinates (x1,y1)
	// f11 is the vector x component at vertex v11
	// g11 is the vector y component at vertex v11
	// ...
	// (continues on next slide) ...

	x1 = INT_MAX;
	x2 = INT_MIN;
	y1 = INT_MAX;
	y2 = INT_MIN;

	for (int j = 0; j < 4; j++)
	{
		Vertex* v = cquad->verts[j];
		// Current x is smaller than x1
		if (v->x < x1) {
			x1 = v->x;
		}

		// Current y is smaller than y1
		if (v->y < y1) {
			y1 = v->y;
		}

		// Current x is bigger than x2
		if (v->x > x2) {
			x2 = v->x;
		}

		// Current y is bigger than y2 
		if (v->y > y2) {
			y2 = v->y;
		}
	}

	for (int j = 0; j < 4; j++)
	{
		Vertex* v = cquad->verts[j];

		// Found v11
		if (v->x == x1 && v->y == y1) {
			v11 = v;
			f11 = v->vx;
			g11 = v->vy;
		}
		// Found v21
		if (v->x == x2 && v->y == y1) {
			v21 = v;
			f21 = v->vx;
			g21 = v->vy;
		}
		// Found v12
		if (v->x == x1 && v->y == y2) {
			v12 = v;
			f12 = v->vx;
			g12 = v->vy;
		}
		// Found v22
		if (v->x == x2 && v->y == y2) {
			v22 = v;
			f22 = v->vx;
			g22 = v->vy;
		}
	}




	double	x0 = cpos.x;
	double	y0 = cpos.y;

	icVector2 vect = icVector2(x0, y0);
	vect.x = (((x2 - x0) * (y2 - y0) * f11) + ((x0 - x1) * (y2 - y0) * f21) +
		((x2 - x0) * (y0 - y1) * f12) + ((x0 - x1) * (y0 - y1) * f22)) /
		((x2 - x1) * (y2 - y1));
	vect.y = (((x2 - x0) * (y2 - y0) * g11) + ((x0 - x1) * (y2 - y0) * g21) +
		((x2 - x0) * (y0 - y1) * g12) + ((x0 - x1) * (y0 - y1) * g22)) /
		((x2 - x1) * (y2 - y1));

	normalize(vect);

	if (!forward)  vect *= -1.0;

	icVector2 divVect = icVector2(vect);
	divVect.x *= STEP;
	divVect.y *= STEP;
	npos = cpos + divVect;

	Quad* nquad = cquad;
	if (npos.x < x1 || npos.x > x2 || npos.y < y1 || npos.y > y2) {
		icVector2 cross_x1, cross_y1, cross_x2, cross_y2;
		double dprod_x1, dprod_y1, dprod_x2, dprod_y2;
		Edge* cross_edge;

		/*
		Use parametric forms of lines to find the intersection along each edge

		parametric form -> l:x = x0 + t*f(x0,y0) and y = y0 + t*g(x0,y0)

		Once you have each crossing point, check whether the vector at cpos is point toward each crossing point
		*/

		double x0 = cpos.x;
		double y0 = cpos.y;
		double xStep = vect.x;
		double yStep = vect.y;

		// Crossing the x1 line
		double x1T = (x1 - x0) / xStep;
		cross_x1 = icVector2(x0 + x1T * xStep, y0 + x1T * yStep);

		// Crossing the x2 line
		double x2T = (x2 - x0) / xStep;
		cross_x2 = icVector2(x0 + x2T * xStep, y0 + x2T * yStep);

		// Crossing the y1 line
		double y1T = (y1 - y0) / yStep;
		cross_y1 = icVector2(x0 + y1T * xStep, y0 + y1T * yStep);

		// Crossing the y2 line
		double y2T = (y2 - y0) / yStep;
		cross_y2 = icVector2(x0 + y2T * xStep, y0 + y2T * yStep);


		dprod_x1 = dot(vect, cross_x1 - cpos);
		dprod_x2 = dot(vect, cross_x2 - cpos);
		dprod_y1 = dot(vect, cross_y1 - cpos);
		dprod_y2 = dot(vect, cross_y2 - cpos);


		// Find the crossing point that is both on the border of cquad and
		// has a dprod value greater than 0.
		if (cross_x1.y >= y1 && cross_x1.y <= y2 && dprod_x1 > 0)
		{
			npos = cross_x1;
			cross_edge = poly->find_edge(v11, v12);
			nquad = poly->other_quad(cross_edge, nquad);
		}
		else if (cross_x2.y >= y1 && cross_x2.y <= y2  && dprod_x2 > 0) {
			npos = cross_x2;
			cross_edge = poly->find_edge(v21, v22);
			nquad = poly->other_quad(cross_edge, cquad);
		}
		else if (cross_y1.x >= x1 && cross_y1.x <= x2 && dprod_y1 > 0) {
			npos = cross_y1;
			cross_edge = poly->find_edge(v11, v21);
			nquad = poly->other_quad(cross_edge, cquad);
		}
		else if (cross_y2.x >= x1 && cross_y2.x <= x2  && dprod_y2 > 0) {
			npos = cross_y2;
			cross_edge = poly->find_edge(v12, v22);
			nquad = poly->other_quad(cross_edge, cquad);
		}
		else {
			nquad = poly->find_quad(npos.x, npos.y);
		}
		// if none of the crossing points meet these conditions, use the poly->find_quad()
		// function with npos to get the appropriate nquad
	}

	double approxSing = sing_prox(npos);

	if (nquad == NULL)
		std::cout << "nquad is null " << std::endl;

	if (approxSing < STEP) return NULL;
	else return nquad;
}

PolyLine build_streamline(double x, double y)
{
	// 1. Initialize current quad, current position, new position, step counter, and PolyLine variables
	Quad* cquad = poly->find_quad(x, y); /*use the find_quad() method of the polyhedron class*/
	PolyLine pline;

	//	if (cquad->singularity == NULL)
	//		return pline;

	icVector2 cpos = icVector2(x, y);
	icVector2 npos;
	int step_counter = 0;


	// 2. Trace the streamline forward until you hit a singularity, the edge of the domain, or the max step limit
	while (cquad != NULL && step_counter < STEP_MAX)
	{
		cquad = streamline_step(cpos, npos, cquad, 1);

		 // std::cout << cpos.x << ", " << cpos.y << std::endl;
		// std::cout << npos.x << ", " << npos.y << std::endl;

		LineSegment line(cpos.x, cpos.y, 0, npos.x, npos.y, 1);
		pline.push_back(line);
		cpos = npos;
		step_counter += 1;
	}

	cquad = poly->find_quad(x, y);
	cpos = icVector2(x, y);
	step_counter = 0;

	while (cquad != NULL && step_counter < STEP_MAX)
	{
		cquad = streamline_step(cpos, npos, cquad, 0);

		LineSegment line(cpos.x, cpos.y, 0, npos.x, npos.y, 0);
		pline.push_back(line);
		cpos = npos;
		step_counter += 1;
	}

	return pline;
}

/******************************************************************************
Process a keyboard action.  In particular, exit the program when an
"escape" is pressed in the window.
******************************************************************************/

void keyboard(unsigned char key, int x, int y) {
	int i;

	// clear out lines and points
	lines.clear();
	points.clear();

	double scale = 0.0025;

	switch (key) {
	case 27:	// set excape key to exit program
		poly->finalize();  // finalize_everything
		exit(0);
		break;

	case '1':	// solid color display with lighting
		display_mode = 1;
		glutPostRedisplay();
		break;

	case '2':	// wireframe display
		display_mode = 2;
		glutPostRedisplay();
		break;

	case '3':	// checkerboard display
	{
		display_mode = 3;

		double L = (poly->radius * 2) / 30;
		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			for (int j = 0; j < 4; j++) {

				Vertex* temp_v = temp_q->verts[j];

				temp_v->R = int(temp_v->x / L) % 2 == 0 ? 1 : 0;
				temp_v->G = int(temp_v->y / L) % 2 == 0 ? 1 : 0;
				temp_v->B = 0.0;
			}
		}
		glutPostRedisplay();
	}
	break;

	case '4':	// Drawing points and lines created by the dots_and_lines_example() function
		display_mode = 4;
		dots_and_lines_example(&points, &lines);
		glutPostRedisplay();
		break;

	case '5':	// IBFV vector field display
		display_mode = 5;
		glutPostRedisplay();
		break;

	case '6':	// add your own display mode
		display_mode = 6;
		{
			find_singularities();
			classify_singularities();
		}
		glutPostRedisplay();
		break;

	case '7':	// add your own display mode
		display_mode = 7;
		{
			streamlines.clear();
			find_singularities();
			classify_singularities();
			double x = 0;
			double y = 0;

			std::cout << "Please enter a x value: ";
			std::cin >> x;
			std::cout << "Please enter a y value: ";
			std::cin >> y;

			PolyLine l = build_streamline(x, y);

			std::cout << "Streamline len: " << l.size() << std::endl;
			streamlines.push_back(l);
		}
		glutPostRedisplay();
		break;

	case '8':	// add your own display mode
		display_mode = 8;
		{
			streamlines.clear();
			find_singularities();
			classify_singularities();

			for (int i = -255; i < 255; i += 3)
			{
				for (int j = -255; j < 255; j += 3)
				{
					PolyLine l = build_streamline(i, j);
					streamlines.push_back(l);
				}
			}

		}
		glutPostRedisplay();
		break;
	case '9':	// add your own display mode
		display_mode = 9;
		{

		}
		glutPostRedisplay();
		break;
	case '0':
		display_mode = 0;
		{
		}
		break;


	case 'r':	// reset rotation and transformation
		mat_ident(rotmat);
		translation[0] = 0;
		translation[1] = 0;
		zoom = 1.0;
		glutPostRedisplay();
		break;
	}
}

/******************************************************************************
Callback function for dragging mouse
******************************************************************************/

void motion(int x, int y) {
	float r[4];
	float s, t;

	s = (2.0 * x - win_width) / win_width;
	t = (2.0 * (win_height - y) - win_height) / win_height;

	if ((s == s_old) && (t == t_old))
		return;

	switch (mouse_mode) {
	case 2:

		Quaternion rvec;

		mat_to_quat(rotmat, rvec);
		trackball(r, s_old, t_old, s, t);
		add_quats(r, rvec, rvec);
		quat_to_mat(rvec, rotmat);

		s_old = s;
		t_old = t;

		display();
		break;

	case 1:

		translation[0] += (s - s_old);
		translation[1] += (t - t_old);

		s_old = s;
		t_old = t;

		display();
		break;
	}
}

/******************************************************************************
Callback function for mouse clicks
******************************************************************************/

void mouse(int button, int state, int x, int y) {

	int key = glutGetModifiers();

	if (button == GLUT_LEFT_BUTTON || button == GLUT_RIGHT_BUTTON) {

		if (state == GLUT_DOWN) {
			float xsize = (float)win_width;
			float ysize = (float)win_height;

			float s = (2.0 * x - win_width) / win_width;
			float t = (2.0 * (win_height - y) - win_height) / win_height;

			s_old = s;
			t_old = t;

			/*translate*/
			if (button == GLUT_LEFT_BUTTON)
			{
				mouse_mode = 1;
			}

			/*rotate*/
			if (button == GLUT_RIGHT_BUTTON)
			{
				mouse_mode = 2;
			}
		}
		else if (state == GLUT_UP) {

			if (button == GLUT_LEFT_BUTTON && key == GLUT_ACTIVE_SHIFT) {  // build up the selection feedback mode

				/*select face*/

				GLuint selectBuf[512];
				GLint hits;
				GLint viewport[4];

				glGetIntegerv(GL_VIEWPORT, viewport);

				glSelectBuffer(win_width, selectBuf);
				(void)glRenderMode(GL_SELECT);

				glInitNames();
				glPushName(0);

				glMatrixMode(GL_PROJECTION);
				glPushMatrix();
				glLoadIdentity();

				/*create 5x5 pixel picking region near cursor location */
				gluPickMatrix((GLdouble)x, (GLdouble)(viewport[3] - y), 1.0, 1.0, viewport);

				set_view(GL_SELECT);
				set_scene(GL_SELECT, poly);
				display_quads(GL_SELECT, poly);

				glMatrixMode(GL_PROJECTION);
				glPopMatrix();
				glFlush();

				glMatrixMode(GL_MODELVIEW);

				hits = glRenderMode(GL_RENDER);
				poly->selected_quad = processHits(hits, selectBuf);
				printf("Selected quad id = %d\n", poly->selected_quad);
				glutPostRedisplay();

				CHECK_GL_ERROR();

			}
			else if (button == GLUT_LEFT_BUTTON && key == GLUT_ACTIVE_CTRL)
			{
				/*select vertex*/

				GLuint selectBuf[512];
				GLint hits;
				GLint viewport[4];

				glGetIntegerv(GL_VIEWPORT, viewport);

				glSelectBuffer(win_width, selectBuf);
				(void)glRenderMode(GL_SELECT);

				glInitNames();
				glPushName(0);

				glMatrixMode(GL_PROJECTION);
				glPushMatrix();
				glLoadIdentity();

				/*  create 5x5 pixel picking region near cursor location */
				gluPickMatrix((GLdouble)x, (GLdouble)(viewport[3] - y), 1.0, 1.0, viewport);

				set_view(GL_SELECT);
				set_scene(GL_SELECT, poly);
				display_vertices(GL_SELECT, poly);

				glMatrixMode(GL_PROJECTION);
				glPopMatrix();
				glFlush();

				glMatrixMode(GL_MODELVIEW);

				hits = glRenderMode(GL_RENDER);
				poly->selected_vertex = processHits(hits, selectBuf);
				printf("Selected vert id = %d\n", poly->selected_vertex);
				glutPostRedisplay();

			}

			mouse_mode = -1;
		}
	}
}

/******************************************************************************
Callback function for mouse wheel scroll
******************************************************************************/

void mousewheel(int wheel, int direction, int x, int y) {
	if (direction == 1) {
		zoom *= zoomspeed;
		glutPostRedisplay();
	}
	else if (direction == -1) {
		zoom /= zoomspeed;
		glutPostRedisplay();
	}
}

/******************************************************************************
Callback function for window reshaping
******************************************************************************/

void reshape(int width, int height)
{
	win_width = width;
	win_height = height;

	aspectRatio = (float)width / (float)height;

	glViewport(0, 0, width, height);

	set_view(GL_RENDER);

	// reset IBFV pixels buffer
	free(pixels);
	initIBFV();
}

/******************************************************************************
Display IBFV vector field visualization (used for Project 3)
******************************************************************************/

void displayIBFV()
{



	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glDisable(GL_BLEND);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

	glEnable(GL_TEXTURE_2D);
	glShadeModel(GL_FLAT);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glClearColor(0.5, 0.5, 0.5, 1.0);  // background for rendering color coding and lighting
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);





	// draw the mesh using pixels and use vector field to advect texture coordinates
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, win_width, win_height, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels);

	double modelview_matrix[16], projection_matrix[16];
	int viewport[4];
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix);
	glGetDoublev(GL_PROJECTION_MATRIX, projection_matrix);
	glGetIntegerv(GL_VIEWPORT, viewport);

	for (int i = 0; i < poly->nquads; i++)
	{
		Quad* qtemp = poly->qlist[i];

		glBegin(GL_QUADS);
		for (int j = 0; j < 4; j++)
		{
			Vertex* vtemp = qtemp->verts[j];

			double tx, ty, dummy;
			gluProject((GLdouble)vtemp->x, (GLdouble)vtemp->y, (GLdouble)vtemp->z,
				modelview_matrix, projection_matrix, viewport, &tx, &ty, &dummy);

			tx = tx / win_width;
			ty = ty / win_height;

			icVector2 dp = icVector2(vtemp->vx, vtemp->vy);
			normalize(dp);
			dp *= dmax;

			double dx = -dp.x;
			double dy = -dp.y;

			float px = tx + dx;
			float py = ty + dy;

			glTexCoord2f(px, py);
			glVertex3d(vtemp->x, vtemp->y, vtemp->z);
		}
		glEnd();
	}



	glEnable(GL_BLEND);

	// blend in noise pattern
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glTranslatef(-1.0, -1.0, 0.0);
	glScalef(2.0, 2.0, 1.0);

	glCallList(1);

	glBegin(GL_QUAD_STRIP);

	glTexCoord2f(0.0, 0.0);  glVertex2f(0.0, 0.0);
	glTexCoord2f(0.0, tmax); glVertex2f(0.0, 1.0);
	glTexCoord2f(tmax, 0.0);  glVertex2f(1.0, 0.0);
	glTexCoord2f(tmax, tmax); glVertex2f(1.0, 1.0);
	glEnd();
	glDisable(GL_BLEND);

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glReadPixels(0, 0, win_width, win_height, GL_RGB, GL_UNSIGNED_BYTE, pixels);

	// draw the mesh using pixels without advecting texture coords
	glClearColor(1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, win_width, win_height, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels);
	for (int i = 0; i < poly->nquads; i++)
	{
		Quad* qtemp = poly->qlist[i];
		glBegin(GL_QUADS);
		for (int j = 0; j < 4; j++)
		{
			Vertex* vtemp = qtemp->verts[j];
			double tx, ty, dummy;
			gluProject((GLdouble)vtemp->x, (GLdouble)vtemp->y, (GLdouble)vtemp->z,
				modelview_matrix, projection_matrix, viewport, &tx, &ty, &dummy);
			tx = tx / win_width;
			ty = ty / win_height;
			glTexCoord2f(tx, ty);
			glVertex3d(vtemp->x, vtemp->y, vtemp->z);
		}
		glEnd();
	}
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);
	glShadeModel(GL_SMOOTH);


	for (int k = 0; k < saddles.size(); k++) {
		icVector3 p = saddles[k];
		drawDot(p.x, p.y, p.z, 0.15, 1, 0, 1);
	}

	for (int k = 0; k < sources.size(); k++) {
		icVector3 p = sources[k];
		drawDot(p.x, p.y, p.z, 0.15, 0, 1, 1);
	}

	for (int k = 0; k < higher_order.size(); k++) {
		icVector3 p = higher_order[k];
		drawDot(p.x, p.y, p.z, 0.15, 1, 1, 0);
	}
}

/******************************************************************************
Callback function for scene display
******************************************************************************/

void display(void)
{
	glClearColor(1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	set_view(GL_RENDER);
	set_scene(GL_RENDER, poly);

	/*display the mesh*/
	display_polyhedron(poly);

	/*display selected elements*/
	display_selected_vertex(poly);
	display_selected_quad(poly);


	glFlush();
	glutSwapBuffers();
	glFinish();

	CHECK_GL_ERROR();
}

/******************************************************************************
Diaplay the polygon with visualization results
******************************************************************************/

void display_polyhedron(Polyhedron* poly)
{
	unsigned int i, j;

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);

	switch (display_mode)
	{
	case 1:	// solid color display with lighting
	{
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_LIGHT1);

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		GLfloat mat_diffuse[4] = { 0.24, 0.4, 0.47, 0.0 };
		GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialf(GL_FRONT, GL_SHININESS, 50.0);

		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}
	}
	break;

	case 2:	// wireframe display
	{
		glDisable(GL_LIGHTING);
		glEnable(GL_LINE_SMOOTH);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glLineWidth(1.0);

		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];

			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_q->normal.entry[0], temp_q->normal.entry[1], temp_q->normal.entry[2]);
				glColor3f(0.0, 0.0, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}

		glDisable(GL_BLEND);
	}
	break;

	case 3:	// checkerboard pattern display
	{
		glDisable(GL_LIGHTING);
		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glColor3f(temp_v->R, temp_v->G, temp_v->B);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}
	}
	break;

	case 4: // points and lines drawing example
	{
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_LIGHT1);

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		GLfloat mat_diffuse[4] = { 0.24, 0.4, 0.47, 0.0 };
		GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialf(GL_FRONT, GL_SHININESS, 50.0);

		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}

		// draw lines
		for (int k = 0; k < lines.size(); k++)
		{
			drawPolyLine(lines[k], 1.0, 1.0, 0.0, 0.0);
		}

		// draw points
		for (int k = 0; k < points.size(); k++)
		{
			icVector3 point = points[k];
			drawDot(point.x, point.y, point.z);
		}
		break;
	}
	break;

	case 5:	// IBFV vector field display
	{
		displayIBFV();
		glutPostRedisplay();
	}
	break;

	case 6: // add your own display mode
	{

		displayIBFV();
		glutPostRedisplay();
		break;
	}
	case 7: // add your own display mode
	{


		displayIBFV();
		glutPostRedisplay();
		// draw lines
		for (int k = 0; k < streamlines.size(); k++)
		{
			drawPolyLine(streamlines[k], 1.0, 0, 0.4, 0.4);
		}

		break;
	}
	break;
	case 8: // add your own display mode
	{
		displayIBFV();
		glutPostRedisplay();
		// draw lines
		for (int k = 0; k < streamlines.size(); k++)
		{
			drawPolyLine(streamlines[k], 1.0, 0, 0.4, 0.4);
		}
		break;
	}
	break;
	default:
	{
		// don't draw anything
	}

	}
}
