/*
Functions for learnply

Eugene Zhang, 2005
*/

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>

#include "glut.h"
#include "ply.h"
#include "icVector.H"
#include "icMatrix.H"
#include "trackball.h"
#include "tmatrix.h"
#include "polyhedron.h"
#include "utility.h"

// Our imports
#include "List.h"

static PlyFile *in_ply;

unsigned char orientation; // 0=ccw, 1=cw

FILE *this_file;
const int win_width = 1024;
const int win_height = 1024;

double radius_factor = 0.9;

int display_mode = 1;
double error_threshold = 1.0e-13;
char reg_model_name[128];
FILE *f;
int ACSIZE = 1;		 // for antialiasing
int view_mode = 0; // 0 = othogonal, 1=perspective
float s_old, t_old;
float rotmat[4][4];
static Quaternion rvec;

int mouse_mode = -2;	 // -2=no action, -1 = down, 0 = zoom, 1 = rotate x, 2 = rotate y, 3 = tranlate x, 4 = translate y, 5 = cull near 6 = cull far
int mouse_button = -1; // -1=no button, 0=left, 1=middle, 2=right
int last_x, last_y;

struct jitter_struct
{
	double x;
	double y;
} jitter_para;

jitter_struct ji1[1] = {{0.0, 0.0}};
jitter_struct ji16[16] = {
		{0.125, 0.125},
		{0.375, 0.125},
		{0.625, 0.125},
		{0.875, 0.125},
		{0.125, 0.375},
		{0.375, 0.375},
		{0.625, 0.375},
		{0.875, 0.375},
		{0.125, 0.625},
		{0.375, 0.625},
		{0.625, 0.625},
		{0.875, 0.625},
		{0.125, 0.875},
		{0.375, 0.875},
		{0.625, 0.875},
		{0.875, 0.875},
};

Polyhedron *poly;

// forward function declarations
void init(void);
void keyboard(unsigned char key, int x, int y);
void motion(int x, int y);
void display(void);
void mouse(int button, int state, int x, int y);
void display_shape(GLenum mode, Polyhedron *poly);

/******************************************************************************
Main program.
******************************************************************************/
int main(int argc, char *argv[])
{
	FILE *this_file = fopen("../data/3D/bunny.ply", "r");
	poly = new Polyhedron(this_file);
	fclose(this_file);
	mat_ident(rotmat);

	poly->initialize(); // initialize everything

	poly->calc_bounding_sphere();
	poly->calc_face_normals_and_area();
	poly->average_normals();
	orientation = poly->orientation;

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowPosition(20, 20);
	glutInitWindowSize(win_width, win_height);
	glutCreateWindow("Geometric Modeling");
	init();
	glutKeyboardFunc(keyboard);
	glutDisplayFunc(display);
	glutMotionFunc(motion);
	glutMouseFunc(mouse);
	glutMainLoop();
	poly->finalize(); // finalize everything

	return 0; /* ANSI C requires main to return int. */
}

void init(void)
{
	/* select clearing color */

	glClearColor(0.0, 0.0, 0.0, 0.0); // background
	glShadeModel(GL_FLAT);
	glPolygonMode(GL_FRONT, GL_FILL);

	glDisable(GL_DITHER);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	// may need it
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glEnable(GL_NORMALIZE);
	if (orientation == 0)
		glFrontFace(GL_CW);
	else
		glFrontFace(GL_CCW);
}

void initMatrix3x3(icMatrix3x3 M)
{
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			M.set(i, j, 0.0);
}

void initMatrix(Matrix M)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			M[i][j] = 0.0;
		}
	}
}

double *normalToTriangle(double **verts)
{
	double x = 0, y = 0, z = 0;
	double *res = new double[3];

	for (int j = 0; j < 2; j++)
	{
		x = x + (verts[j][1] - verts[j + 1][1]) * (verts[j][2] + verts[j + 1][2]);
		y = y + (verts[j][2] - verts[j + 1][2]) * (verts[j][0] + verts[j + 1][0]);
		z = z + (verts[j][0] - verts[j + 1][0]) * (verts[j][1] + verts[j + 1][1]);
	}

	x = x + (verts[2][1] - verts[0][1]) * (verts[2][2] + verts[0][2]);
	y = y + (verts[2][2] - verts[0][2]) * (verts[2][0] + verts[0][0]);
	z = z + (verts[2][0] - verts[0][0]) * (verts[2][1] + verts[0][1]);

	res[0] = x;
	res[1] = y;
	res[2] = z;

	return res;
}

double dotProduct(double *a, double *b, int d)
{
	double res = 0;
	for (int i = 0; i < d; i++)
		res += a[i] * b[i];
	return res;
}

double **transpose(double **A, int m, int n)
{
	double **B = new double *[n];
	for (int i = 0; i < n; i++)
		B[i] = new double[m];

	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			B[j][i] = A[i][j];

	return B;
}

double det3x3(double **A)
{
	double res = 0;
	res = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
	return res;
}

void sumTriangle(double *newC, double *t, Triangle *t_i)
{
	double a[3], b[3], c[3];

	double **M = new double *[3];
	for (int i = 0; i < 3; i++)
	{
		M[i] = new double[3];
		for (int j = 0; j < 3; j++)
		{
			M[i][j] = 0.0;
		}
	}

	double v0_x = t_i->verts[0]->x;
	double v0_y = t_i->verts[0]->y;
	double v0_z = t_i->verts[0]->z;

	double v1_x = t_i->verts[1]->x;
	double v1_y = t_i->verts[1]->y;
	double v1_z = t_i->verts[1]->z;

	double v2_x = t_i->verts[2]->x;
	double v2_y = t_i->verts[2]->y;
	double v2_z = t_i->verts[2]->z;

	a[0] = v0_x;
	a[1] = v0_y;
	a[2] = v0_z;

	b[0] = v1_x;
	b[1] = v1_y;
	b[2] = v1_z;

	c[0] = v2_x;
	c[1] = v2_y;
	c[2] = v2_z;

	double **verts = new double *[3];
	for (int i = 0; i < 3; i++)
	{
		verts[i] = new double[3];
		for (int j = 0; j < 3; j++)
		{
			verts[i][j] = 0.0;
		}
	}

	verts[0][0] = a[0];
	verts[0][1] = a[1];
	verts[0][2] = a[2];
	verts[1][0] = b[0];
	verts[1][1] = b[1];
	verts[1][2] = b[2];
	verts[2][0] = c[0];
	verts[2][1] = c[1];
	verts[2][2] = c[2];

	double *normal = normalToTriangle(verts);

	double a1 = normal[0];
	double b1 = normal[1];
	double c1 = normal[2];

	for (int i = 0; i < 3; i++)
	{
		M[0][i] = a[i];
		M[1][i] = b[i];
		M[2][i] = c[i];
	}

	double *fp = new double[3];
	for (int i = 0; i < 3; i++)
	{
		fp[i] = verts[0][i];
	}

	double d = -1 * dotProduct(normal, fp, 3);
	newC[0] += a1;
	newC[1] += b1;
	newC[2] += c1;

	double **transposed = transpose(M, 3, 3);
	*t = *t - det3x3(transposed);
}

double norm(double *a, int d)
{
	double res = 0;
	for (int i = 0; i < d; i++)
		res += a[i] * a[i];
	return sqrt(res);
}

double cosTwoVecs(double* a, double* b, int n){
	return (dotProduct(a, b, n) / (norm(a, n) * norm(b, n)));
}

double* crossProduct(double* a, double* b, int d){
	double* res = new double[d];
	for (int i = 0; i < d; i++)
	{
		res[i] = a[(i + 1) % d] * b[(i + 2) % d] - a[(i + 2) % d] * b[(i + 1) % d];
	}
	return res;
}

double sinTwoVecs(double* a, double* b, int n){
	double* c = crossProduct(a, b, n);
	return norm(c, n) / (norm(a, n) * norm(b, n));
}

int checkConstraints(icMatrix3x3* A_c, icVector3* b_c, int n){
	double a1[3], a2[3], a3[3];
	double left, right, **matrix;

	switch (n)
	{
		case 1:
			if(A_c->entry[0][0] <= 0.000001 && A_c->entry[0][1] <= 0.000001 && A_c->entry[0][2] <= 0.000001){
				return 0;
			}
			return 1;
			break;
		case 2:
			a1[0] = A_c->entry[0][0];
			a1[1] = A_c->entry[0][1];
			a1[2] = A_c->entry[0][2];

			a2[0] = A_c->entry[1][0];
			a2[1] = A_c->entry[1][1];
			a2[2] = A_c->entry[1][2];

			double cos = cosTwoVecs(a1, a2, 3);
			return cos * cos < 0.97;
			break;
		case 3:
			a1[0] = A_c->entry[0][0];
			a1[1] = A_c->entry[0][1];
			a1[2] = A_c->entry[0][2];

			a2[0] = A_c->entry[1][0];
			a2[1] = A_c->entry[1][1];
			a2[2] = A_c->entry[1][2];

			a3[0] = A_c->entry[2][0];
			a3[1] = A_c->entry[2][1];
			a3[2] = A_c->entry[2][2];

			double* cp = crossProduct(a1, a2, 3);

			left = dotProduct(cp, a3, 3);
			left = left * left;

			right = norm(cp, 3) * norm(a3, 3) * sinTwoVecs(cp, a3, 3);
			right = right * right;

			return left > right;
			break;
	}
	return 0;
}

int addConstraintIfIndep(icMatrix3x3 *A_c, icVector3 *b_c, int n, double* newC,  double b){
	for (int i = 0; i < 3; i++)
	{
		A_c->entry[n][i] = newC[i];
	}

	b_c[n] = b;

	if(checkConstraints(A_c, b_c, n+1)){
		n++;
	} else {
		for (int i = 0; i < 3; i++)
		{
			A_c->entry[n][i] = 0;
		}
	}

	return n;
}

int isBoundaryEdge(Edge *e){
	
	int v1_i = e->verts[0]->index;
	int v2_i = e->verts[1]->index;

	Vertex *v1 = poly->vlist[v1_i];
	Vertex *v2 = poly->vlist[v2_i];

	v1->
}

int getAllConstraints(icMatrix3x3 *A_c, icVector3 *b_c, Matrix *H1, Matrix *H2, int n_constraints, int edge)
{
	icVector3 *a1, *b1;
	double *newC = new double[3];
	double *t = (double *)malloc(sizeof(double));

	/* Set newC[i] = 0 */
	for (int i = 0; i < 3; i++)
	{
		newC[i] = 0.0;
	}

	/* Set t to 0 */
	*t = 0.0;

	/* Get out two vertexes of our edge */
	Vertex *v1 = poly->elist[edge]->verts[0];
	Vertex *v2 = poly->elist[edge]->verts[1];

	double a[3], b[3], c[3];

	double *e1, *e2, *e3;
	e1 = (double *) malloc(sizeof(double) * 3);
	e2 = (double *) malloc(sizeof(double) * 3);
	e3 = (double *) malloc(sizeof(double) * 3);

	/* Create two linked lists */
	LinkedList *list1 = new LinkedList();
	LinkedList *list2 = new LinkedList();

	createList(list1);
	createList(list2);

	switch (n_constraints)
	{
	case 0: // Constraint volume
		// Set the first element of a linked list, l1, to be v1
		addToList(list1, v1);
		// Set the first element of a linked list, l2, to be v2
		addToList(list2, v2);

		// Set the current node to be l1
		Node *current = list1->head;

		/* 
				While current is not null
					- Get the value of some triangle ?
					- Set current to current->next
					- Sum triangle components
			 */
		while (current != NULL)
		{
			Triangle *t_i = (Triangle *)current->value;
			current = current->next;
			sumTriangle(newC, t, t_i);
		}

		// Set the current node to be l2
		current = list2->head;

		/* 
				While current is not null
					If l1 does not contain current
						- Get the value of some triangle ?
						- Sum triangle components
					- Set current to current->next
			 */
		while (current != NULL) {
			if(containsVal(list1, current->value) == 0) {
				Triangle *t_i = (Triangle *)current->value;
				current = current->next;
				sumTriangle(newC, t, t_i);
			}
			else {
				current = current->next;
			}
		}

		// For each element in newC, set it to be 1/6 of its value
		for (int i = 0; i < 3; i++)
		{
			newC[i] = newC[i] / 6.0;
		}

		// Set t to be equal to 1/6 of its current value
		*t = *t / 6.0;

		// Add constraint if independent

		n_constraints = addConstraintIfIndep(A_c, b_c, n_constraints, newC, *t);

		// Set t to 0
		*t = 0.0;

	case 1: // Boundary preservation
			// Initialize e1[i] = e2[i] = e3[i] = 0
		for (int i = 0; i < 3; i++)
		{
			e1[i] = 0.0;
			e2[i] = 0.0;
			e3[i] = 0.0;
		}

		// l1 = v1->edge
		destroyList(list1);
		createList(list1);

		int v1_i = v1->index;
		Edge *e1_i = poly->elist[v1_i];
		addToList(list1, e1_i);

		// l2 = v2->edge
		destroyList(list2);
		createList(list2);
		int v2_i = v2->index;
		Edge *e2_i = poly->elist[v2_i];
		addToList(list2, e2_i);

		/* 
				While l2 is not NULL
					- If isBoundaryEdge && l2->value != current edge
						- counter ++
						- Set current edge to l2->value
						- Set v1_e1 to v1
						- Set v2_e1 to v2
						- For i to 3
							- Set e1[i] to e1[i] + v2_e1[i] - v1_e1[i]
						- Find cross product between v2_e1 and v1_e1
						- For i to 3
							- Set e2[i] to e2[i] + cross_product[i]
					- Set l2 to l2->next
			 */
		while (list2->head != NULL) {

		}
		/* 
				While l1 is not NULL
					- If isBoundaryEdge && l1->value != current edge
						- counter ++
						- Set current edge to l2->value
						- Set v1_e1 to v1
						- Set v2_e1 to v2
						- For i to 3
							- Set e1[i] to e1[i] + v2_e1[i] - v1_e1[i]
						- Find cross product between v2_e1 and v1_e1
						- For i to 3
							- Set e2[i] to e2[i] + cross_product[i]
					- Set l1 to l1->next
			 */

		// If counter > 0
		// Set v1_e1 to v1
		// Set v2_e1 to v2
		// For i to 3
		// Set e1[i] to e1[i] + v2_e1[i] - v1_e1[i]
		// Find cross product between v2_e1 and v1_e1
		// For i to 3
		// Set e2[i] to e2[i] + cross_product[i]
		// Set e3 to corssProduct of e1 and e2
		// Set newC to be scalar vector of the dot product(e1, e1, 3), e3, 3
		// Get newly accepted constraints if independent
		// newC = crossProduct of e1 e3
		// n_constraints updated
	}
}

/******************************************************************************
Process a keyboard action.  In particular, exit the program when an
"escape" is pressed in the window.
******************************************************************************/
void keyboard(unsigned char key, int x, int y)
{
	int i;

	/* set escape key to exit */
	switch (key)
	{
	case 27:
		poly->finalize(); // finalize_everything
		exit(0);
		break;

	case '0':
		display_mode = 0;
		display();
		break;

	case '1':
		display_mode = 1;
		display();
		break;

	case '2':
		display_mode = 2;
		display();
		break;

	case '3':
		display_mode = 3;
		display();
		break;

	case '4':
		display_mode = 4;
		{

			icMatrix3x3 *A = new icMatrix3x3;
			icVector3 *b = new icVector3(0.0, 0.0, 0.0);

			Matrix *H1;
			Matrix *H2;
			icVector3 *c = new icVector3(0.0, 0.0, 0.0);
			double k = 0;

			initMatrix3x3(*A);
			initMatrix(*H1);
			initMatrix(*H2);

			double solution, vv1, vv2, cv;

			/* Get all contrains given some poly p, our constraints A and b, H1, H2, 0, n */
		}

		display();
		break;

	case '5':
		display_mode = 5;
		display();
		break;

	case '6':
		display_mode = 6;
		display();
		break;

	case '7':
		display_mode = 7;
		display();
		break;

	case '8':
		display_mode = 8;
		display();
		break;

	case '9':
		display_mode = 9;
		display();
		break;

	case 'x':
		switch (ACSIZE)
		{
		case 1:
			ACSIZE = 16;
			break;

		case 16:
			ACSIZE = 1;
			break;

		default:
			ACSIZE = 1;
			break;
		}
		fprintf(stderr, "ACSIZE=%d\n", ACSIZE);
		display();
		break;

	case '|':
		this_file = fopen("rotmat.txt", "w");
		for (i = 0; i < 4; i++)
			fprintf(this_file, "%f %f %f %f\n", rotmat[i][0], rotmat[i][1], rotmat[i][2], rotmat[i][3]);
		fclose(this_file);
		break;

	case '^':
		this_file = fopen("rotmat.txt", "r");
		for (i = 0; i < 4; i++)
			fscanf(this_file, "%f %f %f %f ", (&rotmat[i][0]), (&rotmat[i][1]), (&rotmat[i][2]), (&rotmat[i][3]));
		fclose(this_file);
		display();
		break;
	}
}

void multmatrix(const Matrix m)
{
	int i, j, index = 0;

	GLfloat mat[16];

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			mat[index++] = m[i][j];

	glMultMatrixf(mat);
}

void set_view(GLenum mode, Polyhedron *poly)
{
	icVector3 up, ray, view;
	GLfloat light_ambient0[] = {0.3, 0.3, 0.3, 1.0};
	GLfloat light_diffuse0[] = {0.7, 0.7, 0.7, 1.0};
	GLfloat light_specular0[] = {0.0, 0.0, 0.0, 1.0};
	GLfloat light_ambient1[] = {0.0, 0.0, 0.0, 1.0};
	GLfloat light_diffuse1[] = {0.5, 0.5, 0.5, 1.0};
	GLfloat light_specular1[] = {0.0, 0.0, 0.0, 1.0};
	GLfloat light_ambient2[] = {1.0, 1.0, 1.0, 1.0};
	GLfloat light_diffuse2[] = {1.0, 1.0, 1.0, 1.0};
	GLfloat light_specular2[] = {1.0, 1.0, 1.0, 1.0};
	GLfloat light_position[] = {0.0, 0.0, 0.0, 1.0};

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);
	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);

	glMatrixMode(GL_PROJECTION);
	if (mode == GL_RENDER)
		glLoadIdentity();

	if (view_mode == 0)
		glOrtho(-radius_factor, radius_factor, -radius_factor, radius_factor, 0.0, 40.0);
	else
		gluPerspective(45.0, 1.0, 0.1, 40.0);

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

void set_scene(GLenum mode, Polyhedron *poly)
{
	glTranslatef(0.0, 0.0, -3.0);
	multmatrix(rotmat);

	glScalef(1.0 / poly->radius, 1.0 / poly->radius, 1.0 / poly->radius);
	glTranslatef(-poly->center.entry[0], -poly->center.entry[1], -poly->center.entry[2]);
}

void motion(int x, int y)
{
	float r[4];
	float xsize, ysize, s, t;

	switch (mouse_mode)
	{
	case -1:

		xsize = (float)win_width;
		ysize = (float)win_height;

		s = (2.0 * x - win_width) / win_width;
		t = (2.0 * (win_height - y) - win_height) / win_height;

		if ((s == s_old) && (t == t_old))
			return;

		mat_to_quat(rotmat, rvec);
		trackball(r, s_old, t_old, s, t);
		add_quats(r, rvec, rvec);
		quat_to_mat(rvec, rotmat);

		s_old = s;
		t_old = t;

		display();
		break;
	}
}

int processHits(GLint hits, GLuint buffer[])
{
	unsigned int i, j;
	GLuint names, *ptr;
	double smallest_depth = 1.0e+20, current_depth;
	int seed_id = -1;
	unsigned char need_to_update;

	printf("hits = %d\n", hits);
	ptr = (GLuint *)buffer;
	for (i = 0; i < hits; i++)
	{ /* for each hit  */
		need_to_update = 0;
		names = *ptr;
		ptr++;

		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth)
		{
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth)
		{
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		for (j = 0; j < names; j++)
		{ /* for each name */
			if (need_to_update == 1)
				seed_id = *ptr - 1;
			ptr++;
		}
	}
	printf("triangle id = %d\n", seed_id);
	return seed_id;
}

void mouse(int button, int state, int x, int y)
{
	if (button == GLUT_LEFT_BUTTON || button == GLUT_RIGHT_BUTTON)
	{
		switch (mouse_mode)
		{
		case -2: // no action
			if (state == GLUT_DOWN)
			{
				float xsize = (float)win_width;
				float ysize = (float)win_height;

				float s = (2.0 * x - win_width) / win_width;
				float t = (2.0 * (win_height - y) - win_height) / win_height;

				s_old = s;
				t_old = t;

				mouse_mode = -1; // down
				mouse_button = button;
				last_x = x;
				last_y = y;
			}
			break;

		default:
			if (state == GLUT_UP)
			{
				button = -1;
				mouse_mode = -2;
			}
			break;
		}
	}
	else if (button == GLUT_MIDDLE_BUTTON)
	{
		if (state == GLUT_DOWN)
		{ // build up the selection feedback mode

			GLuint selectBuf[win_width];
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
			gluPickMatrix((GLdouble)x, (GLdouble)(viewport[3] - y),
										1.0, 1.0, viewport);

			set_view(GL_SELECT, poly);
			glPushMatrix();
			set_scene(GL_SELECT, poly);
			display_shape(GL_SELECT, poly);
			glPopMatrix();
			glFlush();

			hits = glRenderMode(GL_RENDER);
			poly->seed = processHits(hits, selectBuf);
			display();
		}
	}
}

void display_object()
{
	unsigned int i, j;
	Polyhedron *the_patch = poly;
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	for (i = 0; i < poly->ntris; i++)
	{
		Triangle *temp_t = poly->tlist[i];
		glBegin(GL_POLYGON);
		GLfloat mat_diffuse[] = {1.0, 1.0, 1.0, 1.0};

		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);

		glColor3f(1.0, 1.0, 1.0);
		glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
		for (j = 0; j < 3; j++)
		{
			Vertex *temp_v = temp_t->verts[j];
			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
		}
		glEnd();
	}
}

void display_shape(GLenum mode, Polyhedron *this_poly)
{
	unsigned int i, j;
	GLfloat mat_diffuse[4];

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);

	for (i = 0; i < this_poly->ntris; i++)
	{
		if (mode == GL_SELECT)
			glLoadName(i + 1);

		Triangle *temp_t = this_poly->tlist[i];

		switch (display_mode)
		{
		case 1:
			if (i == this_poly->seed)
			{
				mat_diffuse[0] = 0.0;
				mat_diffuse[1] = 0.0;
				mat_diffuse[2] = 1.0;
				mat_diffuse[3] = 1.0;
			}
			else
			{
				mat_diffuse[0] = 1.0;
				mat_diffuse[1] = 1.0;
				mat_diffuse[2] = 0.0;
				mat_diffuse[3] = 1.0;
			}
			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			glBegin(GL_POLYGON);
			for (j = 0; j < 3; j++)
			{

				Vertex *temp_v = temp_t->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				if (i == this_poly->seed)
					glColor3f(0.0, 0.0, 1.0);
				else
					glColor3f(1.0, 1.0, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;

		case 2:
			glBegin(GL_POLYGON);
			for (j = 0; j < 3; j++)
			{
				Vertex *temp_v = temp_t->verts[j];
				glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
				glColor3f(1.0, 1.0, 1.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;

		case 3:
			glBegin(GL_POLYGON);
			for (j = 0; j < 3; j++)
			{
				mat_diffuse[0] = 1.0;
				mat_diffuse[1] = 0.0;
				mat_diffuse[2] = 0.0;
				mat_diffuse[3] = 1.0;

				glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);

				Vertex *temp_v = temp_t->verts[j];
				glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);

				glColor3f(1.0, 0.0, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;

		case 4:
			glBegin(GL_POLYGON);
			for (j = 0; j < 3; j++)
			{
				Vertex *temp_v = temp_t->verts[j];
				glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
				glColor3f(1.0, 1.0, 1.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}
	}
}

void display(void)
{
	GLint viewport[4];
	int jitter;

	glClearColor(1.0, 1.0, 1.0, 1.0); // background for rendering color coding and lighting
	glGetIntegerv(GL_VIEWPORT, viewport);

	glClear(GL_ACCUM_BUFFER_BIT);
	for (jitter = 0; jitter < ACSIZE; jitter++)
	{
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		set_view(GL_RENDER, poly);
		glPushMatrix();
		switch (ACSIZE)
		{
		case 1:
			glTranslatef(ji1[jitter].x * 2.0 / viewport[2], ji1[jitter].y * 2.0 / viewport[3], 0.0);
			break;

		case 16:
			glTranslatef(ji16[jitter].x * 2.0 / viewport[2], ji16[jitter].y * 2.0 / viewport[3], 0.0);
			break;

		default:
			glTranslatef(ji1[jitter].x * 2.0 / viewport[2], ji1[jitter].y * 2.0 / viewport[3], 0.0);
			break;
		}
		set_scene(GL_RENDER, poly);
		display_shape(GL_RENDER, poly);
		glPopMatrix();
		glAccum(GL_ACCUM, 1.0 / ACSIZE);
	}
	glAccum(GL_RETURN, 1.0);
	glFlush();
	glutSwapBuffers();
	glFinish();
}