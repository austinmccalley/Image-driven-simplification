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
#include <vector>

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

int cases[4];
int edges_simplified = 0;
int trianglesDeleted = 0;

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

void initMatrix(icMatrix3x3 *M)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			M->entry[i][j] = 0.0;
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

icVector3 *normalToTriangle(icMatrix3x3 *M)
{
	double x = 0, y = 0, z = 0;

	icVector3 *res = new icVector3();

	for (int j = 0; j < 2; j++)
	{
		x = x + (M->entry[j][1] - M->entry[j + 1][1]) * (M->entry[j][2] + M->entry[j + 1][2]);
		y = y + (M->entry[j][2] - M->entry[j + 1][2]) * (M->entry[j][0] + M->entry[j + 1][0]);
		z = z + (M->entry[j][0] - M->entry[j + 1][0]) * (M->entry[j][1] + M->entry[j + 1][1]);
	}

	x = x + (M->entry[2][1] - M->entry[0][1]) * (M->entry[2][2] + M->entry[0][2]);
	y = y + (M->entry[2][2] - M->entry[0][2]) * (M->entry[2][0] + M->entry[0][0]);
	z = z + (M->entry[2][0] - M->entry[0][0]) * (M->entry[2][1] + M->entry[0][1]);

	res->x = x;
	res->y = y;
	res->z = z;

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

void sumTriangle(icVector3 *newC, double *t, Triangle *t_i)
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

	newC->entry[0] += a1;
	newC->entry[0] += b1;
	newC->entry[0] += c1;

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

double cosTwoVecs(double *a, double *b, int n)
{
	return (dotProduct(a, b, n) / (norm(a, n) * norm(b, n)));
}

double *crossProduct(double *a, double *b, int d)
{
	double *res = new double[d];
	for (int i = 0; i < d; i++)
	{
		res[i] = a[(i + 1) % d] * b[(i + 2) % d] - a[(i + 2) % d] * b[(i + 1) % d];
	}
	return res;
}

double sinTwoVecs(double *a, double *b, int n)
{
	double *c = crossProduct(a, b, n);
	return norm(c, n) / (norm(a, n) * norm(b, n));
}

int checkConstraints(icMatrix3x3 *A_c, icVector3 *b_c, int n)
{
	double a1[3], a2[3], a3[3];
	double left, right, **matrix;

	switch (n)
	{
	case 1:
	{
		if (A_c->entry[0][0] <= 0.000001 && A_c->entry[0][1] <= 0.000001 && A_c->entry[0][2] <= 0.000001)
		{
			return 0;
		}
		return 1;
		break;
	}
	case 2:
	{
		a1[0] = A_c->entry[0][0];
		a1[1] = A_c->entry[0][1];
		a1[2] = A_c->entry[0][2];

		a2[0] = A_c->entry[1][0];
		a2[1] = A_c->entry[1][1];
		a2[2] = A_c->entry[1][2];

		double cos = cosTwoVecs(a1, a2, 3);
		return cos * cos < 0.97;
		break;
	}
	case 3:
	{
		a1[0] = A_c->entry[0][0];
		a1[1] = A_c->entry[0][1];
		a1[2] = A_c->entry[0][2];

		a2[0] = A_c->entry[1][0];
		a2[1] = A_c->entry[1][1];
		a2[2] = A_c->entry[1][2];

		a3[0] = A_c->entry[2][0];
		a3[1] = A_c->entry[2][1];
		a3[2] = A_c->entry[2][2];

		double *cp = crossProduct(a1, a2, 3);

		left = dotProduct(cp, a3, 3);
		left = left * left;

		right = norm(cp, 3) * norm(a3, 3) * sinTwoVecs(cp, a3, 3);
		right = right * right;

		return left > right;
		break;
	}
	}
	return 0;
}

int addConstraintIfIndep(icMatrix3x3 *A_c, icVector3 *b_c, int n, icVector3 *newC, double b)
{
	for (int i = 0; i < 3; i++)
	{
		A_c->entry[n][i] = newC->entry[i];
	}

	b_c->entry[n] = b;

	if (checkConstraints(A_c, b_c, n + 1))
	{
		n++;
	}
	else
	{
		for (int i = 0; i < 3; i++)
		{
			A_c->entry[n][i] = 0;
		}
	}

	return n;
}

int isBoundaryEdge(Edge *e)
{
	return poly->getBoundingEdges(e) == 1;
}

double *orthogonalVtoV(icVector3 *a)
{
	icVector3 *res = new icVector3();

	if (fabs(a->entry[2]) >= 0.00001)
	{
		res->entry[0] = 1;
		res->entry[1] = 1;
		res->entry[2] = (-a->entry[0] - a->entry[1]) / a->entry[2];
	}
	else if (fabs(a->entry[1]) >= 0.00001)
	{
		res->entry[0] = 1;
		res->entry[1] = (-a->entry[0] - a->entry[2]) / a->entry[1];
		res->entry[2] = 1;
	}
	else if (fabs(a->entry[0]) >= 0.00001)
	{
		res->entry[0] = (-a->entry[1] - a->entry[2]) / a->entry[0];
		res->entry[1] = 1;
		res->entry[2] = 1;
	}
	else
	{
		res->entry[0] = 1;
		res->entry[1] = 1;
		res->entry[2] = 1;
	}

	res->entry[0] *= 10000;
	res->entry[1] *= 10000;
	res->entry[2] *= 10000;

	double *d_res = new double[3];
	d_res[0] = res->entry[0];
	d_res[1] = res->entry[1];
	d_res[2] = res->entry[2];

	return d_res;
}

double *orthogonalVto2V(icVector3 *a, double *b)
{
	icVector3 *b_vec = new icVector3(b[0], b[1], b[2]);

	icVector3 *res = cross(a, b_vec);

	for (int i = 0; i < 3; i++)
	{
		res->entry[i] *= 10000;
	}

	double *d_res = new double[3];
	d_res[0] = res->entry[0];
	d_res[1] = res->entry[1];
	d_res[2] = res->entry[2];

	return d_res;
}

icVector3 *matvet(double **A, int r, int c, icVector3 *v, int d)
{
	if (c != d)
	{
		std::cout << "Error: matvet: c and d must be the same size" << std::endl;
		return NULL;
	}

	icVector3 *res = new icVector3(0, 0, 0);

	for (int i = 0; i < r; i++)
	{
		res->entry[i] = 0;
		for (int j = 0; j < c; j++)
		{
			res->entry[i] += A[i][j];
			res->entry[i] *= v->entry[j];
		}
	}

	return res;
}

double *matvet(double **A, int r, int c, double *v, int d)
{
	if (c != d)
	{
		std::cout << "Error: matvet: c and d must be the same size" << std::endl;
		return NULL;
	}

	double *res = new double[d];

	for (int i = 0; i < r; i++)
	{
		res[i] = 0;
		for (int j = 0; j < c; j++)
		{
			res[i] += A[i][j];
			res[i] *= v[j];
		}
	}

	return res;
}

double **matmat(double **A, int r1, int c1, double **B, int r2, int c2)
{
	if (c1 != r2)
	{
		std::cout << "Error: matmat: c1 and r2 must be the same size" << std::endl;
		return NULL;
	}

	double **res = new double *[r1];
	for (int i = 0; i < r1; i++)
	{
		res[i] = new double[c2];
		for (int j = 0; j < c2; j++)
		{
			res[i][j] = 0;
			for (int k = 0; k < c1; k++)
			{
				res[i][j] += A[i][k] * B[k][j];
			}
		}
	}

	return res;
}

double **MtoD(icMatrix3x3 *A)
{
	double **res = new double *[3];
	for (int i = 0; i < 3; i++)
	{
		res[i] = new double[3];
		for (int j = 0; j < 3; j++)
		{
			res[i][j] = A->entry[i][j];
		}
	}
	return res;
}

int quadraticOpt(icMatrix3x3 *A_c, icVector3 *b_c, int n, double **A)
{
	double **Q = new double *[3 - n];
	for (int i = 0; i < 3 - n; i++)
	{
		Q[i] = new double[3];
	}

	if (n == 0)
	{
		for (int i = 0; i < 3; i++)
		{
			Q[i][i] = 100000.0;
		}
	}

	icVector3 *curr = new icVector3();

	curr->entry[0] = A_c->entry[0][0];
	curr->entry[1] = A_c->entry[0][1];
	curr->entry[2] = A_c->entry[0][2];

	double mult = 10000000.0;

	if (n == 1)
	{
		free(Q[0]);
		free(Q[1]);
		Q[0] = orthogonalVtoV(curr);
		Q[1] = orthogonalVto2V(curr, Q[0]);

		for (int i = 0; i < 3; i++)
			Q[0][i] *= mult;
		for (int i = 0; i < 3; i++)
			Q[1][i] *= mult;
	}

	if (n == 2)
	{
		icVector3 *curr2 = new icVector3();
		for (int i = 0; i < 3; i++)
		{
			curr->entry[i] = A[0][i];
			curr2->entry[i] = A[1][i];
		}
		free(Q[0]);

		double *tmp = new double[3];
		tmp[0] = curr2->entry[0];
		tmp[1] = curr2->entry[1];
		tmp[2] = curr2->entry[2];

		Q[0] = orthogonalVto2V(curr, tmp);

		for (int i = 0; i < 3; i++)
			Q[0][i] *= mult;

		free(curr2);
	}

	icMatrix3x3 *A_red = new icMatrix3x3();
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			A_red->entry[i][j] = A[i][j];
		}
	}

	icVector3 *b_tmp = new icVector3();
	b_tmp->entry[0] = A[0][3];
	b_tmp->entry[1] = A[1][3];
	b_tmp->entry[2] = A[2][3];

	icVector3 *Qb = matvet(Q, 3 - n, 3, b_tmp, 3);
	double **QA = matmat(Q, 3 - n, 3, MtoD(A_red), 3, 3);

	int n_iters = 3 - n;

	for (int i = 0; i < n_iters; i++)
	{
		icVector3 *QA_i = new icVector3();
		QA_i->entry[0] = QA[i][0];
		QA_i->entry[1] = QA[i][1];
		QA_i->entry[2] = QA[i][2];
		n = addConstraintIfIndep(A_c, b_c, n, QA_i, Qb->entry[i]);
	}

	/* TODO: Free memory */
	return n;
}

void calculateBoundaryEdge(double **H, icVector3 *e1, icVector3 *e2)
{
	/* Define a 4x4 matrix x_e1 */
	double **x_e1 = new double *[4];
	for (int i = 0; i < 4; i++)
	{
		x_e1[i] = new double[4];
		for (int j = 0; j < 4; j++)
		{
			x_e1[i][j] = 0;
		}
	}

	double **e1_x;

	x_e1[0][0] = -e1->entry[2];
	x_e1[1][0] = e1->entry[2];
	x_e1[0][2] = e1->entry[1];
	x_e1[2][0] = -e1->entry[1];
	x_e1[1][2] = -e1->entry[0];
	x_e1[2][1] = e1->entry[0];

	e1_x = transpose(x_e1, 3, 3);

	double **p = matmat(x_e1, 3, 3, e1_x, 3, 3);

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			H[i][j] += p[i][j];

	icVector3 *cp = cross(e1, e2);

	for (int i = 0; i < 3; i++)
	{
		H[i][3] += cp->entry[i];

		/* Free cp */

		H[3][3] = 1 / 2 * dot(e1, e2);

		/* Free x_e1, e1_x, p */
	}
	{
		/* code */
	}
}

int getAllConstraints(icMatrix3x3 *A_c, icVector3 *b_c, double **H1, double **H2, int n_constraints, int edge, icVector3 *c1, icVector3 *c2, double k1, double k2)
{
	icVector3 *a1, *b1, *newC;

	double *t = (double *)malloc(sizeof(double));

	/* Set newC[i] = 0 */
	newC = new icVector3(0, 0, 0);

	/* Set t to 0 */
	*t = 0.0;

	/* Get out two vertexes of our edge */
	Vertex *v1 = poly->elist[edge]->verts[0];
	Vertex *v2 = poly->elist[edge]->verts[1];

	double a[3], b[3], c[3];

	icVector3 *e1 = new icVector3();
	icVector3 *e2 = new icVector3();
	icVector3 *e3 = new icVector3();

	icVector3 *v1_e1 = new icVector3();
	icVector3 *v2_e1 = new icVector3();

	icVector3 *crossP = new icVector3();

	/* Counter */
	int counter = 0;

	/* Current edge */
	Edge *currE;

	/* Cross product */
	icVector3 *cp = new icVector3();

	switch (n_constraints)
	{
	case 0:
	{ // Constraint volume

		/*
				While current is not null
					- Get the value of some triangle ?
					- Set current to current->next
					- Sum triangle components
			 */

		std::vector<int> triangles;

		for (int i = 0; i < v1->ntris; i++)
		{
			Triangle *t_i = v1->tris[i];

			triangles.push_back(t_i->index);
			sumTriangle(newC, t, t_i);
		}

		/*
				While current is not null
					If l1 does not contain current
						- Get the value of some triangle ?
						- Sum triangle components
					- Set current to current->next
			 */

		for (int i = 0; i < v2->ntris; i++)
		{
			if (std::find(triangles.begin(), triangles.end(), v2->tris[i]->index) == triangles.end())
			{
				Triangle *t_i = v2->tris[i];

				triangles.push_back(t_i->index);
				sumTriangle(newC, t, t_i);
			}
		}

		// For each element in newC, set it to be 1/6 of its value
		for (int i = 0; i < 3; i++)
		{
			newC->entry[i] /= 6.0;
		}

		// Set t to be equal to 1/6 of its current value
		*t = *t / 6.0;

		// Add constraint if independent

		n_constraints = addConstraintIfIndep(A_c, b_c, n_constraints, newC, *t);

		// Set t to 0
		*t = 0.0;
	}
	case 1:
	{ // Boundary preservation
		// Initialize e1[i] = e2[i] = e3[i] = 0
		e1->entry[0] = 0.0;
		e1->entry[1] = 0.0;
		e1->entry[2] = 0.0;

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

		for (int i = 0; i < v2->ntris; i++)
		{
			Triangle *t_i = v2->tris[i];
			for (int j = 0; j < 3; j++)
			{
				Edge *e = t_i->edges[j];
				if (isBoundaryEdge(e) && e->index != edge)
				{
					counter++;
					currE = e;

					e1->entry[0] = currE->verts[0]->x;
					e1->entry[1] = currE->verts[0]->y;
					e1->entry[2] = currE->verts[0]->z;

					e2->entry[0] = currE->verts[1]->x;
					e2->entry[1] = currE->verts[1]->y;
					e2->entry[2] = currE->verts[1]->z;

					v1_e1->entry[0] = v1->x;
					v1_e1->entry[1] = v1->y;
					v1_e1->entry[2] = v1->z;

					v2_e1->entry[0] = v2->x;
					v2_e1->entry[1] = v2->y;
					v2_e1->entry[2] = v2->z;

					for (int k = 0; k < 3; k++)
					{
						e1->entry[k] += v2_e1->entry[k] - v1_e1->entry[k];
					}

					cp = cross(v2_e1, v1_e1);

					for (int k = 0; k < 3; k++)
					{
						e2->entry[k] += cp->entry[k];
					}
				}
			}
		}

		/*
				While l1 is not NULL
					- If isBoundaryEdge && l1->value != current edge
						- counter ++
						- Set current edge to l1->value
						- Set v1_e1 to v1
						- Set v2_e1 to v2
						- For i to 3
							- Set e1[i] to e1[i] + v2_e1[i] - v1_e1[i]
						- Find cross product between v2_e1 and v1_e1
						- For i to 3
							- Set e2[i] to e2[i] + cross_product[i]
					- Set l1 to l1->next
			 */
		for (int i = 0; i < v1->ntris; i++)
		{
			Triangle *t_i = v1->tris[i];
			for (int j = 0; j < 3; j++)
			{
				Edge *e = t_i->edges[j];
				if (isBoundaryEdge(e) && e->index != edge)
				{
					counter++;
					currE = e;

					v1_e1->entry[0] = v1->x;
					v1_e1->entry[1] = v1->y;
					v1_e1->entry[2] = v1->z;

					v2_e1->entry[0] = v2->x;
					v2_e1->entry[1] = v2->y;
					v2_e1->entry[2] = v2->z;

					e2->operator+=(*v2_e1);
					e2->operator-=(*v1_e1);
				}
			}
		}

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

		if (counter > 0)
		{
			v1_e1->x = v1->x;
			v1_e1->y = v1->y;
			v1_e1->z = v1->z;

			v2_e1->x = v2->x;
			v2_e1->y = v2->y;
			v2_e1->z = v2->z;

			e1->operator+=(*v2_e1);
			e1->operator-=(*v1_e1);

			cp = cross(v2_e1, v1_e1);

			e1->operator+=(*cp);

			e3 = cross(e1, e2);

			newC = scalarVector(dot(e1, e1), e3);

			n_constraints = addConstraintIfIndep(A_c, b_c, n_constraints, newC, *t);
		}
	}
	case 2:
	{
		// Volume Optimization

		// Make a matrix called verts thats 3x3
		icMatrix3x3 *verts = new icMatrix3x3(0, 0, 0, 0, 0, 0, 0, 0, 0);

		// Make a matrix H that is 4x4
		double **H = new double *[4];

		for (int i = 0; i < 4; i++)
		{
			H[i] = new double[4];
			for (int j = 0; j < 4; j++)
			{
				H[i][j] = 0;
			}
		}

		// Set k = 0;
		int k = 0;

		// Set a shadow counter to 0
		int shadowCounter = 0;

		// Set the doubles a,b,c,d
		double a, b, c, d;

		// Make a matrix called transposed
		icMatrix3x3 *transposed = new icMatrix3x3(0, 0, 0, 0, 0, 0, 0, 0, 0);

		// Make a vector called firstPoint of size 3
		icVector3 *firstPoint = new icVector3(0, 0, 0);

		// Initialize some triangle t_i
		Triangle *t_i;

		/* 
			While current is not null
				- Set t_i to the current triangle
				- Set current to current->next
				- Set verts to be the matrix of t_i->v1,v2,v3
				- Find normal to triangle
				- Set a to be the x of the normal
				- Set b to be the y of the normal
				- Set c to be the z of the normal
				- Tranpose is verts transposed
				- d is the negative determinant of verts
		 		- From i..3,
				 	- fp[i] = verts[0][i]
				- Set H to be a,b,c,d ^2
		 */

		std::vector<int> triangles;

		for (int i = 0; i < v1->ntris; i++)
		{
			t_i = v1->tris[i];

			triangles.push_back(t_i->index);

			icVector3 *v1_t_i = new icVector3(t_i->verts[0]->x, t_i->verts[0]->y, t_i->verts[0]->z);
			icVector3 *v2_t_i = new icVector3(t_i->verts[1]->x, t_i->verts[1]->y, t_i->verts[1]->z);
			icVector3 *v3_t_i = new icVector3(t_i->verts[2]->x, t_i->verts[2]->y, t_i->verts[2]->z);

			verts = new icMatrix3x3(v1_t_i, v2_t_i, v3_t_i);

			icVector3 *normal = normalToTriangle(verts);

			a = normal->x;
			b = normal->y;
			c = normal->z;

			transposed = &transpose(verts);

			d = -1 * determinant(*transposed);

			firstPoint->x = verts->entry[0][0];
			firstPoint->y = verts->entry[0][1];
			firstPoint->z = verts->entry[0][2];

			H[0][0] += a * a;
			H[0][1] += a * b;
			H[0][2] += a * c;
			H[0][3] += a * d;

			H[1][0] += b * a;
			H[1][1] += b * b;
			H[1][2] += b * c;
			H[1][3] += b * d;

			H[2][0] += c * a;
			H[2][1] += c * b;
			H[2][2] += c * c;
			H[2][3] += c * d;

			H[3][0] += d * a;
			H[3][1] += d * b;
			H[3][2] += d * c;
			H[3][3] += d * d;
		}

		for (int i = 0; i < v2->ntris; i++)
		{
			t_i = v2->tris[i];
			// If t_i is not in triangles
			if (std::find(triangles.begin(), triangles.end(), t_i->index) == triangles.end())
			{
				icVector3 *v1_t_i = new icVector3(t_i->verts[0]->x, t_i->verts[0]->y, t_i->verts[0]->z);
				icVector3 *v2_t_i = new icVector3(t_i->verts[1]->x, t_i->verts[1]->y, t_i->verts[1]->z);
				icVector3 *v3_t_i = new icVector3(t_i->verts[2]->x, t_i->verts[2]->y, t_i->verts[2]->z);

				verts = new icMatrix3x3(v1_t_i, v2_t_i, v3_t_i);

				icVector3 *normal = normalToTriangle(verts);

				a = normal->x;
				b = normal->y;
				c = normal->z;

				transposed = &transpose(verts);

				d = -1 * determinant(*transposed);

				firstPoint->x = verts->entry[0][0];
				firstPoint->y = verts->entry[0][1];
				firstPoint->z = verts->entry[0][2];

				H[0][0] += a * a;
				H[0][1] += a * b;
				H[0][2] += a * c;
				H[0][3] += a * d;

				H[1][0] += b * a;
				H[1][1] += b * b;
				H[1][2] += b * c;
				H[1][3] += b * d;

				H[2][0] += c * a;
				H[2][1] += c * b;
				H[2][2] += c * c;
				H[2][3] += c * d;

				H[3][0] += d * a;
				H[3][1] += d * b;
				H[3][2] += d * c;
				H[3][3] += d * d;
			}
		}

		double mult = 1;

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				H[i][j] = H1[i][j] = H[i][j] / 18;
			}
		}

		n_constraints = quadraticOpt(A_c, b_c, n_constraints, H);

		// cases[n]++

		if (n_constraints == 3)
			break;
	}
	case 3: // Boundary Optimization
	{
		e1->entry[0] = 0.0;
		e1->entry[1] = 0.0;
		e1->entry[2] = 0.0;

		/* Allocate room for a 4x4 double matrix called H */
		double **H = new double *[4];
		for (int i = 0; i < 4; i++)
		{
			H[i] = new double[4];
			for (int j = 0; j < 4; j++)
			{
				H[i][j] = 0.0;
			}
		}

		int isFirst = 1;

		/* Iterate over all the edges corresponding with  */
		for (int i = 0; i < v2->ntris; i++)
		{
			Triangle *t_i = v2->tris[i];
			for (int j = 0; j < 3; j++)
			{
				Edge *e = t_i->edges[j];
				if (isBoundaryEdge(e) && e->index != edge)
				{
					counter++;
					currE = e;

					v1_e1->entry[0] = v1->x;
					v1_e1->entry[1] = v1->y;
					v1_e1->entry[2] = v1->z;

					v2_e1->entry[0] = v2->x;
					v2_e1->entry[1] = v2->y;
					v2_e1->entry[2] = v2->z;

					e2->operator+=(*v2_e1);
					e2->operator-=(*v1_e1);

					crossP = cross(v1_e1, v2_e1);

					e2->entry[0] = crossP->entry[0];
					e2->entry[1] = crossP->entry[1];
					e2->entry[2] = crossP->entry[2];

					calculateBoundaryEdge(H, e1, e2);
				}
			}
		}

		/* Iterate over all the edges corresponding with  */
		for (int i = 0; i < v1->ntris; i++)
		{
			Triangle *t_i = v1->tris[i];
			for (int j = 0; j < 3; j++)
			{
				Edge *e = t_i->edges[j];
				if (isBoundaryEdge(e) && e->index != edge)
				{
					counter++;
					currE = e;

					v1_e1->entry[0] = v1->x;
					v1_e1->entry[1] = v1->y;
					v1_e1->entry[2] = v1->z;

					v2_e1->entry[0] = v2->x;
					v2_e1->entry[1] = v2->y;
					v2_e1->entry[2] = v2->z;

					e2->operator+=(*v2_e1);
					e2->operator-=(*v1_e1);

					crossP = cross(v1_e1, v2_e1);

					e2->entry[0] = crossP->entry[0];
					e2->entry[1] = crossP->entry[1];
					e2->entry[2] = crossP->entry[2];

					calculateBoundaryEdge(H, e1, e2);
				}
			}
		}

		if (counter > 0)
		{
			v1_e1->entry[0] = v1->x;
			v1_e1->entry[1] = v1->y;
			v1_e1->entry[2] = v1->z;

			v2_e1->entry[0] = v2->x;
			v2_e1->entry[1] = v2->y;
			v2_e1->entry[2] = v2->z;

			e2->operator+=(*v2_e1);
			e2->operator-=(*v1_e1);

			crossP = cross(v1_e1, v2_e1);
			e2->entry[0] = crossP->entry[0];
			e2->entry[1] = crossP->entry[1];
			e2->entry[2] = crossP->entry[2];

			calculateBoundaryEdge(H, e1, e2);
			for (int i = 0; i < 4; i++)
				for (int j = 0; j < 4; j++)
					H[i][j] = H2[i][j] = H[i][j] / 2;
			n_constraints = quadraticOpt(A_c, b_c, n_constraints, H);
			cases[n_constraints]++;
		}

		if (n_constraints == 3)
			break;
	}
	case 4: // Triangle Shape
	{
		Vertex *a1 = new Vertex(0, 0, 0);
		a1->x = v1->x;
		a1->y = v1->y;
		a1->z = v1->z;

		Vertex *b1 = new Vertex(0, 0, 0);
		b1->x = v2->x;
		b1->y = v2->y;
		b1->z = v2->z;

		Edge *currE;
		Vertex *cv;
		int firstIter = 1;

		icVector3 *newC = new icVector3();
		for (int i = 0; i < 3; i++)
		{
			newC->entry[i] = 0.0;
			c1->entry[i] = 0.0;
			/* TODO: Define c1 */
		}

		/* Make H a 4x4 matrix of doubles */
		double **H = new double *[4];
		for (int i = 0; i < 4; i++)
		{
			H[i] = new double[4];
			for (int j = 0; j < 4; j++)
				H[i][j] = 0.0;
		}

		double k = 0;

		/* Look at the current edge of a1 */

		counter = 0;
		for (int i = 0; i < a1->ntris; i++)
		{
			Triangle *t_i = a1->tris[i];
			for (int j = 0; j < 3; j++)
			{
				currE = t_i->edges[j];

				if (edge != currE->index)
				{
					if (currE->verts[0]->index == v1->index || currE->verts[1]->index == v2->index)
					{
						int newCvIdx = v2->index;
						cv = poly->vlist[newCvIdx];
					}
					else
					{
						int newCvIdx = v1->index;
						cv = poly->vlist[newCvIdx];
					}
					counter++;
					newC->entry[0] -= cv->x;
					newC->entry[1] -= cv->y;
					newC->entry[2] -= cv->z;
					k += dot(newC, newC);
				}
			}
		}
		firstIter = 0;

		for (int i = 0; i < b1->ntris; i++)
		{
			Triangle *t_i = b1->tris[i];
			for (int j = 0; j < 3; j++)
			{
				currE = t_i->edges[j];

				if (edge != currE->index)
				{
					if (currE->verts[0]->index == v1->index || currE->verts[1]->index == v2->index)
					{
						int newCvIdx = v2->index;
						cv = poly->vlist[newCvIdx];
					}
					else
					{
						int newCvIdx = v1->index;
						cv = poly->vlist[newCvIdx];
					}
					counter++;
					newC->entry[0] -= cv->x;
					newC->entry[1] -= cv->y;
					newC->entry[2] -= cv->z;
					k += dot(newC, newC);
				}
			}
		}

		for (int i = 0; i < 3; i++)
		{
			H[i][i] = 2 * counter;
			H[i][3] = 2 * newC->entry[i];
		}

		int prev = n_constraints;
		n_constraints = quadraticOpt(A_c, b_c, n_constraints, H);
		cases[n_constraints]++;
		/* TODO: Free H */
		break;
	}
	}

	/* TODO: Free e1, e2, e3, newC, t, v1_e1, v2_e2 */

	return n_constraints;
}

icVector3 *solveLinearSystem(icMatrix3x3 *A, icVector3 *b)
{
	icVector3 *res = new icVector3();

	for (int i = 0; i < 3; i++)
	{
		res->entry[i] = 0;
	}

	double detA = determinant(*A);

	if (fabs(detA) < 0.0001)
	{
		return NULL;
	}

	double currDet = 0;
	icVector3 col;

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			col.entry[j] = A->entry[j][i];
			A->entry[j][i] = b->entry[j];
		}
		currDet = determinant(*A);
		res->entry[i] = currDet / detA;
		for (int j = 0; j < 3; j++)
		{
			A->entry[j][i] = col.entry[j];
		}
	}

	return res;
}

void calculateSolutions(int n)
{
	icMatrix3x3 *A_c = new icMatrix3x3(0.0);
	icVector3 *b_c = new icVector3(0.0, 0.0, 0.0);

	double **fv1_H = new double *[4];
	for (int i = 0; i < 4; i++)
	{
		fv1_H[i] = new double[4];
		for (int j = 0; j < 4; j++)
		{
			fv1_H[i][j] = 0;
		}
	}

	double **fv2_H = new double *[4];
	for (int i = 0; i < 4; i++)
	{
		fv2_H[i] = new double[4];
		for (int j = 0; j < 4; j++)
		{
			fv2_H[i][j] = 0;
		}
	}

	double *vv1;
	double vv2;

	icVector3 *c1 = new icVector3(0.0, 0.0, 0.0);
	icVector3 *c2 = new icVector3(0.0, 0.0, 0.0);

	int k1 = 0;
	int k2 = 0;

	double cv;

	icMatrix3x3 *currA = new icMatrix3x3(0.0);
	icVector3 *currB = new icVector3(0.0, 0.0, 0.0);

	int n_constr = getAllConstraints(A_c, b_c, fv1_H, fv2_H, 0, n, c1, c2, k1, k2);
	if (n_constr != 3)
	{
		poly->solutions[n]->entry[0] = poly->solutions[n]->entry[1] = poly->solutions[n]->entry[2] = 9999;
		poly->elist[n]->cost = 99999999999.99;
	}

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			currA->entry[i][j] = A_c->entry[i][j];
		}
	}

	for (int i = 0; i < 3; i++)
	{
		currB->entry[i] = b_c->entry[i];
	}

	icVector3 *sol = solveLinearSystem(currA, currB);

	if (sol == NULL)
	{
		poly->solutions[n]->entry[0] = poly->solutions[n]->entry[1] = poly->solutions[n]->entry[2] = 9999;
		poly->elist[n]->cost = 99999999999.99;
	}
	else
	{
		poly->solutions[n]->entry[0] = sol->entry[0];
		poly->solutions[n]->entry[1] = sol->entry[1];
		poly->solutions[n]->entry[2] = sol->entry[2];

		double *currSol = new double[4];
		for (int i = 0; i < 3; i++)
		{
			currSol[i] = sol->entry[i];
		}

		currSol[3] = 1;

		vv1 = matvet(fv1_H, 4, 4, currSol, 4);
		vv2 = dotProduct(vv1, currSol, 4);
		icVector3 *pointV1 = new icVector3(0.0, 0.0, 0.0);
		icVector3 *pointV2 = new icVector3(0.0, 0.0, 0.0);

		pointV1->entry[0] = poly->elist[n]->verts[0]->x;
		pointV1->entry[1] = poly->elist[n]->verts[0]->y;
		pointV1->entry[2] = poly->elist[n]->verts[0]->z;

		pointV2->entry[0] = poly->elist[n]->verts[1]->x;
		pointV2->entry[1] = poly->elist[n]->verts[1]->y;
		pointV2->entry[2] = poly->elist[n]->verts[1]->z;

		double len = distance(*pointV1, *pointV2);

		poly->elist[n]->cost += 1 / 2 * len * len * vv2;

		/* 
			TODO: Free vv1, currSol, pointV1, pointV2, sol
		*/
	}

	/* 
	TODO:
		Free currA, currB
		H1
		H2
		c1
		c2
		c
		k1
		k2
	 */

	/* Get all contrains given some poly p, our constraints A and b, H1, H2, 0, n */
}

void swapEdges(int i, int j)
{
	Edge *temp;
	icVector3 *tmpv;
	int tmpPos;

	temp = poly->elist[i];
	poly->elist[i] = poly->elist[j];
	poly->elist[j] = temp;

	tmpPos = poly->elist[i]->index;
	poly->elist[i]->index = poly->elist[j]->index;
	poly->elist[j]->index = tmpPos;

	tmpv = poly->solutions[i];
	poly->solutions[i] = poly->solutions[j];
	poly->solutions[j] = tmpv;
}

void quicksort(int l, int r)
{
	if (l >= r)
		return;

	int i = l;
	int j = r;

	double pivot = poly->elist[i]->cost;

	for (;;)
	{
		while (poly->elist[i]->cost < pivot)
			i++;
		while (poly->elist[j]->cost > pivot)
			j--;
		if (i >= j)
			break;
		// Build out swap edges
		swapEdges(i, j);
		i++;
		j--;
	}

	quicksort(l, i - 1);
	quicksort(j + 1, r);
}

Vertex *getThirdVertex(Triangle *t, Vertex *v1, Vertex *v2)
{
	if (t->verts[0] == v1 || t->verts[0] == v2)
	{
		if (t->verts[1] == v1 || t->verts[1] == v2)
		{
			return t->verts[3];
		}
		else
		{
			return t->verts[1]
		}
	}
	else
	{
		return t->verts[0];
	}
}

Edge *getLinkingEdge(Vertex *v1, Vertex *v2)
{
	Edge *edge;

	for (int i = 0; i < v1->ntris; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			edge = v1->tris[i]->edges[j];

			for (int k = 0; k < v2->ntris; k++)
			{
				for (int l = 0; l < 3; l++)
				{
					if (edge == v2->tris[k]->edges[l])
						return edge;
				}
			}
		}
	}

	return NULL;
}

int leftArray(int p)
{
	if ((2 * p) > poly->nedges)
		return -1;
	else
		return 2 * (p - 1);
}

int rightArray(int p)
{
	if ((2 * p + 1) > poly->nedges)
		return -1;
	else
		return (2 * p);
}

void heapify(int i)
{
	if (i >= poly->nedges)
		return;

	int min, newindex;
	int local = i - 1;
	int l = leftArray(i);
	int r = rightArray(i);

	if (l != -1 && poly->elist[local]->cost > poly->elist[l]->cost)
	{
		min = l;
		newindex = 2 * i;
	}
	else
	{
		min = local;
	}

	if (r != -1 && poly->elist[r]->cost < poly->elist[min]->cost)
	{
		min = r;
		newindex = 2 * i + 1;
	}

	if (min != local)
	{
		swapEdges(local, min);
		heapify(newindex);
	}
}

int simplification(int n)
{
	std::cout << "Simplifying the edges..." << std::endl;

	int start = 0;
	int deleted;

	int nComputations = 0;
	int casesIsZero = 0;

	while (n > 0)
	{
		if (trianglesDeleted >= (poly->ntris / 100 * (100 - 5.0)))
		{
			break;
		}

		for (int i = 0; i < poly->nedges; i++)
		{
			if (poly->elist[i]->cost <= 9999990.99 && poly->elist[i]->deleted == 0)
			{
				Vertex *v1 = poly->elist[i]->verts[0];
				Vertex *v2 = poly->elist[i]->verts[1];
				if (poly->elist[i]->verts[0] == poly->elist[i]->verts[1])
				{
					return 0;
				}

				std::vector<Triangle *> commonTriangles;

				for (int j = 0; j < v1->ntris; j++)
				{
					Triangle *t1 = v1->tris[j];
					for (int k = 0; k < v2->ntris; k++)
					{
						Triangle *t2 = v2->tris[k];
						if (t1 == t2)
						{
							commonTriangles.push_back(t1);
						}
					}
				}

				Triangle *currentT;
				Vertex *v3;
				Edge *e2, *e3;
				Vertex *newV3;

				for (int i = 0; i < commonTriangles.size(); i++)
				{
					Triangle *t_i = commonTriangles[i];

					newV3 = getThirdVertex(t_i, v1, v2);

					t_i->deleted = 1;

					e2 = getLinkingEdge(v1, newV3);
					e3 = getLinkingEdge(v2, newV3);

					v1->tris = NULL;
					v2->tris = NULL;
					newV3->tris = NULL;

					trianglesDeleted++;

					e3->deleted = 1;
					e3->cost = 9999999.99;

					/* Heapify */
					heapify(e3->index + 1);
				}

				Edge *currE;
				for (int i = 0; i < v2->ntris; i++)
				{
					Triangle *t_i = v2->tris[i];
					for (int j = 0; j < 3; j++)
					{
						currE = t_i->edges[j];

						if(currE->verts[0] == v2)
							currE->verts[0] = v1;
						else if(currE->verts[1] == v2)
							currE->verts[1] = v1;
					}
				}
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
					for (int i = 0; i < poly->nedges; i++)
						calculateSolutions(i);
				}

				quicksort(0, poly->nedges - 1); // Quick sort the edges

				edges_simplified = poly->ntris / 12;

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

		void set_view(GLenum mode, Polyhedron * poly)
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

		void set_scene(GLenum mode, Polyhedron * poly)
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

		void display_shape(GLenum mode, Polyhedron * this_poly)
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