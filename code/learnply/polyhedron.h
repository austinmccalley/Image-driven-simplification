/*

Data structures for learnply

Eugene Zhang 2005

*/

#ifndef __LEARNPLY_H__
#define __LEARNPLY_H__


#include "ply.h"
#include "icVector.H"

const double EPS = 1.0e-6;
const double PI = 3.1415926535898;

/* forward declarations */
class Quad;
class Edge;

class Vertex {
public:
	double x, y, z;			/*coordinates*/
	double vx, vy, vz;		/*vector field*/
	double scalar = 0;		/*scalar field*/

	double R = 0, G = 0, B = 0;		/*color*/

	int index;

	/* The highest contour which this belongs to */
	int level;

	int nquads;
	Quad** quads;
	int max_quads;

	int nedges;
	Edge** edges;
	int max_edges;

	icVector3 normal;
	void* other_props = NULL;
public:
	Vertex(double xx, double yy, double zz) { x = xx; y = yy; z = zz; }
};

class Edge {
public:
	int index;

	Vertex* verts[2];

	icVector3* crossing;

	int nquads;
	Quad** quads;

	double length;

};

class Quad {
public:
	int index;

	Vertex* verts[4];
	Edge* edges[4];

	float area;

	icVector3 normal;
	void* other_props;

	icVector3* crtical = NULL;
	icVector2* singularity = NULL;
};


class Polyhedron {
public:

	Quad** qlist;		/* list of quads */
	int nquads;
	int max_quads;

	Vertex** vlist;    /* list of vertices */
	int nverts;
	int max_verts;

	int maxScalar;
	int minScalar;

	Edge** elist;      /* list of edges */
	int nedges;
	int max_edges;

	icVector3 center;
	double radius;
	double area;

	int selected_quad;
	int selected_vertex;
	unsigned char orientation;  // 0=ccw, 1=cw

	PlyOtherProp* vert_other, * face_other;

	/*constructors*/
	Polyhedron();
	Polyhedron(FILE*);

	/*initialization functions*/
	void create_pointers();
	void average_normals();
	void create_edge(Vertex*, Vertex*);
	void create_edges();
	int face_to_vertex_ref(Quad*, Vertex*);
	void order_vertex_to_quad_ptrs(Vertex*);
	void vertex_to_quad_ptrs();
	void vertex_to_edge_ptrs();
	void calc_bounding_sphere();
	void calc_face_normals_and_area();
	void calc_edge_length();

	/*utilties*/
	Quad* find_common_edge(Quad*, Vertex*, Vertex*);
	Quad* other_quad(Edge*, Quad*);

	/*feel free to add more to help youself*/
	void write_info();
	void write_file(FILE*);

	Vertex* Polyhedron::other_vert(Vertex* vert, Edge* edge) {
		if (edge->verts[0] == vert)
			return edge->verts[1];
		else if (edge->verts[1] == vert)
			return edge->verts[0];
		else
			return NULL;
	}

	Edge* Polyhedron::find_edge(Vertex* v1, Vertex* v2) {
		Edge* edge;
		for (int i = 0; i < v1->nedges; i++)
		{
			edge = v1->edges[i];
			if (other_vert(v1, edge) == v2)
				return edge;
		}
		return NULL;
	}

	Quad* Polyhedron::find_quad(double x, double y) {
		for (int i = 0; i < nquads; i++) {
			Quad* qTemp = qlist[i];
			double x1, y1, x2, y2;

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

			if (x >= x1 && x <= x2 && y >= y1 && y <= y2)
				return qTemp;
		}
		return NULL;
	}

	/*initialization and finalization*/
	void initialize();
	void finalize();
};

#endif /* __LEARNPLY_H__ */

