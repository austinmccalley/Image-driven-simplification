/*

Data structure for I/O of polygonal models

Eugene Zhang 2005

*/

#ifndef __LEARNPLY_IO_H__
#define __LEARNPLY_IO_H__

typedef struct Vertex_io {
	double x, y, z;
	double vx, vy, vz;
	double s;
	void *other_props;       /* other properties */
} Vertex_io;

typedef struct Face_io {
	unsigned char nverts;    /* number of vertex indices in list */
	int *verts;              /* vertex index list */
	void *other_props;       /* other properties */
} Face_io;

char *elem_names[] = { /* list of the kinds of elements in the user's object */
	"vertex", "face"
};

PlyProperty vert_props[] = { /* list of property information for a vertex */
	{"x", Float64, Float64, offsetof(Vertex_io,x), 0, 0, 0, 0},
	{"y", Float64, Float64, offsetof(Vertex_io,y), 0, 0, 0, 0},
	{"z", Float64, Float64, offsetof(Vertex_io,z), 0, 0, 0, 0},
	{"vx", Float64, Float64, offsetof(Vertex_io,vx), 0, 0, 0, 0},
	{"vy", Float64, Float64, offsetof(Vertex_io,vy), 0, 0, 0, 0},
	{"vz", Float64, Float64, offsetof(Vertex_io,vz), 0, 0, 0, 0},
	{"s", Float64, Float64, offsetof(Vertex_io,s), 0, 0, 0, 0},
};

PlyProperty face_props[] = { /* list of property information for a face */
	{"vertex_indices", Int32, Int32, offsetof(Face_io,verts), 1, Uint8, Uint8, offsetof(Face_io,nverts)},
};

#endif /* __LEARNPLY_IO_H__ */

