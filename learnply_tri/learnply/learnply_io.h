/*

Data structure for I/O of polygonal models

Eugene Zhang 2005

*/

#ifndef __LEARNPLY_IO_H__
#define __LEARNPLY_IO_H__

typedef struct Vertex_io {
  float x,y,z;
  float vx, vy, vz;
  float t00, t01, t10, t11;
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
  {"x", Float32, Float32, offsetof(Vertex_io,x), 0, 0, 0, 0},
  {"y", Float32, Float32, offsetof(Vertex_io,y), 0, 0, 0, 0},
  {"z", Float32, Float32, offsetof(Vertex_io,z), 0, 0, 0, 0},
  {"vx", Float32, Float32, offsetof(Vertex_io,vx), 0, 0, 0, 0},
  {"vy", Float32, Float32, offsetof(Vertex_io,vy), 0, 0, 0, 0},
  {"vz", Float32, Float32, offsetof(Vertex_io,vz), 0, 0, 0, 0},
  {"t00", Float32, Float32, offsetof(Vertex_io,t00), 0, 0, 0, 0},
  {"t01", Float32, Float32, offsetof(Vertex_io,t01), 0, 0, 0, 0},
  {"t10", Float32, Float32, offsetof(Vertex_io,t10), 0, 0, 0, 0},
  {"t11", Float32, Float32, offsetof(Vertex_io,t11), 0, 0, 0, 0}
};

PlyProperty face_props[] = { /* list of property information for a face */
  {"vertex_indices", Int32, Int32, offsetof(Face_io,verts), 1, Uint8, Uint8, offsetof(Face_io,nverts)},
};

#endif /* __LEARNPLY_IO_H__ */

