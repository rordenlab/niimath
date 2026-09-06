#ifndef QUADRIC_H
#define QUADRIC_H

#include "meshtypes.h"

void quadric_simplify_mesh(vec3d **vs, vec3i **ts, int* nvert, int *ntri, int target_count, double agressiveness, bool verbose, bool finishLossless);
#ifdef HAVE_QUADRIC2
void quadric2_simplify_mesh(vec3d **vs, vec3i **ts, int *nvert, int *ntri, int target_count, bool verbose, bool guard);
#endif
/* the one dispatch: quality 0 fast, 1 no guards, 2 best; engine 1 = quadric2 (HAVE_QUADRIC2 builds only) */
void mesh_simplify(vec3d **vs, vec3i **ts, int *nvert, int *ntri, int target_count, int quality, int engine, bool verbose);
void laplacian_smoothHC(vec3d *verts, vec3i *tris, int nvert, int ntri, double alpha, double beta, int iter, bool lockEdges, bool guard);

#endif /* QUADRIC_H */