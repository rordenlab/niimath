#ifndef MESHIFY_H
#define MESHIFY_H

#include <stdbool.h>
#include <stdint.h>
#include "meshtypes.h"

void strip_ext(char *fname);
int save_mz3(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt, bool isGz, float *perVertexScalar);
int save_mesh(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt, bool isGz);
int meshify(float * img, short dim[3], int originalMC, float isolevel, vec3i **t, vec3d **p, int *nt, int *np, bool preSmooth, bool onlyLargest, bool fillBubbles, bool verbose);
void apply_sform(vec3i *t, vec3d *p, int nt, int np, float srow_x[4], float srow_y[4], float srow_z[4]);
int mesh_report(vec3i *tris, vec3d *pts, int ntri, int npt, const char *label);
/* Moller 1997 triangle-triangle overlap; coplanar pairs are reported as NOT crossing */
int mesh_tri_tri(const double *v0, const double *v1, const double *v2, const double *u0, const double *u1, const double *u2);
int mesh_self_intersections(vec3i *tris, vec3d *pts, int ntri, int npt, uint8_t *hit);
double clockMsec(void);
long timediff(double startTimeMsec, double endTimeMsec);

#endif /* MESHIFY_H */