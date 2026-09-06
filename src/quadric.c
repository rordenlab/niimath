/////////////////////////////////////////////
//
// Mesh Simplification Tutorial
//
// (C) by Sven Forstmann in 2014
//
// License : MIT
// http://opensource.org/licenses/MIT
//
//https://github.com/sp4cerat/Fast-Quadric-Mesh-Simplification
//
// 5/2016: Chris Rorden created minimal version for OSX/Linux/Windows compile
// 1/2022: Chris Rorden ported from C++ to pure C

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#ifdef _MSC_VER

#else
 #include <unistd.h>
#endif
#include <time.h>
#include "meshtypes.h"
#include "quadric.h"
#include "meshify.h"
#include <float.h> //FLT_EPSILON, DBL_EPSILON
#include <limits.h>

typedef double TSymetricMatrix[10];
struct TRef{
	int tid,tvertex;
};
struct TVertex { vec3d p;int tstart,tcount;TSymetricMatrix q;int border;};
struct TTriangle{
	int v[3];
	double err[4];
	bool dirty, deleted;
	vec3d n;
};

void symMat1(TSymetricMatrix ret, double c){
	for (int i = 0; i < 10; i++)
		ret[i] = c;
} // symMat()

void symMat4(TSymetricMatrix ret, double a,double b,double c,double d){
	ret[0] = a*a; ret[1] = a*b; ret[2] = a*c; ret[3] = a*d;
	ret[4] = b*b; ret[5] = b*c; ret[6] = b*d;
	ret[7] = c*c; ret[8] = c*d;
	ret[9] = d*d;
}// symMat2()

void symMat10(TSymetricMatrix ret, double m11, double m12, double m13, double m14, double m22, double m23, double m24, double m33, double m34, double m44){
	ret[0] = m11; ret[1] = m12; ret[2] = m13; ret[3] = m14;
	ret[4] = m22; ret[5] = m23; ret[6] = m24;
	ret[7] = m33; ret[8] = m34;
	ret[9] = m44;
} // symMat3()

void symMatAdd(TSymetricMatrix ret, TSymetricMatrix n, TSymetricMatrix m) {
	symMat10(ret, n[0]+m[0], n[1]+m[1], n[2]+m[2], n[3]+m[3], n[4]+m[4],
	n[5]+m[5], n[6]+m[6], n[7]+m[7], n[8]+m[8], n[9]+m[9]);
} // symMatAdd()

double symMatDet(TSymetricMatrix m, int a11, int a12, int a13, int a21, int a22, int a23, int a31, int a32, int a33) {
	return m[a11]*m[a22]*m[a33] + m[a13]*m[a21]*m[a32] + m[a12]*m[a23]*m[a31]
	- m[a13]*m[a22]*m[a31] - m[a11]*m[a23]*m[a32]- m[a12]*m[a21]*m[a33];
} // symMatDet()

static inline vec3d ptf(double x, double y, double z) {
	return (vec3d){.x = x, .y = y, .z = z};
}// ptf()

static inline vec3d vCross(vec3d v1, vec3d v2) { //cross-product
	return ptf(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z,
		v1.x * v2.y - v1.y * v2.x);
}

static inline  vec3d vSum(vec3d a, vec3d b){ //add two vectors
	return ptf(a.x+b.x, a.y+b.y, a.z+b.z);
}

static inline vec3d vSubtract(vec3d a, vec3d b){
	return ptf(a.x-b.x, a.y-b.y, a.z-b.z);
}

static inline void vNormalize(vec3d *v){ //make vector unit length
	double len = sqrt( (v->x*v->x) + (v->y*v->y) + (v->z*v->z));
	if (len <= 0) len = 0.001;
	v->x = v->x / len;
	v->y = v->y / len;
	v->z = v->z / len;
}

static inline double vDot (vec3d a,vec3d b){ //dot product
	return a.x*b.x + a.y*b.y + a.z*b.z;
} // vDot()

static inline vec3d vMult(vec3d a, double v){ //multiply
	return ptf(a.x*v, a.y*v, a.z*v);
} // vMult()

double vertex_error(TSymetricMatrix q, double x, double y, double z){
	return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x + q[4]*y*y
	+ 2*q[5]*y*z + 2*q[6]*y + q[7]*z*z + 2*q[8]*z + q[9];
} // vertex_error()

double calculate_error(int id_v1, int id_v2, vec3d *p_result, struct TVertex vertices[]) {
	TSymetricMatrix q;
	symMatAdd(q, vertices[id_v1].q, vertices[id_v2].q);
	int border = vertices[id_v1].border + vertices[id_v2].border;
	double det = symMatDet(q, 0, 1, 2, 1, 4, 5, 2, 5, 7);
	if (( det != 0.0) && ( border == 0)) {
		// q_delta is invertible
		p_result->x = -1.0/det*(symMatDet(q,1, 2, 3, 4, 5, 6, 5, 7 , 8));  // vx = A41/det(q_delta)
		p_result->y =  1.0/det*(symMatDet(q,0, 2, 3, 1, 5, 6, 2, 7 , 8));  // vy = A42/det(q_delta)
		p_result->z = -1.0/det*(symMatDet(q,0, 1, 3, 1, 4, 6, 2, 5,  8));  // vz = A43/det(q_delta)
		return vertex_error(q, p_result->x, p_result->y, p_result->z);
	}
	// det = 0 -> try to find best result
	vec3d p1 = vertices[id_v1].p;
	vec3d p2 = vertices[id_v2].p;
	vec3d p3 = vMult(vSum(p1, p2), 0.5);
	double error1 = vertex_error(q, p1.x,p1.y,p1.z);
	double error2 = vertex_error(q, p2.x,p2.y,p2.z);
	double error3 = vertex_error(q, p3.x,p3.y,p3.z);
	double error = fmin(error1, fmin(error2, error3));
	if (error1 == error) *p_result = p1;
	if (error2 == error) *p_result = p2;
	if (error3 == error) *p_result = p3;
	return error;
}

#define loopi(start_l,end_l) for ( int i=start_l;i<end_l;++i )
#define loopj(start_l,end_l) for ( int j=start_l;j<end_l;++j )
#define loopk(start_l,end_l) for ( int k=start_l;k<end_l;++k )

#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif

void update_mesh(int iteration, struct TTriangle triangles[], struct TVertex vertices[], struct TRef refs[], int *nrefs, int* nTri, int* nVert) {
	if (iteration>0) { // compact triangles
		int dst = 0;
		for (int i = 0; i < *nTri; i++) {
			if(!triangles[i].deleted) {
				triangles[dst] = triangles[i];
				dst = dst + 1;
			}//if not deleted
		} //for each triangle
		*nTri = dst;
		//realloc(triangles, dst); //<- we never resize
	} //if iteration > 0
	loopi(0,*nVert) {
		vertices[i].tstart=0;
		vertices[i].tcount=0;
	}
	loopi(0,*nTri)
		loopj(0,3) vertices[triangles[i].v[j]].tcount++;
	int tstart=0;
	loopi(0,*nVert) {
		vertices[i].tstart=tstart;
		tstart+=vertices[i].tcount;
		vertices[i].tcount=0;
	}
	*nrefs = tstart;   /* the total: appends start past the LAST fan, not inside it (that overwrote
	                      the highest-id vertex's refs and opened holes when it collapsed late) */
	loopi(0,*nTri) {
		struct TTriangle *t=&triangles[i];
		loopj(0,3) {
			struct TVertex* v=&vertices[t->v[j]];
			refs[v->tstart+v->tcount].tid=i;
			refs[v->tstart+v->tcount].tvertex=j;
			v->tcount++;
		}
	}
	if( iteration != 0 ) return;
	// Init Quadrics by Plane & Edge Errors
	//
	// required at the beginning ( iteration == 0 )
	// recomputing during the simplification is not required,
	// but mostly improves the result for closed meshes
	//
	// Identify boundary : vertices[].border=0,1
	//std::vector<int> vcount,vids;
	loopi(0,*nVert)
		vertices[i].border=0;
	int *vids = (int *)malloc(*nVert * sizeof(int));
	int *vcount = (int *)malloc(*nVert * sizeof(int));
	loopi(0,*nVert) {
		int nvcount = 0;
		struct TVertex* v=&vertices[i];
		loopj(0,v->tcount) {
			int k=refs[v->tstart+j].tid;
			struct TTriangle *t=&triangles[k];
			loopk(0,3) {
				int ofs=0,id=t->v[k];
				while(ofs < nvcount) {
					if(vids[ofs]==id)break;
					ofs++;
				}
				if(ofs == nvcount && nvcount < *nVert) {
					vcount[nvcount] = 1;
					vids[nvcount] = id;
					nvcount++;
				}
				else
					vcount[ofs]++;
			}
		}
		loopj(0,nvcount) if(vcount[j]==1)
			vertices[vids[j]].border=1;
	}
	free(vcount);
	free(vids);
	//initialize errors
	loopi(0,*nVert)
		symMat1(vertices[i].q, 0.0);
	loopi(0,*nTri) {
		struct TTriangle *t=&triangles[i];
		vec3d n,p[3];
		loopj(0,3) p[j]=vertices[t->v[j]].p;
		n = vCross(vSubtract(p[1],p[0]),vSubtract(p[2],p[0]));
		vNormalize(&n);
		t->n=n;
		loopj(0,3) {
			TSymetricMatrix q;
			symMat4(q, n.x,n.y,n.z,-vDot(n,p[0]));
			symMatAdd(vertices[t->v[j]].q, vertices[t->v[j]].q, q);
		}
	}
	loopi(0,*nTri) {
		// Calc Edge Error
		struct TTriangle *t=&triangles[i];
		vec3d p;
		loopj(0,3) t->err[j]=calculate_error(t->v[j],t->v[(j+1)%3],&p, vertices);
		t->err[3]=fmin(t->err[0],fmin(t->err[1],t->err[2]));
	}
}

void compact_mesh(struct TTriangle triangles[], struct TVertex vertices[], int* nTri, int* nVert){
		int dst=0;
		loopi(0,*nVert)
			vertices[i].tcount=0;
		loopi(0,*nTri) {
			if(!triangles[i].deleted){
				struct TTriangle t=triangles[i];
				triangles[dst++]=t;
				loopj(0,3)vertices[t.v[j]].tcount=1;
			}
		}
		* nTri = dst;
		dst=0;
		loopi(0, *nVert) {
			if(vertices[i].tcount) {
				vertices[i].tstart=dst;
				vertices[dst].p=vertices[i].p;
				dst++;
			}
		}
		loopi(0,*nTri){
			struct TTriangle *t=&triangles[i];
			loopj(0,3)t->v[j]=vertices[t->v[j]].tstart;
		}
		* nVert = dst;
}

void update_triangles(int i0, struct TVertex* v, bool *deleted, int* deleted_triangles, struct TTriangle triangles[], struct TRef refs[], struct TVertex vertices[], int * nrefs){
	vec3d p;
	loopk(0,v->tcount) {
		struct TRef r=refs[v->tstart+k];
		struct TTriangle *t=&triangles[r.tid];
		if(t->deleted)continue;
		if(deleted[k]) {
			t->deleted=1;
			*deleted_triangles = *deleted_triangles + 1;
			continue;
		}
		{	/* Only the two edges touching the moved vertex change; the opposite edge's endpoints
			   and quadrics are untouched, so its error is still valid (sp4cerat/#46).  One third
			   fewer calculate_error calls in the hot path, bit-identical output. */
			int s = r.tvertex, s2 = (s + 2) % 3;
			t->v[s]=i0;
			t->dirty=1;
			t->err[s]=calculate_error(t->v[s],t->v[(s+1)%3],&p, vertices);
			t->err[s2]=calculate_error(t->v[s2],t->v[s],&p, vertices);
			t->err[3]=fmin(t->err[0],fmin(t->err[1],t->err[2]));
		}
		refs[*nrefs] = r;
		*nrefs = *nrefs + 1;
	}
}

static bool flipped(vec3d p, int i1, const struct TVertex *v0, bool *deleted, struct TTriangle triangles[], struct TRef refs[], struct TVertex vertices[]) {
	loopk(0,v0->tcount) {
		struct TTriangle *t=&triangles[refs[v0->tstart+k].tid];
		if(t->deleted)continue;
		int s=refs[v0->tstart+k].tvertex;
		int id1=t->v[(s+1)%3];
		int id2=t->v[(s+2)%3];
		if(id1==i1 || id2==i1) {// delete ?
			deleted[k]=1;
			continue;
		}
		vec3d d1 = vSubtract(vertices[id1].p, p);
		vNormalize(&d1);
		vec3d d2 = vSubtract(vertices[id2].p, p);
		vNormalize(&d2);
		if(fabs(vDot(d1, d2))>0.999) return true;
		vec3d n = vCross(d1,d2);
		vNormalize(&n);
		deleted[k]=0;
		if (vDot(n, t->n)<0.2) return true;
	}
	return false;
}

/* Link condition (Dey et al. 1999): collapsing edge (i0,i1) keeps the mesh manifold only if the
   one-rings of i0 and i1 share EXACTLY the apex vertices of the faces on that edge -- two for an
   interior edge, one on a boundary.  Any other common neighbour is a thin bridge, and the
   collapse pinches it into a non-manifold edge.  Rings are small (~6), so this is a few dozen
   integer compares per candidate. */
#define RING_MAX 64
static int ring(int v, struct TTriangle triangles[], struct TRef refs[], struct TVertex vertices[], int *out) {
	int n = 0;
	loopk(0,vertices[v].tcount) {
		struct TTriangle *t = &triangles[refs[vertices[v].tstart+k].tid];
		if (t->deleted) continue;
		loopj(0,3) {
			int a = t->v[j], seen = 0;
			if (a == v) continue;
			for (int m = 0; m < n; m++) if (out[m] == a) { seen = 1; break; }
			if (!seen) { if (n == RING_MAX) return -1; out[n++] = a; }
		}
	}
	return n;
}
static bool link_ok(int i0, int i1, struct TTriangle triangles[], struct TRef refs[], struct TVertex vertices[]) {
	int r0[RING_MAX], r1[RING_MAX], n0 = ring(i0, triangles, refs, vertices, r0), n1 = ring(i1, triangles, refs, vertices, r1);
	int common = 0, apex = 0;
	if (n0 < 0 || n1 < 0) return false;   /* absurd valence: refuse rather than guess */
	for (int a = 0; a < n0; a++) for (int b = 0; b < n1; b++) if (r0[a] == r1[b]) { common++; break; }
	loopk(0,vertices[i0].tcount) {
		struct TTriangle *t = &triangles[refs[vertices[i0].tstart+k].tid];
		if (!t->deleted && (t->v[0]==i1 || t->v[1]==i1 || t->v[2]==i1)) apex++;
	}
	if (apex == 2 && vertices[i0].border && vertices[i1].border) return false;   /* a chord: collapsing pinches the boundary loop */
	return common == apex;
}
/* ---- self-intersection guard -----------------------------------------------------------------
 * An edge collapse can drag a triangle through a non-adjacent sheet (thin gyri, touching folds);
 * flipped() bounds normal rotation, which does not see that.  Before accepting a collapse, the
 * new fan is tested against every live triangle sharing a grid cell with it, using Moller's
 * triangle-triangle test; a hit rejects the collapse.  The grid is rebuilt whenever update_mesh
 * compacts the triangle ids and APPENDED to whenever a fan moves in between (stale entries are
 * harmless: deleted faces are skipped, moved faces are just re-tested).  Triangles that share a
 * vertex are adjacent, not crossing, and are skipped. */
/* Cells are CSR (contiguous) at build time, so a walk is sequential; a triangle whose fan moved
   into a NEW cell goes into a small per-cell overflow list (`ohead`/`oitems`, (tri,next) pairs).
   Moves are short -- a collapse point sits within an edge of the fan, cells are several edges
   wide -- so the overflow stays small between rebuilds. */
/* A cell entry carries the triangle's box, so a walk is sequential memory and rejects most entries
   before touching anything per-triangle.  MEASURED (MNI, -q 2): 130 entries visited per query, 22
   pass the box, 6 reach a triangle test; with the box looked up per triangle the walk was three
   scattered loads per entry and the guard cost 1.9 s instead of 1.1 s. */
typedef struct { int t; float bx[6]; } XEntry;
typedef struct { int t, next; float bx[6]; } XOEntry;   /* overflow chain: next is 1-based, 0 ends */
typedef struct { double lo[3], cell; int g[3]; long ncell; int *start; XEntry *items; int *ohead; XOEntry *oitems; int no, ocap; int *stamp, epoch; uint8_t *dead; float *box; } XGrid;

static void xg_cells(XGrid *G, const float *bx, int c0[3], int c1[3]) {
	for (int k = 0; k < 3; k++) {
		int i0 = (int)((bx[k] - G->lo[k]) / G->cell), i1 = (int)((bx[3 + k] - G->lo[k]) / G->cell);
		c0[k] = i0 < 0 ? 0 : (i0 >= G->g[k] ? G->g[k] - 1 : i0);
		c1[k] = i1 < 0 ? 0 : (i1 >= G->g[k] ? G->g[k] - 1 : i1);
	}
}

static void xg_free(XGrid *G) { free(G->start); free(G->items); free(G->ohead); free(G->oitems); free(G->stamp); free(G->dead); free(G->box); memset(G, 0, sizeof *G); }
/* Boxes are float, rounded OUTWARD: a to-nearest cast can move a bound across a cell boundary
   and leave a face unregistered in the one cell it shares with a query. */
static inline float xg_lo(double v) { return nextafterf((float)v, -FLT_MAX); }
static inline float xg_hi(double v) { return nextafterf((float)v, FLT_MAX); }
static void xg_setbox(XGrid *G, struct TTriangle *T, struct TVertex *V, int t) {
	const double *a = &V[T[t].v[0]].p.x, *b = &V[T[t].v[1]].p.x, *c = &V[T[t].v[2]].p.x;
	float *bx = G->box + 6 * t;
	for (int d = 0; d < 3; d++) { bx[d] = xg_lo(fmin(a[d], fmin(b[d], c[d]))); bx[3 + d] = xg_hi(fmax(a[d], fmax(b[d], c[d]))); }
}

static int xg_build(XGrid *G, struct TTriangle *T, struct TVertex *V, int ntri) {
	double hi[3], elen = 0.0;
	int nlive = 0;
	memset(G, 0, sizeof *G);
	if (ntri < 1) return 1;
	G->lo[0] = hi[0] = V[T[0].v[0]].p.x; G->lo[1] = hi[1] = V[T[0].v[0]].p.y; G->lo[2] = hi[2] = V[T[0].v[0]].p.z;
	for (int t = 0; t < ntri; t++) {
		vec3d e;
		if (T[t].deleted) continue;
		nlive++;
		for (int k = 0; k < 3; k++) {
			const double *p = &V[T[t].v[k]].p.x;
			for (int d = 0; d < 3; d++) { if (!(fabs(p[d]) <= DBL_MAX)) return 1; if (p[d] < G->lo[d]) G->lo[d] = p[d]; if (p[d] > hi[d]) hi[d] = p[d]; }   /* NaN/Inf: run unguarded */
		}
		e = vSubtract(V[T[t].v[1]].p, V[T[t].v[0]].p); elen += sqrt(vDot(e, e));
	}
	G->cell = nlive ? 2.0 * elen / nlive : 1.0;
	if (!(G->cell > 0.0)) G->cell = 1.0;
	for (int k = 0; k < 3; k++) { double c = (hi[k] - G->lo[k]) / 256 + 1e-9; if (c > G->cell) G->cell = c; }
	for (int k = 0; k < 3; k++) G->g[k] = (int)((hi[k] - G->lo[k]) / G->cell) + 1;
	G->ncell = (long)G->g[0] * G->g[1] * G->g[2];
	G->start = (int *)calloc((size_t)G->ncell + 1, sizeof(int));
	G->ohead = (int *)calloc((size_t)G->ncell, sizeof(int));
	G->stamp = (int *)calloc((size_t)ntri, sizeof(int));
	G->dead = (uint8_t *)calloc((size_t)ntri, 1);   /* memo of T[u].deleted, so a dead entry costs one byte load */
	G->box = (float *)malloc((size_t)ntri * 6 * sizeof(float));   /* per-triangle AABB, so the walk never touches vertices */
	G->ocap = 1 << 16;
	G->oitems = (XOEntry *)malloc((size_t)G->ocap * sizeof(XOEntry));
	if (!G->start || !G->ohead || !G->stamp || !G->dead || !G->box || !G->oitems) { xg_free(G); return 1; }
	for (int t = 0; t < ntri; t++) if (!T[t].deleted) xg_setbox(G, T, V, t);
	{	long long total = 0;
		for (int t = 0; t < ntri; t++) {
			int c0[3], c1[3];
			if (T[t].deleted) continue;
			xg_cells(G, G->box + 6 * t, c0, c1);
			total += (long long)(c1[0] - c0[0] + 1) * (c1[1] - c0[1] + 1) * (c1[2] - c0[2] + 1);
			for (int z = c0[2]; z <= c1[2]; z++) for (int y = c0[1]; y <= c1[1]; y++) for (int x = c0[0]; x <= c1[0]; x++)
				G->start[1 + x + (long)G->g[0] * (y + (long)G->g[1] * z)]++;
		}
		if (total > INT_MAX) { xg_free(G); return 1; }   /* only a hostile mesh: faces spanning the whole grid */
	}
	for (long c = 0; c < G->ncell; c++) G->start[c + 1] += G->start[c];
	G->items = (XEntry *)malloc((size_t)G->start[G->ncell] * sizeof(XEntry));
	if (!G->items) { xg_free(G); return 1; }
	{	int *fill = (int *)calloc((size_t)G->ncell, sizeof(int));
		if (!fill) { xg_free(G); return 1; }
		for (int t = 0; t < ntri; t++) {
			int c0[3], c1[3];
			const float *bx = G->box + 6 * t;
			if (T[t].deleted) continue;
			xg_cells(G, bx, c0, c1);
			for (int z = c0[2]; z <= c1[2]; z++) for (int y = c0[1]; y <= c1[1]; y++) for (int x = c0[0]; x <= c1[0]; x++) {
				long c = x + (long)G->g[0] * (y + (long)G->g[1] * z);
				XEntry *e = &G->items[G->start[c] + fill[c]++];
				e->t = t; memcpy(e->bx, bx, sizeof e->bx);
			}
		}
		free(fill);
	}
	return 0;
}

/* Is triangle t registered in cell c?  If so refresh its entry's box (the triangle moved) and
   return 1.  Checks the CSR run and the overflow chain. */
static int xg_refresh(XGrid *G, long c, int t, const float *bx) {
	for (int i = G->start[c]; i < G->start[c + 1]; i++) if (G->items[i].t == t) { memcpy(G->items[i].bx, bx, 6 * sizeof(float)); return 1; }
	for (int it = G->ohead[c]; it; it = G->oitems[it - 1].next) if (G->oitems[it - 1].t == t) { memcpy(G->oitems[it - 1].bx, bx, 6 * sizeof(float)); return 1; }
	return 0;
}

/* Would collapsing i1 into i0 at position p make any fan triangle cross a non-adjacent one?
   One grid walk per COLLAPSE over the box that bounds the whole moved fan, collecting candidates
   into `cand`; then every new fan face is tested against that short list.  Walking per fan face
   instead cost 100M list visits for 2.8M real tests on the MNI mesh. */
#define XG_MAXCAND 4096
static bool xg_crosses(XGrid *G, struct TTriangle *T, struct TVertex *V, struct TRef *refs, int i0, int i1, vec3d p, int *cand, float *cbox) {
	int ids[2] = { i0, i1 }, ncand = 0, c0[3], c1[3];
	double bl[3] = { p.x, p.y, p.z }, bh[3] = { p.x, p.y, p.z };
	/* box of the moved fan: p plus every ring vertex of i0 and i1 */
	for (int s = 0; s < 2; s++) loopk(0,V[ids[s]].tcount) {
		struct TTriangle *t = &T[refs[V[ids[s]].tstart+k].tid];
		if (t->deleted) continue;
		for (int j = 0; j < 3; j++) {
			const double *q = &V[t->v[j]].p.x;
			for (int d = 0; d < 3; d++) { if (q[d] < bl[d]) bl[d] = q[d]; if (q[d] > bh[d]) bh[d] = q[d]; }
		}
	}
	{	float fb[6] = { xg_lo(bl[0]), xg_lo(bl[1]), xg_lo(bl[2]), xg_hi(bh[0]), xg_hi(bh[1]), xg_hi(bh[2]) };
		xg_cells(G, fb, c0, c1);
	}
	G->epoch++;
	for (int z = c0[2]; z <= c1[2]; z++) for (int y = c0[1]; y <= c1[1]; y++) for (int x = c0[0]; x <= c1[0]; x++) {
		long c = x + (long)G->g[0] * (y + (long)G->g[1] * z);
		int i = G->start[c], it = G->ohead[c];
		for (;;) {
			int u;
			struct TTriangle *tu;
			const float *ebx, *bx;
			if (i < G->start[c + 1]) { ebx = G->items[i].bx; u = G->items[i++].t; }
			else if (it) { ebx = G->oitems[it - 1].bx; u = G->oitems[it - 1].t; it = G->oitems[it - 1].next; }
			else break;
			/* the entry's box is current in every cell the triangle covers NOW (xg_moved refreshes
			   those); a stale one elsewhere can only reject, and the fresh entry is then found */
			if (ebx[3] < bl[0] || ebx[0] > bh[0] || ebx[4] < bl[1] || ebx[1] > bh[1] || ebx[5] < bl[2] || ebx[2] > bh[2]) continue;
			if (G->dead[u] || G->stamp[u] == G->epoch) continue;
			G->stamp[u] = G->epoch;
			bx = G->box + 6 * u;   /* the true box, for the per-face prefilter below */
			if (bx[3] < bl[0] || bx[0] > bh[0] || bx[4] < bl[1] || bx[1] > bh[1] || bx[5] < bl[2] || bx[2] > bh[2]) continue;
			tu = &T[u];
			if (tu->deleted) { G->dead[u] = 1; continue; }
			if (tu->v[0]==i0 || tu->v[1]==i0 || tu->v[2]==i0 || tu->v[0]==i1 || tu->v[1]==i1 || tu->v[2]==i1) continue;   /* the fan itself */
			if (ncand == XG_MAXCAND) return true;   /* absurdly dense: refuse the collapse rather than miss one */
			memcpy(cbox + 6 * ncand, bx, 6 * sizeof(float));
			cand[ncand++] = u;
		}
	}
	if (!ncand) return false;
	/* Every new fan face against the short candidate list.  The face box is rounded OUTWARD to
	   float so the branch-free compare can only admit more pairs than the exact one, never fewer;
	   this loop is ~100M box tests on MNI and was the guard's largest single cost. */
	for (int s = 0; s < 2; s++) loopk(0,V[ids[s]].tcount) {
		struct TTriangle *t = &T[refs[V[ids[s]].tstart+k].tid];
		int v[3];
		vec3d q[3];
		float fl[3], fh[3];
		if (t->deleted) continue;
		for (int j = 0; j < 3; j++) { v[j] = t->v[j] == i1 ? i0 : t->v[j]; q[j] = (v[j] == i0) ? p : V[v[j]].p; }
		if (v[0] == v[1] || v[1] == v[2] || v[0] == v[2]) continue;   /* the two faces being deleted */
		for (int d = 0; d < 3; d++) {
			const double *a = &q[0].x, *b = &q[1].x, *c = &q[2].x;
			fl[d] = xg_lo(fmin(a[d], fmin(b[d], c[d]))); fh[d] = xg_hi(fmax(a[d], fmax(b[d], c[d])));
		}
		for (int n = 0; n < ncand; n++) {
			const float *cb = cbox + 6 * n;
			struct TTriangle *tu;
			int adj = 0;
			if ((cb[3] < fl[0]) | (cb[0] > fh[0]) | (cb[4] < fl[1]) | (cb[1] > fh[1]) | (cb[5] < fl[2]) | (cb[2] > fh[2])) continue;
			tu = &T[cand[n]];
			for (int a = 0; a < 3 && !adj; a++) for (int b = 0; b < 3; b++) if (tu->v[a] == v[b]) { adj = 1; break; }
			if (adj) continue;
			if (mesh_tri_tri(&q[0].x, &q[1].x, &q[2].x, &V[tu->v[0]].p.x, &V[tu->v[1]].p.x, &V[tu->v[2]].p.x)) return true;
		}
	}
	return false;
}

/* After a collapse the fan around i0 moved: register its faces in any cell they now cover and
   were not in before. */
static int xg_moved(XGrid *G, struct TTriangle *T, struct TVertex *V, struct TRef *refs, int i0) {
	loopk(0,V[i0].tcount) {
		int tid = refs[V[i0].tstart+k].tid, c0[3], c1[3];
		struct TTriangle *t = &T[tid];
		if (t->deleted) continue;
		const float *bx = G->box + 6 * tid;
		xg_setbox(G, T, V, tid);
		xg_cells(G, bx, c0, c1);
		for (int z = c0[2]; z <= c1[2]; z++) for (int y = c0[1]; y <= c1[1]; y++) for (int x = c0[0]; x <= c1[0]; x++) {
			long c = x + (long)G->g[0] * (y + (long)G->g[1] * z);
			XOEntry *o;
			if (xg_refresh(G, c, tid, bx)) continue;
			if (G->no == G->ocap) {
				XOEntry *grown = (XOEntry *)realloc(G->oitems, (size_t)G->ocap * 2 * sizeof(XOEntry));
				if (!grown) return 1;
				G->oitems = grown; G->ocap *= 2;
			}
			o = &G->oitems[G->no]; o->t = tid; o->next = G->ohead[c]; memcpy(o->bx, bx, sizeof o->bx); G->ohead[c] = ++G->no;
		}
	}
	return 0;
}

void laplacian_smooth(vec3d *restrict verts, vec3i *restrict tris, int nvert, int ntri) {
	vec3d* sum = (vec3d*) malloc(nvert * sizeof(vec3d));
	memset(sum, 0, nvert * sizeof(vec3d));
	int* num = (int*) malloc(nvert * sizeof(int));
	memset(num, 0, nvert * sizeof(int));
	loopi(0, ntri) {
		//each point of a triangle has two neighbors:
		int p0 = tris[i].x;
		int p1 = tris[i].y;
		int p2 = tris[i].z;
		num[p0] += 2;
		sum[p0] = vSum(sum[p0], vSum(verts[p1], verts[p2]));
		num[p1] += 2;
		sum[p1] = vSum(sum[p1], vSum(verts[p0], verts[p2]));
		num[p2] += 2;
		sum[p2] = vSum(sum[p2], vSum(verts[p0], verts[p1]));
	}
	loopi(0, nvert) { //mean location of neighbors
		if (num[i] <= 0) continue;
		verts[i].x = sum[i].x / num[i];
		verts[i].y = sum[i].y / num[i];
		verts[i].z = sum[i].z / num[i];
	}
	free(sum);
	free(num);
}

static vec3d face_normal(const vec3d *p, vec3i t) { return vCross(vSubtract(p[t.y], p[t.x]), vSubtract(p[t.z], p[t.x])); }

/* One smoothing step moved every vertex; put back any whose fan folded (a face normal reversed)
   or now crosses a face it does not touch.  A reverted vertex beside moved neighbours can expose
   a new crossing, so repeat; the reverted set only grows, so it ends (MEASURED: 1-2 passes).
   A displacement gate was tried and rejected: crossings come from sub-edge moves where two sheets
   are close, indistinguishable from the bulk. */
static int smooth_revert(vec3d *p, const vec3d *prev, vec3i *tris, int nvert, int ntri) {   /* 0, or -1 when the scan cannot run */
	uint8_t *hit = (uint8_t *)malloc((size_t)ntri), *back = (uint8_t *)calloc((size_t)nvert, 1);
	int changed = 1;
	if (!hit || !back) { free(hit); free(back); return -1; }
	while (changed) {
		/* rescan every pass: a face with one reverted and two moved vertices is a geometry neither
		   mesh had, and a fold between ADJACENT faces is invisible to the crossing test */
		if (mesh_self_intersections(tris, p, ntri, nvert, hit) < 0) { free(hit); free(back); return -1; }
		loopi(0,ntri) if (vDot(face_normal(p, tris[i]), face_normal(prev, tris[i])) <= 0.0) hit[i] = 1;
		changed = 0;
		loopi(0,ntri) if (hit[i]) {
			int v[3] = { tris[i].x, tris[i].y, tris[i].z };
			loopk(0,3) if (!back[v[k]]) { back[v[k]] = 1; p[v[k]] = prev[v[k]]; changed = 1; }
		}
	}
	free(hit); free(back);
	return 0;
}

/* Border vertices (an edge with one face), as a byte per vertex; NULL on OOM. */
static uint8_t *mesh_border(vec3d *verts, vec3i *tris, int nvert, int ntri) {
	struct TVertex* vertices = (struct TVertex*) malloc(nvert * sizeof(struct TVertex));
	struct TTriangle* triangles = (struct TTriangle*) calloc(ntri, sizeof(struct TTriangle));
	struct TRef* refs = (struct TRef*) malloc(ntri * 3 * sizeof(struct TRef));
	uint8_t *border = (uint8_t *)malloc((size_t)nvert);
	int nref = 0, ntriOK = ntri, vertex_count = nvert;
	if (!vertices || !triangles || !refs || !border) { free(vertices); free(triangles); free(refs); free(border); return NULL; }
	loopi(0, nvert) vertices[i].p = verts[i];
	loopi(0, ntri) { triangles[i].v[0] = tris[i].x; triangles[i].v[1] = tris[i].y; triangles[i].v[2] = tris[i].z; }
	update_mesh(0, triangles, vertices, refs, &nref, &ntriOK, &vertex_count);
	loopi(0, nvert) border[i] = (uint8_t)vertices[i].border;
	free(vertices); free(triangles); free(refs);
	return border;
}

void laplacian_smoothHC(vec3d *verts, vec3i *tris, int nvert, int ntri, double alpha, double beta, int iter, bool lockEdges, bool guard) {
	// Laplacian smooth with Humphrey’s Classes to preserve volume: https://doi.org/10.1111/1467-8659.00334
	//  trimesh.smoothing.filter_humphrey(mesh, alpha=0.1, beta=0.5, iterations=10) https://trimsh.org/trimesh.smoothing.html
	double alpha1 = 1.0 - alpha;
	double beta1 = 1.0 - beta;
	vec3d* p = (vec3d*) malloc(nvert * sizeof(vec3d));
	vec3d* q = (vec3d*) malloc(nvert * sizeof(vec3d));
	vec3d* b = (vec3d*) malloc(nvert * sizeof(vec3d));
	vec3d* prev = guard ? (vec3d*) malloc(nvert * sizeof(vec3d)) : NULL;
	/* Locked border vertices are pinned INSIDE every iteration, so their neighbours smooth against
	   the position that ships and the guard checks the geometry that ships.  (Snapping them back
	   only at the end let a fold along a clipped boundary through the guard.) */
	uint8_t *border = lockEdges ? mesh_border(verts, tris, nvert, ntri) : NULL;
	if (!p || !q || !b || (lockEdges && !border)) { free(p); free(q); free(b); free(prev); free(border); return; }
	memcpy(p, verts, nvert * sizeof(vec3d)); //dst,src,n
	loopj(0,iter) {
		memcpy(q, p, nvert * sizeof(vec3d)); //dst,src,n
		if (prev) memcpy(prev, p, nvert * sizeof(vec3d));
		laplacian_smooth(p, tris, nvert, ntri);
		loopi(0,nvert)
			b[i] = vSubtract(p[i], vSum(vMult(verts[i], alpha), vMult(q[i], alpha1)));
		memcpy(q, b, nvert * sizeof(vec3d));
		laplacian_smooth(q, tris, nvert, ntri);
		loopi(0,nvert)
			p[i] = vSubtract(p[i], vSum(vMult(b[i], beta), vMult(q[i], beta1)));
		if (border) loopi(0,nvert) if (border[i]) p[i] = verts[i];
		if (prev && smooth_revert(p, prev, tris, nvert, ntri) < 0) {
			fprintf(stderr, "smooth: the intersection guard cannot run (no memory, or a non-finite coordinate); continuing unguarded\n");
			free(prev); prev = NULL;
		}
	}
	memcpy(verts, p, nvert * sizeof(vec3d));
	free(q); free(b); free(prev); free(border); free(p);
}

void quadric_simplify_mesh(vec3d **vs, vec3i **ts, int* nvert, int *ntri, int target_count, double aggressiveness, bool verbose, bool finishLossless) {
	// init: load vertices
	vec3d *verts = *vs;
	struct TVertex* vertices = (struct TVertex*) malloc(*nvert * sizeof(struct TVertex));
	loopi(0,*nvert)
		vertices[i].p = verts[i];
	free(*vs);
	//init: load triangle faces
	vec3i *tris = *ts;
	struct TTriangle* triangles = (struct TTriangle*) calloc(*ntri, sizeof(struct TTriangle));
	loopi(0,*ntri) {
		triangles[i].v[0] = tris[i].x;
		triangles[i].v[1] = tris[i].y;
		triangles[i].v[2] = tris[i].z;
	}
	free(*ts);
	int nref = 0;
	int refcap = *ntri * 6;   /* each collapse appends the merged fan; update_mesh compacts every 5th lossy iteration, every lossless one */
	struct TRef* refs = (struct TRef*) malloc(refcap * sizeof(struct TRef));
	//init other structures
	bool* deleted0 = (bool*) malloc(*ntri * 3 * sizeof(bool)); //overprovision so we never need to realloc
	bool* deleted1 = (bool*) malloc(*ntri * 3 * sizeof(bool)); //overprovision so we never need to realloc
	XGrid xgrid; memset(&xgrid, 0, sizeof xgrid);
	int *xcand = (int *)malloc(XG_MAXCAND * sizeof(int));
	float *xcbox = (float *)malloc(XG_MAXCAND * 6 * sizeof(float));
	/* The self-intersection guard is a broad-phase walk per collapse, for a handful of rejected
	   collapses (MNI at r=0.2: 11 self-intersecting triangles without it, 0 with it).  It rides the
	   "best" quality level (finishLossless), which is the default; -q 1 drops it for speed. */
	bool guardCrossings = finishLossless && xcand && xcbox;   /* no scratch: run unguarded */
	// main iteration loop
	int deleted_triangles=0;
	int vertex_count = *nvert;
	int triangle_count=*ntri;
	int ntriOK = triangle_count;
	int max_iter = 100;
	bool lossy = true;
	double threshold = DBL_EPSILON, scale;
	{	/* The threshold schedule is a quadric error (length^2); Simplify.h's constants assume
		   unit-scale models.  Scaled by the mean edge length squared it is unit-invariant --
		   MEASURED: the MNI mesh in microns stopped at 1045008 of a 226336-face target at -q 2. */
		double e = 0.0;
		loopi(0, ntriOK) { vec3d d = vSubtract(vertices[triangles[i].v[1]].p, vertices[triangles[i].v[0]].p); e += sqrt(vDot(d, d)); }
		e = ntriOK ? e / ntriOK : 1.0;
		scale = (e > 0.0) ? e * e : 1.0;
	}
	if (target_count >= ntriOK) {
		lossy = false;
		max_iter = 1000;
	}
	int iterationStartCount = 0;
	for (int iteration = 0; iteration < max_iter; iteration ++) {
		if ((lossy) && ((triangle_count-deleted_triangles)<=target_count)) {
			if (!finishLossless) break;
			lossy = false;
			threshold = DBL_EPSILON;
			max_iter = 1000;
		}
		{	/* update_mesh COMPACTS triangle ids, which the intersection grid stores, so the
			   grid is rebuilt whenever it runs: cheap (one pass) beside the collapse loop. */
			bool updated = false;
			if (lossy) {
				//lossy: update mesh once in a while
				if(iteration%5==0) { update_mesh(iteration, triangles, vertices, refs, &nref, &ntriOK, &vertex_count); updated = true; }
				threshold = scale*0.000000001*pow((double)(iteration+3.0), aggressiveness);
			} else {
				if (iterationStartCount == (triangle_count-deleted_triangles)) break;
				//lossless: update mesh constantly
				update_mesh(iteration, triangles, vertices, refs, &nref, &ntriOK, &vertex_count); updated = true;
			}
			if (guardCrossings && updated) {
				xg_free(&xgrid);
				xg_build(&xgrid, triangles, vertices, ntriOK);   /* frees itself on failure: run unguarded */
			}
		}
		iterationStartCount = triangle_count-deleted_triangles;
		// clear dirty flag
		loopi(0,ntriOK)
			triangles[i].dirty=0;
		//
		// All triangles with edges below the threshold will be removed
		//
		// The following numbers works well for most models.
		// If it does not, try to adjust the 3 parameters
		//
		// target number of triangles reached ? Then break
		if ((verbose) && (iteration%5==0))
			printf(" iteration %d - triangles %d threshold %g\n",iteration,triangle_count-deleted_triangles, threshold);
		// remove vertices & mark deleted triangles
		loopi(0,ntriOK) {
			struct TTriangle *t=&triangles[i];
			if(t->err[3]>threshold) continue;
			if(t->deleted) continue;
			if(t->dirty) continue;
			loopj(0,3)if(t->err[j]<threshold) {
				int i0=t->v[j];
				struct TVertex *v0 = &vertices[i0];
				int i1=t->v[(j+1)%3];
				struct TVertex *v1 = &vertices[i1];
				if(v0->border != v1->border)  continue;
				/* refs is a fixed buffer: a fan wider than the slack (valence ~70 on 140 faces) overran
				   it and opened a hole.  Defer the collapse to after the next compaction instead. */
				if(nref + v0->tcount + v1->tcount > refcap) continue;
				// Compute vertex to collapse to
				vec3d p;
				calculate_error(i0,i1,&p, vertices);
				if( !link_ok(i0,i1,triangles,refs,vertices) ) continue;
				if( flipped(p,i1,v0,deleted0, triangles, refs, vertices) ) continue;
				if( flipped(p,i0,v1,deleted1, triangles, refs, vertices) ) continue;
				if( xgrid.start && xg_crosses(&xgrid, triangles, vertices, refs, i0, i1, p, xcand, xcbox) ) continue;
				// not flipped, so remove edge
				v0->p=p;
				symMatAdd(v0->q, v1->q, v0->q); //v0.q=v1.q+v0.q;
				int tstart=nref;
				update_triangles(i0,v0,deleted0,&deleted_triangles, triangles, refs, vertices, &nref);
				update_triangles(i0,v1,deleted1,&deleted_triangles, triangles, refs, vertices, &nref);
				int tcount=nref-tstart;
				if(tcount<=v0->tcount) {
					// save ram
						if(tcount)memcpy(&refs[v0->tstart],&refs[tstart],tcount*sizeof(struct TRef));
				}
				else
					// append
					v0->tstart=tstart;
				v0->tcount=tcount;
				if (xgrid.start && xg_moved(&xgrid, triangles, vertices, refs, i0)) xg_free(&xgrid);   /* OOM: run unguarded */
				break;
			}
			// done?
			//if(triangle_count-deleted_triangles<=target_count) threshold = DBL_EPSILON;
			if((lossy) && ((triangle_count-deleted_triangles)<=target_count)) break;
		} //for each triangle
	} //for each iteration
	free(refs);
	free(deleted0);
	free(deleted1);
	triangle_count = ntriOK;
	// clean up mesh
	xg_free(&xgrid); free(xcand); free(xcbox);
	compact_mesh(triangles, vertices, &triangle_count, &vertex_count);
	*ntri = triangle_count;
	*ts= (vec3i *) malloc(*ntri * sizeof(vec3i));
	tris = *ts;
	loopi(0,*ntri)
		tris[i] = (vec3i){ .x = triangles[i].v[0], .y = triangles[i].v[1], .z = triangles[i].v[2]};
	*nvert = vertex_count;
	*vs = (vec3d *) malloc(*nvert * sizeof(vec3d));
	verts = *vs;
	loopi(0,*nvert)
		verts[i] = (vec3d){ .x = vertices[i].p.x, .y = vertices[i].p.y, .z = vertices[i].p.z};
	free(triangles);
	free(vertices);
} //quadric_simplify_mesh()

void mesh_simplify(vec3d **vs, vec3i **ts, int *nvert, int *ntri, int target_count, int quality, int engine, bool verbose) {
	double aggressiveness = quality == 0 ? 8.0 : quality == 2 ? 5.0 : 7.0;   /* Simplify.h's default is 7 */
	#ifdef HAVE_QUADRIC2
	if (engine) { quadric2_simplify_mesh(vs, ts, nvert, ntri, target_count, verbose, quality > 1); return; }
	#endif
	(void)engine;
	quadric_simplify_mesh(vs, ts, nvert, ntri, target_count, aggressiveness, verbose, quality > 1);
}
