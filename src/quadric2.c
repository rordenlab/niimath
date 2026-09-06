// quadric2.c - quadric edge-collapse simplification on a half-edge mesh
//
// A ground-up alternative to quadric.c's threshold-sweep simplifier (opt-in build), selectable with
// `-mesh -n 1` so the two can be compared on the same input.  Same quadric error metric
// (Garland & Heckbert 1997); what differs is the machinery around it:
//
//   * a HALF-EDGE mesh, so one-rings, fans and the link condition are pointer walks, and no
//     reference pool ever needs compacting;
//   * an INDEXED PRIORITY QUEUE of edges keyed on collapse cost (updated in place, not lazily:
//     a lazy heap was measured at 87% stale pops), so collapses happen in true cost order -- no
//     threshold schedule, no aggressiveness knob; the run stops within one face of the target (an
//     interior collapse removes two), then only zero-cost edges (exact duplicates) still go;
//   * the link condition (Dey et al. 1999) and a normal-flip check on every collapse, so the
//     output has the input's topology and no inverted faces;
//   * a self-intersection guard: a uniform grid of face boxes updated per collapse and rebuilt
//     when the face count halves, tested with Moller's triangle-triangle test.
//
// Deterministic: single-threaded, one canonical collapse order; the heap breaks ties on edge id.
// Comparator-free heap, so `Q2=1 make wasm` stays fast.

#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "quadric.h"
#include "meshify.h"


typedef double sym10[10];   /* symmetric 4x4: 00 01 02 03 11 12 13 22 23 33 */

typedef struct {
	vec3d p;
	sym10 q;
	int he;            /* one outgoing half-edge, or -1 once dead */
	uint8_t border;
} Vtx;

/* half-edge k of face f is 3f+k, runs v[k] -> v[(k+1)%3]; `to` is stored, `from` is implicit */
typedef struct { int to, twin; } HE;   /* next = 3*(h/3) + (h%3+1)%3, face = h/3 */

static inline int he_next(int h) { return 3 * (h / 3) + (h % 3 + 1) % 3; }
static inline int he_prev(int h) { return 3 * (h / 3) + (h % 3 + 2) % 3; }

typedef struct {
	Vtx *v;  int nv;
	HE *he;  int nf;          /* 3*nf half-edges */
	uint8_t *falive;
	int nlive;                /* live faces */
} Mesh;

static inline int he_from(const Mesh *m, int h) { return m->he[he_prev(h)].to; }

/* ---------------------------------------------------------------------------- geometry */
static inline vec3d vsub(vec3d a, vec3d b) { vec3d r = { a.x - b.x, a.y - b.y, a.z - b.z }; return r; }
static inline vec3d vcross(vec3d a, vec3d b) { vec3d r = { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x }; return r; }
static inline double vdot(vec3d a, vec3d b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static inline vec3d vnorm(vec3d a) { double l = sqrt(vdot(a, a)); if (l > 0) { a.x /= l; a.y /= l; a.z /= l; } return a; }

static void sym_add(sym10 o, const sym10 a, const sym10 b) { for (int i = 0; i < 10; i++) o[i] = a[i] + b[i]; }
static double sym_err(const sym10 q, double x, double y, double z) {
	return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x + q[4]*y*y + 2*q[5]*y*z + 2*q[6]*y + q[7]*z*z + 2*q[8]*z + q[9];
}
static double det3(double a, double b, double c, double d, double e, double f, double g, double h, double i) {
	return a * e * i + d * h * c + g * b * f - a * h * f - g * e * c - d * b * i;
}

/* collapse cost and the point it collapses to; identical rule to quadric.c */
static double q2_cost(const Mesh *m, int a, int b, vec3d *out) {
	sym10 q;
	sym_add(q, m->v[a].q, m->v[b].q);
	if (!m->v[a].border && !m->v[b].border) {
		double det = det3(q[0], q[1], q[2], q[1], q[4], q[5], q[2], q[5], q[7]);
		if (det != 0.0) {
			out->x = -1.0 / det * det3(q[1], q[2], q[3], q[4], q[5], q[6], q[5], q[7], q[8]);
			out->y =  1.0 / det * det3(q[0], q[2], q[3], q[1], q[5], q[6], q[2], q[7], q[8]);
			out->z = -1.0 / det * det3(q[0], q[1], q[3], q[1], q[4], q[6], q[2], q[5], q[8]);
			return sym_err(q, out->x, out->y, out->z);
		}
	}
	{	vec3d p1 = m->v[a].p, p2 = m->v[b].p, p3 = { 0.5 * (p1.x + p2.x), 0.5 * (p1.y + p2.y), 0.5 * (p1.z + p2.z) };
		double e1 = sym_err(q, p1.x, p1.y, p1.z), e2 = sym_err(q, p2.x, p2.y, p2.z), e3 = sym_err(q, p3.x, p3.y, p3.z);
		double e = fmin(e1, fmin(e2, e3));
		*out = (e == e1) ? p1 : (e == e2) ? p2 : p3;
		return e;
	}
}

/* ---------------------------------------------------------------------------- fans */
/* Outgoing half-edges around v, in order; walks both ways so a boundary vertex is complete.
   Returns the count, or -1 if `cap` is too small. */
static int fan(const Mesh *m, int v, int *out, int cap) {
	int n = 0, h = m->v[v].he, start = h;
	if (h < 0) return 0;
	for (;;) {                                   /* counter-clockwise: twin(prev(h)) */
		int t;
		if (n == cap) return -1;
		out[n++] = h;
		t = m->he[he_prev(h)].twin;
		if (t < 0) break;
		h = t;
		if (h == start) return n;
	}
	h = m->he[start].twin;                       /* boundary hit: sweep clockwise from the start */
	while (h >= 0) {
		h = he_next(h);
		if (n == cap) return -1;
		out[n++] = h;
		h = m->he[h].twin;
	}
	return n;
}

#define FAN_MAX 128

/* ---------------------------------------------------------------------------- heap */
/* INDEXED binary heap: one slot per live undirected edge (keyed by its canonical half-edge),
   updated in place.  A lazy-invalidation heap was tried first and 87% of its pops were stale --
   every collapse re-pushes a dozen edges and invalidates a dozen -- which made the heap the
   single largest cost.  A bucket queue (quarter-octave bins, O(1) update) was slower still:
   the cost is the scattered memory traffic of any global order, not the sift.  Comparator-free,
   for the wasm build. */
typedef struct { double cost; int he; } HeapItem;
typedef struct { HeapItem *a; int *pos; int n, cap; } Heap;   /* pos[he] = slot, or -1 */

static inline bool heap_less(const HeapItem *x, const HeapItem *y) {
	return x->cost < y->cost || (x->cost == y->cost && x->he < y->he);   /* tie on id: deterministic */
}
static void heap_up(Heap *h, int i) {
	HeapItem it = h->a[i];
	while (i > 0) {
		int p = (i - 1) / 2;
		if (!heap_less(&it, &h->a[p])) break;
		h->a[i] = h->a[p]; h->pos[h->a[i].he] = i; i = p;
	}
	h->a[i] = it; h->pos[it.he] = i;
}
static void heap_down(Heap *h, int i) {
	HeapItem it = h->a[i];
	for (;;) {
		int c = 2 * i + 1;
		if (c >= h->n) break;
		if (c + 1 < h->n && heap_less(&h->a[c + 1], &h->a[c])) c++;
		if (!heap_less(&h->a[c], &it)) break;
		h->a[i] = h->a[c]; h->pos[h->a[i].he] = i; i = c;
	}
	h->a[i] = it; h->pos[it.he] = i;
}
static int heap_set(Heap *h, int he, double cost) {   /* insert or update */
	int i = h->pos[he];
	if (i < 0) {
		if (h->n == h->cap) return 1;
		i = h->n++;
		h->a[i].he = he; h->a[i].cost = cost; h->pos[he] = i;
		heap_up(h, i);
	} else {
		double old = h->a[i].cost;
		h->a[i].cost = cost;
		if (cost < old) heap_up(h, i); else heap_down(h, i);
	}
	return 0;
}
static void heap_remove(Heap *h, int he) {
	int i = h->pos[he], last;
	if (i < 0) return;
	last = --h->n;
	h->pos[he] = -1;
	if (i == last) return;
	h->a[i] = h->a[last]; h->pos[h->a[i].he] = i;
	if (i > 0 && heap_less(&h->a[i], &h->a[(i - 1) / 2])) heap_up(h, i); else heap_down(h, i);
}
static HeapItem heap_pop(Heap *h) {
	HeapItem top = h->a[0];
	heap_remove(h, top.he);
	return top;
}
static int heap_init(Heap *h, int nhe, int cap) {
	h->a = (HeapItem *)malloc((size_t)cap * sizeof(HeapItem));
	h->pos = (int *)malloc((size_t)nhe * sizeof(int));
	h->n = 0; h->cap = cap;
	if (!h->a || !h->pos) return 1;
	for (int i = 0; i < nhe; i++) h->pos[i] = -1;
	return 0;
}
static void heap_free(Heap *h) { free(h->a); free(h->pos); }

/* canonical half-edge of an undirected edge: the smaller id of the pair */
static inline int canon(const Mesh *m, int h) { int t = m->he[h].twin; return (t >= 0 && t < h) ? t : h; }

static int push_edge(Mesh *m, Heap *heap, int h) {
	vec3d p;
	h = canon(m, h);
	return heap_set(heap, h, q2_cost(m, he_from(m, h), m->he[h].to, &p));
}

/* ------------------------------------------------------------- self-intersection grid */
/* Per-face AABB in a uniform grid; each cell is a small growable array so a face can be
   REMOVED when its fan moves.  Moller's triangle-triangle test decides. */
/* A cell entry carries the face's box: the walk over a cell is then sequential memory.  Looking
   the box up by face id instead was a cache miss per visit and half the guard's whole cost. */
typedef struct { int f; float bx[6]; } Entry;
typedef struct { Entry *a; int n, cap; } Cell;
typedef struct { double lo[3], cell; int g[3]; long ncell; Cell *c; float *box; int *stamp, epoch; } Grid;   /* stamp: a face is tested once per query */

static void grid_free(Grid *G) {
	if (G->c) for (long i = 0; i < G->ncell; i++) free(G->c[i].a);
	free(G->c); free(G->box); free(G->stamp); memset(G, 0, sizeof *G);
}
static void grid_range(const Grid *G, const float *bx, int c0[3], int c1[3]) {
	for (int k = 0; k < 3; k++) {
		int i0 = (int)((bx[k] - G->lo[k]) / G->cell), i1 = (int)((bx[3 + k] - G->lo[k]) / G->cell);
		c0[k] = i0 < 0 ? 0 : (i0 >= G->g[k] ? G->g[k] - 1 : i0);
		c1[k] = i1 < 0 ? 0 : (i1 >= G->g[k] ? G->g[k] - 1 : i1);
	}
}
static void face_box(const Mesh *m, int f, float *bx) {
	const double *p[3] = { &m->v[m->he[3*f].to].p.x, &m->v[m->he[3*f+1].to].p.x, &m->v[m->he[3*f+2].to].p.x };
	for (int d = 0; d < 3; d++) { bx[d] = nextafterf((float)fmin(p[0][d], fmin(p[1][d], p[2][d])), -FLT_MAX); bx[3+d] = nextafterf((float)fmax(p[0][d], fmax(p[1][d], p[2][d])), FLT_MAX); }   /* outward */
}
static int cell_add(Cell *c, int f, const float *bx) {
	if (c->n == c->cap) { int cap = c->cap ? 2 * c->cap : 8; Entry *a = (Entry *)realloc(c->a, (size_t)cap * sizeof(Entry)); if (!a) return 1; c->a = a; c->cap = cap; }
	c->a[c->n].f = f; memcpy(c->a[c->n].bx, bx, sizeof c->a[c->n].bx); c->n++;
	return 0;
}
static void cell_del(Cell *c, int f) { for (int i = 0; i < c->n; i++) if (c->a[i].f == f) { c->a[i] = c->a[--c->n]; return; } }
static int grid_insert(Grid *G, const Mesh *m, int f) {
	int c0[3], c1[3];
	face_box(m, f, G->box + 6 * f);
	grid_range(G, G->box + 6 * f, c0, c1);
	for (int z = c0[2]; z <= c1[2]; z++) for (int y = c0[1]; y <= c1[1]; y++) for (int x = c0[0]; x <= c1[0]; x++)
		if (cell_add(&G->c[x + (long)G->g[0] * (y + (long)G->g[1] * z)], f, G->box + 6 * f)) return 1;
	return 0;
}
static void grid_remove(Grid *G, int f) {
	int c0[3], c1[3];
	grid_range(G, G->box + 6 * f, c0, c1);
	for (int z = c0[2]; z <= c1[2]; z++) for (int y = c0[1]; y <= c1[1]; y++) for (int x = c0[0]; x <= c1[0]; x++)
		cell_del(&G->c[x + (long)G->g[0] * (y + (long)G->g[1] * z)], f);
}
static int grid_build(Grid *G, const Mesh *m) {
	double hi[3], elen = 0.0;
	memset(G, 0, sizeof *G);
	if (m->nv < 1) return 1;
	G->lo[0] = hi[0] = m->v[0].p.x; G->lo[1] = hi[1] = m->v[0].p.y; G->lo[2] = hi[2] = m->v[0].p.z;
	for (int i = 0; i < m->nv; i++) {
		const double *p = &m->v[i].p.x;
		for (int d = 0; d < 3; d++) { if (!(fabs(p[d]) <= DBL_MAX)) return 1; if (p[d] < G->lo[d]) G->lo[d] = p[d]; if (p[d] > hi[d]) hi[d] = p[d]; }
	}
	for (int f = 0; f < m->nf; f++) if (m->falive[f]) { vec3d e = vsub(m->v[m->he[3*f].to].p, m->v[m->he[3*f+2].to].p); elen += sqrt(vdot(e, e)); }
	G->cell = m->nlive ? 3.0 * elen / m->nlive : 1.0;
	if (!(G->cell > 0.0)) G->cell = 1.0;
	for (int k = 0; k < 3; k++) { double c = (hi[k] - G->lo[k]) / 256 + 1e-9; if (c > G->cell) G->cell = c; }
	for (int k = 0; k < 3; k++) G->g[k] = (int)((hi[k] - G->lo[k]) / G->cell) + 1;
	G->ncell = (long)G->g[0] * G->g[1] * G->g[2];
	G->c = (Cell *)calloc((size_t)G->ncell, sizeof(Cell));
	G->box = (float *)malloc((size_t)m->nf * 6 * sizeof(float));
	G->stamp = (int *)calloc((size_t)m->nf, sizeof(int));
	if (!G->c || !G->box || !G->stamp) { grid_free(G); return 1; }
	for (int f = 0; f < m->nf; f++) if (m->falive[f] && grid_insert(G, m, f)) { grid_free(G); return 1; }
	return 0;
}

/* Would the fan around `a` and `b`, with both at p, cross any face outside that fan?
   Candidates come from the cells under the union box of the moved fan; each surviving fan face
   is then tested only against candidates whose box overlaps ITS box.  Without that per-face
   prefilter the guard did twenty times the triangle tests it needed. */
typedef struct { vec3d q[3]; int v[3]; float bx[6]; } FanFace;
static bool grid_crosses(const Grid *G, const Mesh *m, int a, int b, vec3d p, const int *fa, int na, const int *fb, int nb) {
	FanFace ff[2 * FAN_MAX];
	int nff = 0, c0[3], c1[3];
	float box[6] = { (float)p.x, (float)p.y, (float)p.z, (float)p.x, (float)p.y, (float)p.z };
	for (int s = 0; s < 2; s++) {
		const int *fh = s ? fb : fa; int n = s ? nb : na;
		for (int k = 0; k < n; k++) {
			int f = fh[k] / 3;
			FanFace *F = &ff[nff];
			for (int j = 0; j < 3; j++) { F->v[j] = m->he[3*f+j].to; if (F->v[j] == b) F->v[j] = a; F->q[j] = (F->v[j] == a) ? p : m->v[F->v[j]].p; }
			if (F->v[0]==F->v[1]||F->v[1]==F->v[2]||F->v[0]==F->v[2]) continue;   /* a face being deleted */
			for (int d = 0; d < 3; d++) {
				const double *x = &F->q[0].x, *y = &F->q[1].x, *z = &F->q[2].x;
				F->bx[d] = nextafterf((float)fmin(x[d], fmin(y[d], z[d])), -FLT_MAX); F->bx[3+d] = nextafterf((float)fmax(x[d], fmax(y[d], z[d])), FLT_MAX);
				if (F->bx[d] < box[d]) box[d] = F->bx[d];
				if (F->bx[3+d] > box[3+d]) box[3+d] = F->bx[3+d];
			}
			nff++;
		}
	}
	grid_range(G, box, c0, c1);
	((Grid *)G)->epoch++;
	for (int z = c0[2]; z <= c1[2]; z++) for (int y = c0[1]; y <= c1[1]; y++) for (int x = c0[0]; x <= c1[0]; x++) {
		const Cell *c = &G->c[x + (long)G->g[0] * (y + (long)G->g[1] * z)];
		for (int i = 0; i < c->n; i++) {
			int u = c->a[i].f, uv[3];
			const float *ub = c->a[i].bx;
			if (G->stamp[u] == G->epoch) continue;
			G->stamp[u] = G->epoch;
			if (ub[3] < box[0] || ub[0] > box[3] || ub[4] < box[1] || ub[1] > box[4] || ub[5] < box[2] || ub[2] > box[5]) continue;
			uv[0] = m->he[3*u].to; uv[1] = m->he[3*u+1].to; uv[2] = m->he[3*u+2].to;
			if (uv[0]==a||uv[1]==a||uv[2]==a||uv[0]==b||uv[1]==b||uv[2]==b) continue;   /* the fan itself */
			for (int k = 0; k < nff; k++) {
				const FanFace *F = &ff[k];
				int adj = 0;
				if (ub[3] < F->bx[0] || ub[0] > F->bx[3] || ub[4] < F->bx[1] || ub[1] > F->bx[4] || ub[5] < F->bx[2] || ub[2] > F->bx[5]) continue;
				for (int j = 0; j < 3 && !adj; j++) for (int l = 0; l < 3; l++) if (uv[l] == F->v[j]) { adj = 1; break; }
				if (adj) continue;
				if (mesh_tri_tri(&F->q[0].x, &F->q[1].x, &F->q[2].x, &m->v[uv[0]].p.x, &m->v[uv[1]].p.x, &m->v[uv[2]].p.x)) return true;
			}
		}
	}
	return false;
}

/* ---------------------------------------------------------------------------- checks */
/* Link condition: the one-rings of a and b share exactly the apexes of the faces on edge (a,b). */
static bool link_ok(const Mesh *m, int a, int b, const int *fa, int na, const int *fb, int nb) {
	int ra[2 * FAN_MAX], rb[2 * FAN_MAX], n_a = 0, n_b = 0, common = 0, apex = 0;
	/* ring = both other vertices of every fan face (an outgoing edge alone misses the last ring
	   vertex of a boundary vertex); duplicates are skipped below */
	for (int k = 0; k < na; k++) { ra[n_a++] = m->he[fa[k]].to; ra[n_a++] = m->he[he_next(fa[k])].to; }
	for (int k = 0; k < nb; k++) { rb[n_b++] = m->he[fb[k]].to; rb[n_b++] = m->he[he_next(fb[k])].to; }
	for (int i = 0; i < n_a; i++) {
		int dup = 0;
		if (ra[i] == a || ra[i] == b) continue;
		for (int j = 0; j < i; j++) if (ra[j] == ra[i]) { dup = 1; break; }
		if (dup) continue;
		for (int j = 0; j < n_b; j++) if (rb[j] == ra[i]) { common++; break; }
	}
	for (int k = 0; k < na; k++) if (m->he[fa[k]].to == b) apex += (m->he[fa[k]].twin >= 0) ? 2 : 1;
	return common == apex;
}

/* Would moving a (and b) to p flip or degenerate any surviving face of the fan? */
static bool flips(const Mesh *m, int a, int b, vec3d p, const int *fh, int n) {
	for (int k = 0; k < n; k++) {
		int f = fh[k] / 3, fv[3];
		vec3d q[3], d1, d2, nn, old;
		for (int j = 0; j < 3; j++) fv[j] = m->he[3*f+j].to;
		if (fv[0]==a||fv[1]==a||fv[2]==a) { if (fv[0]==b||fv[1]==b||fv[2]==b) continue; }   /* deleted with the edge */
		old = vnorm(vcross(vsub(m->v[fv[1]].p, m->v[fv[0]].p), vsub(m->v[fv[2]].p, m->v[fv[0]].p)));
		for (int j = 0; j < 3; j++) q[j] = (fv[j] == a || fv[j] == b) ? p : m->v[fv[j]].p;
		{	int j = (fv[0]==a||fv[0]==b) ? 0 : (fv[1]==a||fv[1]==b) ? 1 : 2;
			d1 = vnorm(vsub(q[(j+1)%3], q[j])); d2 = vnorm(vsub(q[(j+2)%3], q[j]));
		}
		if (fabs(vdot(d1, d2)) > 0.999) return true;       /* needle */
		nn = vnorm(vcross(d1, d2));
		if (vdot(nn, old) < 0.2) return true;                /* flipped or nearly */
	}
	return false;
}

/* ---------------------------------------------------------------------------- collapse */
/* Merge b into a at p along half-edge h (a->b). */
static void collapse(Mesh *m, int h, vec3d p, const int *fb, int nb) {
	int a = he_from(m, h), b = m->he[h].to, t = m->he[h].twin, removed = 0;
	int h1 = he_next(h), h2 = he_prev(h), t1 = m->he[h1].twin, t2 = m->he[h2].twin;   /* face f1: a->b->c->a */
	int c = m->he[h1].to;
	/* every half-edge that pointed at b now points at a */
	for (int k = 0; k < nb; k++) m->he[he_prev(fb[k])].to = a;
	/* f1 goes; its two outer half-edges become each other's twins */
	if (t1 >= 0) m->he[t1].twin = t2;
	if (t2 >= 0) m->he[t2].twin = t1;
	m->falive[h / 3] = 0; removed++;
	if (m->v[c].he == h2) m->v[c].he = t1 >= 0 ? t1 : (t2 >= 0 ? he_next(t2) : -1);   /* c's outgoing h2 (c->a) died; on a boundary c still has faces across t2 */
	if (t >= 0) {                              /* face f2: b->a->d->b */
		int h3 = he_next(t), h4 = he_prev(t), t3 = m->he[h3].twin, t4 = m->he[h4].twin, d = m->he[h3].to;
		if (t3 >= 0) m->he[t3].twin = t4;
		if (t4 >= 0) m->he[t4].twin = t3;
		m->falive[t / 3] = 0; removed++;
		if (m->v[d].he == h4) m->v[d].he = t3 >= 0 ? t3 : (t4 >= 0 ? he_next(t4) : -1);
		m->v[a].he = (t4 >= 0) ? t4 : (t2 >= 0 ? t2 : -1);   /* t4 is a->d after relabelling */
	} else
		m->v[a].he = (t2 >= 0) ? t2 : -1;
	if (m->v[a].he < 0) {   /* a's remaining faces are all in b's old fan */
		for (int k = 0; k < nb && m->v[a].he < 0; k++)
			if (m->falive[fb[k] / 3]) m->v[a].he = fb[k];
	}
	m->v[a].p = p;
	sym_add(m->v[a].q, m->v[a].q, m->v[b].q);
	m->v[a].border |= m->v[b].border;
	m->v[b].he = -1;
	m->nlive -= removed;
}

/* ---------------------------------------------------------------------------- build */
/* Returns 1 on OOM, 2 if the input is not a manifold triangle mesh: an edge used by more than two
   faces or twice in the same direction, a vertex whose faces do not form one fan (or a fan wider
   than FAN_MAX), or a zero-area face.  The half-edge structure cannot represent those, and a
   collapse over them would corrupt the connectivity silently. */
typedef struct { uint64_t key; int he; } EdgeSlot;
static int build(Mesh *m, const vec3d *pts, const vec3i *tris, int npt, int ntri) {
	size_t cap = 1;
	EdgeSlot *tab;
	m->nv = npt; m->nf = ntri; m->nlive = ntri;
	m->v = (Vtx *)calloc((size_t)npt, sizeof(Vtx));
	m->he = (HE *)malloc((size_t)ntri * 3 * sizeof(HE));
	m->falive = (uint8_t *)malloc((size_t)ntri);
	while (cap < (size_t)ntri * 6) cap <<= 1;
	tab = (EdgeSlot *)calloc(cap, sizeof(EdgeSlot));
	if (!m->v || !m->he || !m->falive || !tab) { free(tab); return 1; }
	for (int i = 0; i < npt; i++) { m->v[i].p = pts[i]; m->v[i].he = -1; }
	for (int f = 0; f < ntri; f++) {
		int v[3] = { tris[f].x, tris[f].y, tris[f].z };
		m->falive[f] = 1;
		for (int k = 0; k < 3; k++) {
			int h = 3 * f + k, from = v[k], to = v[(k + 1) % 3];
			uint64_t key = ((uint64_t)from << 32) | (uint32_t)to, rkey = ((uint64_t)to << 32) | (uint32_t)from;
			size_t i;
			m->he[h].to = to; m->he[h].twin = -1;
			if (from == to || (unsigned)from >= (unsigned)npt || (unsigned)to >= (unsigned)npt) { free(tab); return 2; }
			if (m->v[from].he < 0) m->v[from].he = h;
			/* pair with the reverse edge if it is already in; otherwise register ours */
			i = (size_t)((rkey * 0x9E3779B97F4A7C15ull) >> 20) & (cap - 1);
			while (tab[i].key && tab[i].key != rkey + 1) i = (i + 1) & (cap - 1);
			if (tab[i].key) { if (m->he[tab[i].he].twin >= 0) { free(tab); return 2; } m->he[h].twin = tab[i].he; m->he[tab[i].he].twin = h; continue; }
			i = (size_t)((key * 0x9E3779B97F4A7C15ull) >> 20) & (cap - 1);
			while (tab[i].key && tab[i].key != key + 1) i = (i + 1) & (cap - 1);
			if (tab[i].key) { free(tab); return 2; }
			tab[i].key = key + 1; tab[i].he = h;
		}
	}
	free(tab);
	/* plane quadrics, border flags, and a boundary-first outgoing edge for boundary vertices */
	for (int f = 0; f < ntri; f++) {
		int v[3] = { tris[f].x, tris[f].y, tris[f].z };
		vec3d n = vnorm(vcross(vsub(pts[v[1]], pts[v[0]]), vsub(pts[v[2]], pts[v[0]])));
		double d = -vdot(n, pts[v[0]]);
		if (!(vdot(n, n) > 0.0)) return 2;   /* zero-area face: no plane, no quadric */
		sym10 q = { n.x*n.x, n.x*n.y, n.x*n.z, n.x*d, n.y*n.y, n.y*n.z, n.y*d, n.z*n.z, n.z*d, d*d };
		for (int k = 0; k < 3; k++) {
			sym_add(m->v[v[k]].q, m->v[v[k]].q, q);
			if (m->he[3 * f + k].twin < 0) { m->v[v[k]].border = 1; m->v[v[(k + 1) % 3]].border = 1; m->v[v[k]].he = 3 * f + k; }
		}
	}
	/* vertex manifoldness: every face at v must be reachable by walking around v from v.he.  A
	   bow-tie (two fans touching at one vertex) passes the edge tests above and would leave the
	   collapse machinery working on half a neighbourhood. */
	{	int *cnt = (int *)calloc((size_t)npt, sizeof(int)), fh[FAN_MAX], bad = 0;
		if (!cnt) return 1;
		for (int f = 0; f < ntri; f++) { cnt[tris[f].x]++; cnt[tris[f].y]++; cnt[tris[f].z]++; }
		for (int i = 0; i < npt && !bad; i++) if (cnt[i] && fan(m, i, fh, FAN_MAX) != cnt[i]) bad = 1;
		free(cnt);
		if (bad) return 2;
	}
	return 0;
}

/* ---------------------------------------------------------------------------- driver */
void quadric2_simplify_mesh(vec3d **vs, vec3i **ts, int *nvert, int *ntri, int target_count, bool verbose, bool guard) {
	Mesh m;
	Heap heap;
	Grid grid;
	int fa[FAN_MAX], fb[FAN_MAX], oom = 0, rejected = 0, collapsed = 0;
	memset(&grid, 0, sizeof grid);
	{	int rc = build(&m, *vs, *ts, *nvert, *ntri);
		if (rc) { fprintf(stderr, rc == 2 ? "quadric2: input is not a manifold triangle mesh (or a vertex has over 128 faces); left unsimplified\n" : "quadric2: out of memory\n"); free(m.v); free(m.he); free(m.falive); return; }
	}
	{	int nedge = 0;   /* one queue slot per undirected edge; pos stays per half-edge */
		for (int h = 0; h < 3 * m.nf; h++) if (canon(&m, h) == h) nedge++;
		memset(&heap, 0, sizeof heap);
		if (heap_init(&heap, 3 * m.nf, nedge)) oom = 1;
	}
	for (int h = 0; !oom && h < 3 * m.nf; h++) if (canon(&m, h) == h && push_edge(&m, &heap, h)) { oom = 1; break; }
	int grid_faces = m.nlive;   /* live faces when the grid was last built */
	if (guard && grid_build(&grid, &m)) { if (verbose) fprintf(stderr, "quadric2: no memory for the intersection grid, running unguarded\n"); }
	if (target_count < 0) target_count = 0;
	while (!oom && heap.n) {
		HeapItem it;
		/* The cell size is tied to the mean edge length at build time; as the mesh coarsens the
		   faces outgrow their cells and every fan spans dozens.  Rebuild once the face count has
		   halved -- a handful of O(n) passes over the whole run, each cheaper than the walks it
		   saves.  MEASURED: 54M cell visits per run without this. */
		if (grid.c && 2 * m.nlive < grid_faces) {
			grid_free(&grid);
			if (grid_build(&grid, &m)) { if (verbose) fprintf(stderr, "quadric2: grid rebuild failed, continuing unguarded\n"); }
			grid_faces = m.nlive;
		}
		it = heap_pop(&heap);
		if (!(it.cost <= DBL_MAX)) continue;   /* a NaN cost (1/det overflow) would pass every later test */
		int h = it.he, a, b, na, nb;
		vec3d p;
		if (m.nlive <= target_count && it.cost > 0.0) break;       /* target reached; zero-cost edges still go */
		a = he_from(&m, h); b = m.he[h].to;
		if (m.v[a].border != m.v[b].border) continue;
		if (m.v[a].border && m.he[h].twin >= 0) continue;   /* a chord between two boundary points: collapsing pinches the loop */
		na = fan(&m, a, fa, FAN_MAX); nb = fan(&m, b, fb, FAN_MAX);
		if (na < 0 || nb < 0) continue;                              /* absurd valence: leave it */
		if (!link_ok(&m, a, b, fa, na, fb, nb)) { rejected++; continue; }
		q2_cost(&m, a, b, &p);
		if (flips(&m, a, b, p, fa, na) || flips(&m, a, b, p, fb, nb)) { rejected++; continue; }
		/* Not gated on cost or displacement: MEASURED on mni.mz3, two of three crossings arise from
		   moves under 5% of an edge, indistinguishable from the cheap bulk.  Every collapse is tested (coplanar overlap excepted). */
		if (grid.c && grid_crosses(&grid, &m, a, b, p, fa, na, fb, nb)) { rejected++; continue; }
		if (grid.c) { for (int k = 0; k < na; k++) grid_remove(&grid, fa[k] / 3); for (int k = 0; k < nb; k++) grid_remove(&grid, fb[k] / 3); }
		/* Only the edges that DIE leave the heap: the collapsed edge and the two edges of each
		   deleted face, whose survivors get re-paired under a possibly different canonical id.
		   Every other edge around a or b keeps its id and is updated in place afterwards. */
		{	int t = m.he[h].twin, dying[6] = { h, he_next(h), he_prev(h), -1, -1, -1 }, nd = 3;
			if (t >= 0) { dying[3] = t; dying[4] = he_next(t); dying[5] = he_prev(t); nd = 6; }
			for (int k = 0; k < nd; k++) heap_remove(&heap, canon(&m, dying[k]));
		}
		collapse(&m, h, p, fb, nb);
		collapsed++;
		na = fan(&m, a, fa, FAN_MAX);
		if (na < 0) { fprintf(stderr, "quadric2: a merged fan exceeds %d faces; result is partially simplified\n", FAN_MAX); break; }
		for (int k = 0; k < na; k++) {
			if (grid.c && grid_insert(&grid, &m, fa[k] / 3)) { oom = 1; break; }
			if (push_edge(&m, &heap, fa[k])) { oom = 1; break; }
		}
	}
	if (oom) fprintf(stderr, "quadric2: out of memory; result is partially simplified\n");
	/* compact */
	{	int *remap = (int *)malloc((size_t)m.nv * sizeof(int)), nv2 = 0, nf2 = 0;
		vec3d *pts = (vec3d *)malloc((size_t)m.nv * sizeof(vec3d));
		vec3i *tris = (vec3i *)malloc((size_t)m.nlive * sizeof(vec3i));
		if (remap && pts && tris) {
			for (int i = 0; i < m.nv; i++) remap[i] = -1;
			for (int f = 0; f < m.nf; f++) if (m.falive[f]) for (int k = 0; k < 3; k++) {
				int v = m.he[3 * f + k].to;
				if (remap[v] < 0) { remap[v] = nv2; pts[nv2++] = m.v[v].p; }
			}
			for (int f = 0; f < m.nf; f++) if (m.falive[f]) {
				vec3i t = { remap[m.he[3*f+2].to], remap[m.he[3*f].to], remap[m.he[3*f+1].to] };   /* from, to order */
				tris[nf2++] = t;
			}
			free(*vs); free(*ts);
			*vs = pts; *ts = tris; *nvert = nv2; *ntri = nf2;
		} else { free(pts); free(tris); fprintf(stderr, "quadric2: out of memory compacting\n"); }
		free(remap);
	}
	if (verbose) printf(" quadric2: %d collapses, %d rejected (link/flip%s)\n", collapsed, rejected, grid.c ? "/crossing" : "");
	grid_free(&grid);
	heap_free(&heap); free(m.v); free(m.he); free(m.falive);
}
