#ifndef SKULLSTRIP_H
#define SKULLSTRIP_H

#include "nifti_io.h"

#ifdef __cplusplus
extern "C" {
#endif

// AFNI-style surface skull stripping (-skullstrip).
//
// LICENCE BOUNDARY -- read this before touching anything here.
//
// Spatial normalisation is ADAPTED from AFNI's thd_brainormalize.c + the mask helpers in
// thd_automask.c. Those files are NIH public domain: no MCW GPL-2 header, and they do not
// use 3DEdge. Attribution is preserved below even though it is legally optional.
//
// The deformation and touchup stages are ALSO adaptations of public-domain AFNI
// (SUMA_BrainWrap.c). The basis is affirmative, not merely the absence of a notice:
// AFNI's LICENSE.txt states the whole tree is a "United States Government Work ... cannot
// be copyrighted" apart from a listed set of exceptions, states that "contributions
// without explicit licensing will be assumed to be entered into the public domain", and
// doc/README/README.copyright dates the rule to work "after 15 Jan 2001" -- this file's
// first commit is 2004-12-30 by an NIH author. A US Government work is not copyrightable
// (17 U.S.C. 105), and GPL is a copyright licence, so there is nothing there to license.
// The full citation list and the per-file hashes are in the manifest.
// The 3DEdge dependency is three call sites, ALL inside SUMA_3dedge3, and OUR contract
// (-no_use_edge) never calls it.
//
// THE ONE CARVE-OUT: SUMA_3dedge3 (SUMA_BrainWrap.c, ~line 5739 onward) wraps
// Malandain's GPL-3.0 Extract_Gradient_Maxima_3D. Do not read or adapt it. Anything
// requiring the edge volume (Milestone 4, the -use_edge default) is out of scope.
//
// The surface PRIMITIVES below (icosphere, adjacency, intersection, rasterisation) are
// genuinely clean-room -- written before the provenance reversal and owing nothing to
// that file.
//
// Every measured convention, the provenance table with file hashes, and the
// experiments behind each constant live in skullstrip_bench's
// test/skullstrip_reference_manifest.md. Read it before changing a constant here.
//
// Attribution: spatial normalisation adapted from AFNI (Robert W. Cox / NIMH,
// public domain), thd_brainormalize.c and thd_automask.c, AFNI rev 506e48403.

// The fixed human working grid. NOT derived from the input: for specie HUMAN,
// thd_brainormalize.c line 106 overwrites the min-voxel-size it just computed with
// THD_BN_DXYZ, so the grid is always 1 mm and the input resolution never reaches it.
// (AFNI's "SpatNorm resolution" debug line reports a pre-normalisation resample and
// is NOT this grid -- that misreading cost a round; see the manifest.)
#define SS_NX 167
#define SS_NY 212
#define SS_NZ 175
#define SS_DXYZ 1.0f
#define SS_XORG (-83.0f)
#define SS_YORG (-89.0f)
#define SS_ZORG (-82.0f)
#define SS_XCM 0.0f
#define SS_YCM 20.0f
// TRAP: thd_brainormalize.c guards ZCM with `#ifdef THD_BN_CMTOP`, and reading only
// that file suggests ZCM is 0. It is NOT -- thd_brainormalize.h line 1 is literally
// `#define THD_BN_CMTOP`, and the .c includes it at line 2, before the guard. So the
// shipped build takes the CMTOP branch: ZCM is 20, and the centre of mass is measured
// over only the top SS_CM_DEPTH mm rather than the whole volume. Getting this wrong
// puts the brain ~11 voxels too superior in the box and correlation against AFNI
// stalls near 0.47.
#define SS_ZCM 20.0f
#define SS_ZHEIGHT 170.0f
#define SS_CM_DEPTH 110.0f

typedef struct {
	// The normalised volume: SS_NX*SS_NY*SS_NZ bytes, 0..255, x fastest.
	// AFNI's working volume really is 8-bit -- mri_brainormalize returns MRI_byte
	// and 3dSkullStrip uses that directly (verified: the -write_spatnorm dataset
	// reports datum 'byte').
	unsigned char *vol;

	// Index warp, output grid -> RAI source grid: src = a*out + b, per axis.
	// Pure diagonal scale+shift; there is no rotation anywhere in this pipeline.
	float ai, bi, aj, bj, ak, bk;

	// The RAI permutation applied to the caller's image, in AFNI's signed-axis
	// encoding (+/-1,2,3 meaning source axis 1..3, negative = reversed).
	int fi, fj, fk;

	// Diagnostics, for the Milestone 1 stage gate against AFNI.
	float icm, jcm, kcm; // centre of mass, RAI source voxel indices
	int ktop, kbot;      // superior / inferior clip slices
	int support;         // mask voxel count after clipping
	int clip99;          // top-1% clip value
} ss_norm;

// Normalise nim (scalar 3D, any datatype niimath presents as float32) onto the fixed
// working grid. Returns 0 on success with *out populated and out->vol malloc'd;
// nonzero on failure with *out zeroed and nim untouched.
// in_datatype is the STORED NIfTI datatype before niimath's float32 promotion; see
// skullstrip_run. DT_NONE (0) means unknown and falls back to inspecting the values.
int ss_normalize(const nifti_image *nim, int in_datatype, ss_norm *out);

// Pull a working-grid volume back onto the caller's ORIGINAL grid and voxel order.
// This is the inverse of the index warp composed with the inverse of the RAI
// permutation; the caller's header is never touched, so qform/sform/pixdim are
// preserved by construction rather than by being rebuilt.
//
//   wvol   SS_NX*SS_NY*SS_NZ bytes on the working grid
//   dst    nim->nx*ny*nz floats, caller-allocated, fully overwritten
//   nearest  1 = nearest neighbour (use for masks: keeps labels exact)
//            0 = trilinear (use for intensities)
// Samples falling outside the working grid are written as 0.
// Returns 0 on success.
int ss_restore(const ss_norm *n, const unsigned char *wvol, const nifti_image *nim,
		float *dst, int nearest);

void ss_norm_free(ss_norm *n);

// ---------------------------------------------------------------------------
// Surface primitives. CLEAN-ROOM: ordinary computational geometry, written before the
// provenance reversal above and owing nothing to SUMA_BrainWrap.c. Validated against
// analytic answers in test_skullstrip_mesh.c.
// ---------------------------------------------------------------------------

typedef struct {
	int nv, nt;      // vertices, triangles
	float *v;        // 3*nv, xyz
	int *t;          // 3*nt, indices, consistently wound outward
	int *nbr;        // flattened neighbour lists
	int *nbr_off;    // nv+1 offsets into nbr
	float *nrm;      // 3*nv unit vertex normals (mean of UNIT face normals, NOT area-weighted)
} ss_mesh;

// Geodesic icosphere of frequency ld: 10*ld^2+2 vertices, 20*ld^2 triangles.
// AFNI's -ld/Icold is exactly this frequency -- measured, not assumed:
// `CreateIcosahedron -ld 20` emits 4002 nodes and 8000 triangles.
// Vertices are built by exact combinatorial identity (corner / edge / interior),
// never by welding coordinates, so topology is exact regardless of rounding.
int ss_icosphere(int ld, float radius, const float centre[3], ss_mesh *m);
void ss_mesh_free(ss_mesh *m);

// Recompute unit vertex normals from the current coordinates: each incident face normal is
// normalised FIRST, then averaged, so every face counts equally regardless of area (AFNI's
// SUMA_SurfNorm). Area weighting is measurably wrong here -- see ss_mesh_normals.
void ss_mesh_normals(ss_mesh *m);

// One neighbour-average smoothing pass. lambda in (0,1); reads a complete prior
// state and writes the next, so the result is independent of vertex visit order.
// scratch must hold 3*nv floats.
void ss_mesh_smooth(ss_mesh *m, float lambda, float *scratch);

// Count pairs of intersecting triangles, ignoring pairs that share a vertex.
// Uses a uniform spatial grid, so it is not the O(n^2) all-pairs test.
// Returns the count; negative on allocation failure.
long long ss_mesh_self_intersections(const ss_mesh *m);

// Rasterise a closed surface into a binary mask on a grid of nx*ny*nz unit voxels
// whose voxel (i,j,k) has its CENTRE at (i,j,k) in mesh coordinates. Parity fill
// along x with a half-open edge rule, so shared triangle edges are counted exactly
// once and the surface stays watertight. Returns 0 on success.
int ss_mesh_rasterize(const ss_mesh *m, int nx, int ny, int nz, unsigned char *mask);


// ---------------------------------------------------------------------------
// Deformation core. ADAPTED from AFNI's SUMA_BrainWrap.c (SUMA_StretchToFitLeCerveau,
// SUMA_LoadPrepInVol, SUMA_Find_IminImax), a non-copyrightable US Government work; the
// option defaults come from the public block in SUMA_3dSkullStrip.c. The underlying
// method is Smith 2002 (HBM 17:143-155), which AFNI's -help names as what it modifies.
// ---------------------------------------------------------------------------

typedef struct {
	float t2, t98, t, tm;   // percentiles and the median inside the brain sphere
	float cog[3];           // centre of gravity, working-grid voxel coords
	float radius;           // equivalent-sphere radius, mm
	long long nabove;       // voxels above t
} ss_stats;

// Defaults are AFNI's, read from the public-domain option block:
//   Zt (shrink_fac) 0.6, exp_frac 0.1, niter 250, d1 20 mm, d4 15 mm,
//   shrink_fac_bot_lim 0.65 with -no_use_edge (0.4 with edges), Icold/ld 20.
#define SS_SHRINK_FAC 0.6f
#define SS_SHRINK_BOT_NOEDGE 0.65f
#define SS_EXP_FRAC 0.1f
#define SS_NITER 250
#define SS_D1_MM 20.0f
#define SS_LD_NOEDGE 20
// AFNI -NNsmooth (72 geometric passes) and -max_inter_iter (4). On a self-intersecting
// surface AFNI discards it, adds SS_NNSMOOTH_STEP passes, and restarts from a fresh
// icosphere; smstep = 12*(Icold/25)^2 in INTEGER arithmetic is 0 at Icold 20, so the
// `if (smstep < 12) smstep = 12` floor is what actually applies -- the step is 12.
#define SS_NNSMOOTH 72
#define SS_NNSMOOTH_STEP 12
#define SS_MAX_INTER_ITER 4
#define SS_SMOOTH_END 20        // AFNI -smooth_final: Taubin passes after the repositions

// Compute the BET intensity statistics on the byte working volume.
// Returns 0 on success; nonzero if the contrast is degenerate.
int ss_intensity_stats(const unsigned char *vol, int nx, int ny, int nz, ss_stats *st);



// Touchup: push back out any vertex that stopped short of brain tissue, adapted from
// AFNI's SUMA_Reposition_Touchup + SUMA_Suggest_Touchup. A SINGLE pass: nodes move
// outward along their normals by at most limtouch mm, frozen nodes do not move at all,
// and in the inferior zone a node only moves if its neighbours want to move too.
// Returns the number of troubled nodes, or -1 on failure.
int ss_reposition_touchup(const unsigned char *vol, int nx, int ny, int nz,
		const ss_stats *st, ss_mesh *m, float *stop, float limtouch);

// Taubin lambda/mu smoothing, equal neighbour weights (AFNI -smooth_final, 20 passes
// at lambda 0.6307 / mu -0.6732). scratch must hold 6*nv floats.
void ss_taubin_smooth(ss_mesh *m, int niter, float lambda, float mu, float *scratch);

// AFNI's THD_mask_fillin_once (thd_automask.c): set a background voxel that has mask
// voxels on BOTH sides within nside steps along any ONE axis. This is -fill_hole's
// actual algorithm (default nside 10 when touchup is on) -- a directional gap filler,
// not a connected-component hole fill.
// Returns 0 on success, -1 on allocation failure -- it MUST be checked: silently skipping
// this step yields a plausible but materially different mask at exit 0.
int ss_mask_fillin_once(int nx, int ny, int nz, unsigned char *mmm, int nside);



// Build the final surface: deform, RETRY on self-intersection, then touch up.
//
// The retry decision uses ss_afni_self_intersect -- AFNI's own test, quirk included -- and
// NOT the correct ss_mesh_self_intersections; see the note at that function for why. On a
// hit the surface is DISCARDED, the smoothing rate rises by SS_NNSMOOTH_STEP, and expansion
// restarts from a fresh icosphere, up to max_retry times (max_retry <= 0 means
// SS_MAX_INTER_ITER, 4). The check sits on the surface straight out of the expansion loop,
// BEFORE any repositioning, matching AFNI: folds the later touchups introduce are accepted.
//
// If no attempt is clean the LAST attempt is returned. There is no best-of tracking, and the
// fold count across attempts is NOT monotone -- measured on T1w with smoothing forced off,
// successive attempts gave 47, 22, 36, 4 folds, so `max_retry` 2 returned a worse surface
// than `max_retry` 1. That matches AFNI, which also keeps whatever the last attempt produced,
// but it means more retries is not the same as a better surface.
// *n_inter reports the returned attempt's ss_mesh_self_intersections count -- the CORRECT
// test -- so *n_inter and the retry criterion are deliberately different numbers.
// *n_tries reports how many attempts ran.
// `fast` picks the deformation kernel: 0 is the faithful specialisation, bit-identical to the
// pre-optimisation release and the one the manifest's parity numbers are stated against; 1 is
// the default optimised one (1.8x, OpenMP-parallel over nodes). Same algorithm either way --
// see the note above the two #includes in skullstrip.c for why it is a build of the kernel
// rather than a runtime branch.
int ss_build_surface(const unsigned char *vol, int nx, int ny, int nz,
		const ss_stats *st, int ld, int niter, int max_retry,
		ss_mesh *m, long long *n_inter, int *n_tries, int fast);


// ---------------------------------------------------------------------------
// Top-level entry point (-skullstrip).
// ---------------------------------------------------------------------------
//
// Strips nim IN PLACE: normalise -> intensity statistics -> expand a surface ->
// touch up -> rasterise -> fill small holes -> pull the mask back onto the caller's
// grid -> apply it to the ORIGINAL intensities.
//
// Output contract (settled, and a deliberate divergence from AFNI's default, which
// rescales intensities): in-mask voxels keep their original values, out-of-mask
// voxels become the image MINIMUM -- a convention that was measured, not assumed
// rather than assumed (a -sub 1000 copy of T1w came back with range -1000..-315,
// so background sits at the minimum, not at zero).
//
// Fail-atomic: nim->data is only written once everything has succeeded, so a
// failure anywhere leaves the caller's image untouched.
// Requires a scalar 3D float32 image. Returns 0 on success.
// in_datatype is the image's STORED NIfTI datatype BEFORE niimath promoted it to float32
// (niimath's `ihdr.datatype`). It is not cosmetic: AFNI branches on the stored type when it
// converts to its internal short volume -- an integer source is copied verbatim while a float
// source is rescaled to fill the short range -- and that changes the histogram, hence every
// clip level downstream. Pass DT_NONE (0) if genuinely unknown; the heuristic below is then
// used as a fallback, with the caveat recorded at ss_to_rai_short.
// `fast` as for ss_build_surface: pass 0 for the bit-identical faithful kernel.
int skullstrip_run(nifti_image *nim, int in_datatype, int fast);

#ifdef __cplusplus
}
#endif

#endif // SKULLSTRIP_H
