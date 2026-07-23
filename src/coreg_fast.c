/* coreg_fast.c — fast multiresolution affine coregistration
 * (`-cost fast`/`-cost fastx`/`-cost fasthel`/`-cost fastcr`).
 *
 * Independent implementation of SPM/FLIRT-inspired ideas (clean-room: BSD niimath coreFLT.c ports, see AGENTS.md;
 * clean-room: never derived from src/GPL/). Owns its parameterization, cost
 * evaluators, level builder, and optimizer orchestration. Reuses only the pure BSD
 * demo helpers: al_image_xform() (sform/qform policy), nii_reslice_affine() (final
 * warp), powell_newuoa() (optimizer), mat44 math, and nifti_smooth_gauss_f32()
 * (the float32 core's raw-buffer Gaussian blur).
 *
 * Serial-only at the host-call level. The NEWUOA callback context is thread-local so
 * `fastx` can run its independent coarse candidates concurrently, but two host calls
 * must not enter the estimator at the same time.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <stdint.h>
#include <limits.h>
#include <time.h>

#include "coreg_fast.h"
#include "core32.h"         /* nifti_smooth_gauss_f32 */
#include "allineate.h"      /* al_image_xform, nii_reslice_affine, AL_INTERP_* */
#include "al_thread_local.h"
#include "al_size_guard.h"
#include "coreg_fast_nms.h"

extern int powell_newuoa(int ndim, double *x, double rstart, double rend,
                         int maxcall, double (*ufunc)(int, double *));
extern void powell_newuoa_free_threadlocal(void);
extern void powell_set_mfac(float mm, float aa);
extern void powell_get_mfac(float *mm, float *aa);

#if defined(_OPENMP)
#include <omp.h>
#endif

/* Below which sample count the CR evaluator stays serial (parallel setup/reduction
 * overhead exceeds the benefit at the coarse levels; the 2 mm level dominates). */
#define CF_CR_PAR_MIN 30000
/* The CR reduction uses a FIXED number of contiguous sample chunks, independent of
 * the thread count: each chunk is summed in order into its own accumulator and the
 * chunks are combined in a fixed order, so the double-precision result is identical
 * for any thread count (the NEWUOA path — and thus the fit — is bit-identical at
 * -p 1 and -p N). Threads are mapped onto these chunks; which thread runs which
 * chunk does not affect the result. Must be >= the max useful worker count. */
#define CF_CR_NCHUNK 64
#define CF_HEL_NBIN  32   /* bins per axis for the CF_COST_HEL 2D joint histogram */

/* HEL joint-foreground policy: a fixed sample whose corresponding MOVING intensity is at
 * background/air (<= this fraction of the moving level's dynamic range above its minimum)
 * carries no cross-modal information and is dropped from the joint histogram. The fixed side
 * is already foreground-thresholded (cf_make_samples); this is its moving-side partner, so the
 * statistics use only the INTERSECTION of both foregrounds. It deepens the otherwise-shallow
 * HEL well produced by a full-head moving whose air/sinuses are hard-zeroed (masked scans) or
 * by a skull-stripped/hard-zeroed base (e.g. an @SSwarper SSW template) — those zero-mapped
 * pairs pile into moving-bin-0 as an ~independence floor that flattens the landscape and lets
 * the sparse coarse search lock a wrong pitch. Validated across the full benchmark modality
 * sweep (T1/T2/FLAIR/EPI x plain + SSW bases): large gains on hard-zeroed bases, pose unchanged
 * (<1 deg) elsewhere. The threshold must stay SMALL (only true background) — >=2% re-eats real
 * low-signal cortex and can flip a fit. Deterministic per-sample test -> -p1==-pN stays
 * bit-identical. HEL only (the CR/LS reductions are not wired for it).
 *
 * GATED ON A HARD-ZEROED BASE (cf_base_is_hard_zeroed) — do not apply it unconditionally.
 * The drop is NOT free: `nin` (the overlap floor) still counts a dropped sample, and only the
 * weighted branch has a kept-mass floor, so an unweighted pose pays nothing for discarding base
 * foreground onto moving air. That is safe only when the discard is SURGICAL, which is exactly
 * what a masked/skull-stripped base guarantees: its 5%-of-range foreground (cf_make_samples) is
 * compact real tissue, so at any sane pose it lands inside the moving head and <10% drops.
 * A whole-head base is the opposite — MNI152_T1_1mm's foreground is 68% of its volume (a broad
 * low-intensity halo of scalp/skull/neck) and legitimately extends past a subject's imaged
 * extent, so on a head that does not fill its FOV the drop silently deletes ~38% of the
 * constraint AT EVERY POSE. That was a real, severe regression: M2017 (head = 21% of a
 * 160x256x256 FOV) converged on a pose retaining only 52% of base foreground, masked-Hellinger
 * 0.1814 -> 0.1233. It is the FLIRT Fig-1 / partial-FOV degeneracy this engine documents
 * elsewhere: in-FOV air is OBSERVED evidence ("no tissue here"), not missing data, and treating
 * it as missing lets a pose choose its own denominator. Gating restores every whole-head-base
 * fit BYTE-FOR-BYTE to its pre-drop result while keeping the SSW gains bit-identical (validated:
 * 14/14 whole-head cases byte-identical, 6/6 SSW cases unchanged). Do NOT ungate it; a smooth
 * "charge for the discard" blend was tried and fails because the discarded fraction is nearly
 * pose-invariant (0.59->0.64 across the whole coarse sweep), so it is only a constant offset. */
#define CF_HEL_MDARK_FRAC 0.01

/* Fraction of the base that must sit exactly AT its minimum for the base to count as
 * masked/skull-stripped. Well separated in practice: MNI152_T1_1mm 0.061, avg152T1 0.0001,
 * MNI152_2009_SSW 0.622 — 0.25 sits an order of magnitude clear of both sides. */
#define CF_HEL_MDARK_BASEFRAC 0.25

/* Per-stage profiling, enabled by -DAL_PROFILE (same flag as allineate.c / `make
 * profile`). Output goes to stderr as " [coreg profile] ...". The timing points are
 * per-level (a handful of calls), never per cost evaluation, so the release build is
 * untouched and even the profiling build adds negligible overhead. The final warp is
 * applied by main.c (nii_apply_affine), so it is NOT part of this estimate profile. */
/* Wall-clock timer in seconds (OpenMP-aware; mirrors allineate.c's al_wtime, including the
   Windows fallback so MSVC — which lacks clock_gettime(CLOCK_MONOTONIC) — still links). */
static inline double cf_wtime(void) {
#ifdef _OPENMP
    return omp_get_wtime();
#elif defined(_WIN32)
    return (double)clock() / CLOCKS_PER_SEC; /* fallback: CPU time on Windows without OpenMP */
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
#endif
}

#ifdef AL_PROFILE
static double g_prof_blur = 0.0, g_prof_resample = 0.0;
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define CF_DEG2RAD (M_PI / 180.0)

/* Normalized-parameter -> physical scale (1 normalized unit ~ 1 mm displacement at a
 * 100 mm radius, so translation/rotation/scale/shear steps are comparable to NEWUOA). */
#define CF_PS_TRANS 1.0    /* mm per unit */
#define CF_PS_ROT   0.01   /* rad per unit  (~0.573 deg) */
#define CF_PS_SCALE 0.01   /* fractional scale per unit (1%) */
#define CF_PS_SHEAR 0.01   /* shear per unit */

/* Coarse multi-start (FLIRT-inspired; Jenkinson & Smith 2001 -- from the PUBLISHED
 * description only, never AFNI/GPL source). The 8 mm level fits RIGID, which is right for a
 * whole-head base but CANNOT represent the true pose when the base is skull-stripped: with
 * scale locked at 1 the CORRECT orientation genuinely scores WORSE than a wrong rotated one
 * (measured, M2017 -> MNI152_2009_SSW at 8 mm: true orientation rigid 0.9433 vs the spurious
 * 30-deg-pitch pose 0.9241), so no amount of searching at that level can find it. Conversely,
 * freeing coarse scale unconditionally re-opens the spurious-shrink basin the rigid-coarse
 * rule was introduced to close (measured: T1w1mm -> a stripped base collapses, NCC 0.48 -> 0.10).
 *
 * No single coarse rule is right for every pair, so run all three as independent starts through
 * the whole pyramid and keep whichever ends with the lower FINE-level cost. That is exactly the
 * multi-start philosophy FLIRT uses to escape local minima, and it is safe here because the
 * 2 mm cost is a reliable arbiter even when the 8 mm cost is not (verified on the four decisive
 * cases: the lower final cost picked the better NCC every time, both directions).
 * Pass 0 == the historical rigid-coarse trajectory, so a pair that pass 0 already wins is
 * bit-identical to before. */
#define CF_NSTRAT 3   /* 0: rigid coarse (historical)  1: scale-bracketed coarse  2: CR coarse */
/* Pass 1 coarse scale bracket. Isotropic factors x an extra z-ratio: a skull-stripped base
 * needs a markedly ANISOTROPIC fit (M2017 -> SSW: 0.93/0.90/0.81), and an isotropic-only
 * bracket leaves the true pose tied with the spurious one (0.9224 vs 0.9241 -- too thin to
 * survive refinement), whereas allowing the z ratio wins by a decisive margin (0.8937). */
#define CF_CSC_N 3
#define CF_CZR_N 2
static const double CF_ANG[5] = {-30,-15,0,15,30};

/* -weight uses AFNI 3dAllineate's scheme: a base(fixed)-space GRADED weight, normalized to
 * [0,1] over the whole fixed grid (divide by its max), applied per fixed sample. It is NOT a
 * soft-focus floor — a voxel weighted 0 is excluded and one weighted near 1 dominates, exactly
 * as AFNI weights each base voxel. A whole-head weight that keeps the scalp attenuated (nonzero)
 * still anchors global scale; a weight that fully zeroes everything outside the brain removes the
 * scale anchor (the FLIRT Fig-1 shrink risk), so — like AFNI — supply an attenuated (not zeroed)
 * out-of-ROI weight. HEL/CR are scale-invariant, so the fit is invariant to a global weight
 * scaling. See the weight build below. */

/* Guardrails (normalized units), wider than the supported capture envelope so a
 * valid fit is never clipped; a fit landing on a guardrail is reported as failure.
 * Translation needs more room than the nominal capture range: an oblique/reoriented
 * but substantially overlapping head can require an 80+ mm COM seed. Rejecting that
 * seed at 64 mm forces the search to mimic translation with a wrong rotation. */
static const double CF_LO[12] = {-128,-128,-128, -80,-80,-80,  -45,-45,-45,  -22,-22,-22 };
static const double CF_HI[12] = { 128, 128, 128,  80, 80, 80,   55, 55, 55,   22, 22, 22 };

/*==========================================================================*/
/* mat44 helpers                                                             */
/*==========================================================================*/

static mat44 cf_ident(void) {
    mat44 m; memset(&m, 0, sizeof m);
    m.m[0][0] = m.m[1][1] = m.m[2][2] = m.m[3][3] = 1.0f;
    return m;
}
static void cf_apply(const mat44 *M, double x, double y, double z,
                     double *ox, double *oy, double *oz) {
    *ox = M->m[0][0]*x + M->m[0][1]*y + M->m[0][2]*z + M->m[0][3];
    *oy = M->m[1][0]*x + M->m[1][1]*y + M->m[1][2]*z + M->m[1][3];
    *oz = M->m[2][0]*x + M->m[2][1]*y + M->m[2][2]*z + M->m[2][3];
}
static double cf_colnorm(const mat44 *M, int c) {
    return sqrt((double)M->m[0][c]*M->m[0][c] + (double)M->m[1][c]*M->m[1][c] +
                (double)M->m[2][c]*M->m[2][c]);
}
static int cf_finite(double v) { return v >= -DBL_MAX && v <= DBL_MAX; }
static int cf_affine_finite(mat44 m) {
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 4; j++)
            if (!cf_finite(m.m[i][j])) return 0;
    return 1;
}

/* Validate a single-volume 3D input once on entry (matches al_dims_ok semantics:
 * reject 4D/multi-volume, bound dims to int range, require nvox == nx*ny*nz with an
 * overflow-safe product). Lower layers then trust the shape. Returns 0 if OK. */
static int cf_dims_ok(const nifti_image *n, const char *who) {
    if (n->nz < 2 || n->ny < 1 || n->nx < 1) {
        fprintf(stderr, "coreg fast: %s must be a 3D volume (nz>=2)\n", who); return 1; }
    if (n->nt > 1 || n->nu > 1 || n->nv > 1 || n->nw > 1) {
        fprintf(stderr, "coreg fast: %s 4D/multi-volume input not supported (3D only)\n", who); return 1; }
    if (n->nx > INT_MAX || n->ny > INT_MAX || n->nz > INT_MAX) {
        fprintf(stderr, "coreg fast: %s dimension exceeds int range\n", who); return 1; }
    size_t a = (size_t)n->nx, b = (size_t)n->ny, c = (size_t)n->nz;
    if ((b && a > SIZE_MAX / b) || (c && a*b > SIZE_MAX / c)) {
        fprintf(stderr, "coreg fast: %s voxel count overflow\n", who); return 1; }
    if ((size_t)n->nvox != a*b*c) {
        fprintf(stderr, "coreg fast: %s nvox does not match nx*ny*nz (4D/multi-volume?)\n", who); return 1; }
    /* The engine and the blur/transpose helpers use `int` for the voxel count, sample
       count, and row/transpose offsets, so the spatial product must fit `int` — bound
       it once here (matches al_safe_nvox), not scattered through inner loops. */
    size_t nvox = a*b*c;
    if (nvox > (size_t)INT_MAX || !al_float_nvox_fits((uint64_t)nvox)) {
        fprintf(stderr, "coreg fast: %s voxel count exceeds supported range\n", who); return 1; }
    return 0;
}

/*==========================================================================*/
/* float extraction (does not mutate the input image)                        */
/*==========================================================================*/

/* Extract nim->data into a fresh float buffer. Scaling matches the allineate
 * nii_to_float() convention exactly: scl_slope/scl_inter are applied ONLY when
 * scl_slope != 0 and not the identity (slope 1, inter 0). A zero slope means "no
 * scaling" per the NIfTI spec (inter is NOT added) — so the fast path sees the same
 * intensities (and the same COM seed) as the allineate engine for any valid NIfTI.
 * Returns malloc'd buffer (caller frees) or NULL. Non-finite -> 0. */
static float *cf_extract_float(const nifti_image *nim, size_t *nvox_out) {
    size_t nv = (size_t)nim->nx * nim->ny * nim->nz;
    if (nv == 0 || nim->data == NULL) return NULL;
    float *out = (float *)malloc(nv * sizeof(float));
    if (!out) return NULL;
    int do_scale = (nim->scl_slope != 0.0f &&
                    !(nim->scl_slope == 1.0f && nim->scl_inter == 0.0f));
    double sl = do_scale ? (double)nim->scl_slope : 1.0;
    double in = do_scale ? (double)nim->scl_inter : 0.0;
    const void *d = nim->data;
    for (size_t i = 0; i < nv; i++) {
        double v;
        switch (nim->datatype) {
            case DT_UINT8:   v = ((const uint8_t  *)d)[i]; break;
            case DT_INT8:    v = ((const int8_t   *)d)[i]; break;
            case DT_INT16:   v = ((const int16_t  *)d)[i]; break;
            case DT_UINT16:  v = ((const uint16_t *)d)[i]; break;
            case DT_INT32:   v = ((const int32_t  *)d)[i]; break;
            case DT_UINT32:  v = ((const uint32_t *)d)[i]; break;
            case DT_INT64:   v = (double)((const int64_t  *)d)[i]; break;
            case DT_UINT64:  v = (double)((const uint64_t *)d)[i]; break;
            case DT_FLOAT32: v = ((const float    *)d)[i]; break;
            case DT_FLOAT64: v = ((const double   *)d)[i]; break;
            default: free(out); return NULL;
        }
        v = v * sl + in;
        out[i] = (v >= -FLT_MAX && v <= FLT_MAX) ? (float)v : 0.0f;
    }
    if (nvox_out) *nvox_out = nv;
    return out;
}

/*==========================================================================*/
/* pyramid level (built on demand, released after its stage)                 */
/*==========================================================================*/

typedef struct { float *data; int nx, ny, nz; mat44 v2w; mat44 w2v; } cf_level;

static void cf_free_level(cf_level *L) { if (L) { free(L->data); L->data = NULL; } }

/* Trilinear sample of a level buffer at fractional voxel coord (fx,fy,fz).
 * Returns 0 (and *ok=0) when outside the field of view. */
static inline float cf_trilerp(const cf_level *L, double fx, double fy, double fz, int *ok) {
    if (fx < -0.5 || fy < -0.5 || fz < -0.5 ||
        fx > L->nx - 0.5 || fy > L->ny - 0.5 || fz > L->nz - 0.5) { *ok = 0; return 0.0f; }
    int x0 = (int)floor(fx), y0 = (int)floor(fy), z0 = (int)floor(fz);
    double dx = fx - x0, dy = fy - y0, dz = fz - z0;
    int x1 = x0 + 1, y1 = y0 + 1, z1 = z0 + 1;
    if (x0 < 0) x0 = 0; if (y0 < 0) y0 = 0; if (z0 < 0) z0 = 0;
    if (x1 > L->nx - 1) x1 = L->nx - 1; if (y1 > L->ny - 1) y1 = L->ny - 1;
    if (z1 > L->nz - 1) z1 = L->nz - 1;
    if (x0 > L->nx - 1) x0 = L->nx - 1; if (y0 > L->ny - 1) y0 = L->ny - 1;
    if (z0 > L->nz - 1) z0 = L->nz - 1;
    int nxy = L->nx * L->ny;
    const float *d = L->data;
    #define V(a,b,c) d[(a) + (b)*L->nx + (size_t)(c)*nxy]
    double c00 = V(x0,y0,z0)*(1-dx) + V(x1,y0,z0)*dx;
    double c10 = V(x0,y1,z0)*(1-dx) + V(x1,y1,z0)*dx;
    double c01 = V(x0,y0,z1)*(1-dx) + V(x1,y0,z1)*dx;
    double c11 = V(x0,y1,z1)*(1-dx) + V(x1,y1,z1)*dx;
    #undef V
    double c0 = c00*(1-dy) + c10*dy, c1 = c01*(1-dy) + c11*dy;
    *ok = 1;
    return (float)(c0*(1-dz) + c1*dz);
}

/* Build an approximately-isotropic level of spacing sep_mm from a full-resolution
 * float source (buffer `src` with dims/pixdim and its own vox->world S). The source
 * is blurred (Gaussian sigma from smooth_fwhm_mm) into a scratch copy, then resampled
 * onto the level grid. The level grid shares the source's orientation and world
 * centre; its FOV covers the source. Returns 0 on success. */
static int cf_build_level(const float *src, int snx, int sny, int snz,
                          const mat44 *S, double sep_mm, double smooth_fwhm_mm,
                          cf_level *out) {
    double vx = cf_colnorm(S, 0), vy = cf_colnorm(S, 1), vz = cf_colnorm(S, 2);
    if (vx <= 0 || vy <= 0 || vz <= 0) return 1;
    /* Level dims must cover the source FOV: round UP so the level field spans the
     * full source field (plain rounding could drop up to ~half a voxel of foreground
     * at the edges). Bound the world-FOV BEFORE the int cast — cf_dims_ok bounds only
     * the source voxel count, not this derived level grid, so a finite-but-extreme
     * sform/pixdim could otherwise overflow the (int) conversion (UB) and the int
     * index math (nx*ny) downstream. */
    double dLx = ceil(snx * vx / sep_mm), dLy = ceil(sny * vy / sep_mm), dLz = ceil(snz * vz / sep_mm);
    if (dLx < 1) dLx = 1; if (dLy < 1) dLy = 1; if (dLz < 1) dLz = 1;
    if (dLx > INT_MAX || dLy > INT_MAX || dLz > INT_MAX ||
        dLx*dLy*dLz > (double)INT_MAX ||
        dLx*dLy*dLz > (double)(AL_SIZE_MAX / sizeof(float)))
        return 1;
    /* Resource bound (not just index): a pyramid level must never GROSSLY up-sample the source
       (that only happens for out-of-envelope inputs — a source whose pixdim is far coarser than
       sep_mm — and a level < INT_MAX voxels can still be ~8 GB). Reject a level exceeding 8x the
       source voxel count; a real anisotropic->isotropic level stays below the source count. */
    if (dLx*dLy*dLz > 8.0 * (double)snx*sny*snz) return 1;
    int Lx = (int)dLx, Ly = (int)dLy, Lz = (int)dLz;
    /* level vox->world: same direction cosines as S, isotropic spacing sep_mm */
    mat44 v2w = cf_ident();
    for (int r = 0; r < 3; r++) {
        v2w.m[r][0] = (float)(S->m[r][0] / vx * sep_mm);
        v2w.m[r][1] = (float)(S->m[r][1] / vy * sep_mm);
        v2w.m[r][2] = (float)(S->m[r][2] / vz * sep_mm);
    }
    /* place level centre at source world centre */
    double scx, scy, scz;
    cf_apply(S, (snx-1)*0.5, (sny-1)*0.5, (snz-1)*0.5, &scx, &scy, &scz);
    double lcx, lcy, lcz;
    cf_apply(&v2w, (Lx-1)*0.5, (Ly-1)*0.5, (Lz-1)*0.5, &lcx, &lcy, &lcz);
    v2w.m[0][3] += (float)(scx - lcx);
    v2w.m[1][3] += (float)(scy - lcy);
    v2w.m[2][3] += (float)(scz - lcz);

    /* blur a scratch copy of the source in its own frame */
    size_t snv = (size_t)snx * sny * snz;
    float *blur = (float *)malloc(snv * sizeof(float));
    if (!blur) return 1;
    memcpy(blur, src, snv * sizeof(float));
    double sig = smooth_fwhm_mm / 2.354820045; /* FWHM -> sigma (mm) */
#ifdef AL_PROFILE
    double _tb = cf_wtime();
#endif
    if (sig > 0.0)
        if (nifti_smooth_gauss_f32(blur, snx, sny, snz, 1,
                                   (float)vx, (float)vy, (float)vz,
                                   (float)sig, (float)sig, (float)sig, -6.0f)) {
            free(blur);
            return 1;
        }
#ifdef AL_PROFILE
    g_prof_blur += cf_wtime() - _tb;
    double _tr = cf_wtime();
#endif

    size_t lnv = (size_t)Lx * Ly * Lz;
    float *ld = (float *)malloc(lnv * sizeof(float));
    if (!ld) { free(blur); return 1; }
    mat44 Sinv = nifti_mat44_inverse(*S);
    mat44 gam  = nifti_mat44_mul(Sinv, v2w);   /* level index -> source index */
    cf_level tmp = { blur, snx, sny, snz, *S, Sinv };
    for (int z = 0; z < Lz; z++)
        for (int y = 0; y < Ly; y++)
            for (int x = 0; x < Lx; x++) {
                double fx, fy, fz; cf_apply(&gam, x, y, z, &fx, &fy, &fz);
                int ok; float v = cf_trilerp(&tmp, fx, fy, fz, &ok);
                ld[x + (size_t)y*Lx + (size_t)z*Lx*Ly] = ok ? v : 0.0f;
            }
    free(blur);
#ifdef AL_PROFILE
    g_prof_resample += cf_wtime() - _tr;
#endif
    out->data = ld; out->nx = Lx; out->ny = Ly; out->nz = Lz;
    out->v2w = v2w; out->w2v = nifti_mat44_inverse(v2w);
    return 0;
}

/*==========================================================================*/
/* affine parameterization (pure function of params + world centre)          */
/*==========================================================================*/

typedef struct {
    /* current stage context (read by the NEWUOA callback via g_cf) */
    const cf_level *Lf, *Lm;   /* fixed & moving levels */
    mat44 Sf_lvl;              /* == Lf->v2w (cached; the cost maps fixed-level index) */
    double cx, cy, cz;         /* world rotation centre (fixed FOV centre) */
    double base[12];           /* frozen params for non-free dims (normalized) */
    int    free_idx[12], nfree;
    int    global_scale;       /* if 1, free dim mapping to idx6 drives sx=sy=sz */
    int    cost;               /* CF_COST_* */
    int    serial_cost;        /* suppress inner cost OpenMP while fastx parallelizes candidates */
    /* precomputed fixed sample set for the current level */
    int    ns; int32_t *sx, *sy, *sz; float *fval; int *fbin; int K;
    const cf_level *Lw;        /* optional fixed-space weight level for THIS stage (shares Lf's
                                  grid); NULL at the coarse level and when no -weight is given */
    float *fwt;                /* per-sample weight (Lw sampled at each fixed voxel), or NULL for
                                  unweighted — NULL makes every cost multiply-by-1.0 (bit-identical) */
    double fwt_total;          /* sum of fwt over all samples (total weighted mass); used for the
                                  weighted in-FOV support floor. 0 when unweighted. */
    double fvar, fmean;        /* fixed-sample variance/mean (LS) */
    double m_bg;               /* moving min: LS out-of-FOV fill + HEL binning origin */
    double m_top;              /* moving-level max (for CF_COST_HEL moving-axis binning) */
    double thr_frac;           /* foreground threshold (fraction of dynamic range), set once */
    int    hel_mdark;          /* 1: apply the HEL joint-foreground drop. Set once from the
                                  full-res base (cf_base_is_hard_zeroed); 0 (memset default)
                                  disables it, which is the fail-safe whole-head behavior. */
    int    opt_err;            /* set if any NEWUOA call returned a negative error code */
    int    dof_run;            /* highest DOF actually fitted (tracks real execution, so a
                                  shortened debug schedule serializes the true DOF) */
    double *cr_acc;            /* CR parallel scratch: CF_CR_NCHUNK blocks of (3K+4)
                                  doubles (per-chunk bn/bs/bq histograms + N/Sy/Syy/nin),
                                  reused across evals; combined in fixed chunk order for a
                                  thread-count-independent (deterministic) result */
    int    evals;
} cf_ctx;

static AL_THREAD_LOCAL cf_ctx *g_cf = NULL;

#ifdef COREG_FAST_TEST_ALLOC
/* Test-only fault injection: fail the second sample-array allocation on the
 * requested cf_make_samples call. This exercises a partial allocation after an
 * earlier multi-start strategy has completed without changing release behavior. */
static int g_cf_test_fail_sample_call = 0;
static int g_cf_test_sample_calls = 0;
void coreg_fast_test_fail_samples_on_call(int call) {
    g_cf_test_fail_sample_call = call;
    g_cf_test_sample_calls = 0;
}
int coreg_fast_test_sample_call_count(void) { return g_cf_test_sample_calls; }
#endif

/* Build the world-mm FIXED->MOVING affine from normalized params about centre c. */
static mat44 cf_affine(const cf_ctx *c, const double p[12]) {
    double tx = p[0]*CF_PS_TRANS, ty = p[1]*CF_PS_TRANS, tz = p[2]*CF_PS_TRANS;
    double rx = p[3]*CF_PS_ROT,   ry = p[4]*CF_PS_ROT,   rz = p[5]*CF_PS_ROT;
    double sx = 1.0 + p[6]*CF_PS_SCALE, sy = 1.0 + p[7]*CF_PS_SCALE, sz = 1.0 + p[8]*CF_PS_SCALE;
    double hxy = p[9]*CF_PS_SHEAR, hxz = p[10]*CF_PS_SHEAR, hyz = p[11]*CF_PS_SHEAR;
    double cxr = cos(rx), sxr = sin(rx), cyr = cos(ry), syr = sin(ry), czr = cos(rz), szr = sin(rz);
    /* R = Rx*Ry*Rz */
    double R[3][3];
    double Rz[3][3] = {{czr,-szr,0},{szr,czr,0},{0,0,1}};
    double Ry[3][3] = {{cyr,0,syr},{0,1,0},{-syr,0,cyr}};
    double Rx[3][3] = {{1,0,0},{0,cxr,-sxr},{0,sxr,cxr}};
    double RyRz[3][3];
    for (int i=0;i<3;i++) for (int j=0;j<3;j++){ double s=0; for(int k=0;k<3;k++) s+=Ry[i][k]*Rz[k][j]; RyRz[i][j]=s; }
    for (int i=0;i<3;i++) for (int j=0;j<3;j++){ double s=0; for(int k=0;k<3;k++) s+=Rx[i][k]*RyRz[k][j]; R[i][j]=s; }
    /* Scale*Shear: Sh upper-triangular unit-diagonal, then diag scale */
    double Sh[3][3] = {{1,hxy,hxz},{0,1,hyz},{0,0,1}};
    double SS[3][3] = {{sx,0,0},{0,sy,0},{0,0,sz}};
    double SShr[3][3];
    for (int i=0;i<3;i++) for (int j=0;j<3;j++){ double s=0; for(int k=0;k<3;k++) s+=SS[i][k]*Sh[k][j]; SShr[i][j]=s; }
    double L[3][3];
    for (int i=0;i<3;i++) for (int j=0;j<3;j++){ double s=0; for(int k=0;k<3;k++) s+=R[i][k]*SShr[k][j]; L[i][j]=s; }
    mat44 M = cf_ident();
    for (int i=0;i<3;i++) for (int j=0;j<3;j++) M.m[i][j]=(float)L[i][j];
    /* translation: moving = c + t + L*(fixed - c) */
    double cc[3] = { c->cx, c->cy, c->cz };
    double t[3]  = { tx, ty, tz };
    for (int i=0;i<3;i++){ double Lc=0; for(int k=0;k<3;k++) Lc+=L[i][k]*cc[k]; M.m[i][3]=(float)(cc[i]+t[i]-Lc); }
    return M;
}

/* Expand the NEWUOA free vector x[nfree] into a full 12-param vector. */
static void cf_expand(const cf_ctx *c, const double *x, double p[12]) {
    memcpy(p, c->base, sizeof(double)*12);
    for (int i = 0; i < c->nfree; i++) p[c->free_idx[i]] = x[i];
    if (c->global_scale) { p[7] = p[6]; p[8] = p[6]; }
}

/*==========================================================================*/
/* cost evaluators                                                           */
/*==========================================================================*/

#define CF_PENALTY 1.0e9


static double cf_cost_eval(const double p[12]) {
    cf_ctx *c = g_cf;
    c->evals++;   /* count EVERY cost evaluation (coarse grid, seeds, refine, final) */
    /* bounds guardrail */
    for (int i = 0; i < 12; i++) if (p[i] < CF_LO[i] || p[i] > CF_HI[i]) return CF_PENALTY;
    mat44 F2M = cf_affine(c, p);
    if (!cf_affine_finite(F2M)) return CF_PENALTY;
    /* fixed-level index -> moving-level index: w2v_m * F2M * v2w_f */
    mat44 gam = nifti_mat44_mul(c->Lm->w2v, nifti_mat44_mul(F2M, c->Sf_lvl));
    int ns = c->ns;
    const cf_level *Lm = c->Lm;

    /* LS retains the fixed-background-fill policy below. The two public fast costs
       (HEL/CR) instead form their statistics from the actual image intersection:
       treating an acquisition boundary as moving background biases partial-FOV fits
       toward scale/shear transforms that pull every fixed sample inside the moving
       box. A minimum in-FOV floor still rejects a degenerate near-zero intersection. */
    double m_bg = c->m_bg;

    if (c->cost == CF_COST_LS) {
        double Sm=0, Smm=0, Sfm=0; long nin=0;
        double fmean = c->fmean;
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+:Sm,Smm,Sfm,nin) if(!c->serial_cost && ns>20000)
        #endif
        for (int s = 0; s < ns; s++) {
            double ix, iy, iz; cf_apply(&gam, c->sx[s], c->sy[s], c->sz[s], &ix, &iy, &iz);
            int ok; float v = cf_trilerp(Lm, ix, iy, iz, &ok);
            double mv = ok ? v : m_bg;
            double fv = c->fval[s] - fmean;
            Sm += mv; Smm += (double)mv*mv; Sfm += fv*mv; if (ok) nin++;
        }
        if (nin < (long)(0.10*ns) || nin < 16) return CF_PENALTY;
        double invn = 1.0/ns;
        double mm = Sm*invn;
        double varM = Smm*invn - mm*mm;
        double cov  = Sfm*invn;                /* fixed already mean-centred over full set */
        if (varM <= 1e-12 || c->fvar <= 1e-12) return CF_PENALTY;
        double r = cov / sqrt(c->fvar * varM);
        if (r < -1) r = -1; if (r > 1) r = 1;
        return 1.0 - r;   /* minimize */
    }

    if (c->cost == CF_COST_HEL) {
        /* CF_COST_HEL: Hellinger-affinity (mutual-information family) — minimize the Bhattacharyya
           coefficient BC = sum_ij sqrt(P(i,j)*Pf(i)*Pm(j)) between the joint (fixed,moving) 2D
           histogram and the product of its marginals. BC->1 when fixed/moving are independent
           (misaligned) and is minimized when they are maximally dependent (aligned). Far more
           robust than CR for cross-modal. Clean-room from the standard definition (NOT from AFNI
           / GPL). NB x NB bins (NB = CF_HEL_NBIN): fixed reuses fbin (mapped to NB), moving binned over
           [m_bg, m_top]. Deterministic fixed-chunk reduction like CR -> bit-identical -p1==-pN. */
        const int NB = CF_HEL_NBIN;
        const int NC = CF_CR_NCHUNK;
        int hstride = NB*NB + 1;                 /* joint[NB*NB] + nin */
        double *acc = c->cr_acc;
        memset(acc, 0, (size_t)NC * hstride * sizeof(double));
        double mspan = c->m_top - c->m_bg; if (!(mspan > 0.0)) mspan = 1.0;
        double minv = (NB - 1) / mspan;
        /* Joint-foreground cutoff: drop pairs whose moving value is background/air (m_bg + 1% of
           range). Local like mspan/minv above — a HEL-only derived constant, not a shared ctx field.
           Only for a masked/skull-stripped base; -DBL_MAX makes the per-sample test below always
           false (never drops), so a whole-head base is bit-identical to no cutoff at all. */
        double m_dark_thr = c->hel_mdark ? c->m_bg + CF_HEL_MDARK_FRAC * (c->m_top - c->m_bg)
                                         : -DBL_MAX;
        int Kf = c->K;
        const float *fwt = c->fwt;   /* NULL -> weight 1.0 (default, bit-identical) */
#ifdef _OPENMP
        int par = (!c->serial_cost && ns > CF_CR_PAR_MIN);
        #pragma omp parallel for schedule(static) if(par)
#endif
        for (int cix = 0; cix < NC; cix++) {
            int lo = (int)((int64_t)cix * ns / NC);
            int hi = (int)((int64_t)(cix + 1) * ns / NC);
            double *H = acc + (size_t)cix * hstride;
            long lnin = 0;
            for (int s = lo; s < hi; s++) {
                double ix, iy, iz; cf_apply(&gam, c->sx[s], c->sy[s], c->sz[s], &ix, &iy, &iz);
                int ok; float v = cf_trilerp(Lm, ix, iy, iz, &ok);
                if (!ok) continue;
                double mv = v;
                lnin++;   /* raw in-FOV count for the overlap floor (weight-independent) */
                if (mv <= m_dark_thr) continue;  /* joint-foreground: drop moving-background pairs (see CF_HEL_MDARK_FRAC) */
                int fb = (int)((long)c->fbin[s] * NB / Kf); if (fb < 0) fb = 0; if (fb >= NB) fb = NB-1;
                int mb = (int)((mv - c->m_bg) * minv + 0.5); if (mb < 0) mb = 0; if (mb >= NB) mb = NB-1;
                H[fb*NB + mb] += fwt ? fwt[s] : 1.0;
            }
            H[NB*NB] = (double)lnin;
        }
        /* combine chunk joints into block 0 in fixed order (deterministic) */
        long nin = (long)acc[NB*NB];
        for (int cix = 1; cix < NC; cix++) {
            const double *Ht = acc + (size_t)cix * hstride;
            for (int i = 0; i < NB*NB; i++) acc[i] += Ht[i];
            nin += (long)Ht[NB*NB];
        }
        if (nin < (long)(0.10*ns) || nin < 16) return CF_PENALTY;
        double Pf[CF_HEL_NBIN], Pm[CF_HEL_NBIN], tot = 0.0;
        for (int i = 0; i < NB; i++) { Pf[i] = 0.0; Pm[i] = 0.0; }
        for (int i = 0; i < NB; i++)
            for (int j = 0; j < NB; j++) { double h = acc[i*NB+j]; Pf[i] += h; Pm[j] += h; tot += h; }
        if (tot <= 0.0) return CF_PENALTY;
        /* Weighted ROI support floor: `nin` above only counts in-FOV samples regardless of
           weight, so a pose that pushes the high-weight ROI out of FOV could pass it on
           zero-weight head voxels. Require the in-FOV weighted mass (tot) to be >=10% of the
           total weighted mass. Only active when weighted (fwt_total>0); no effect unweighted. */
        if (c->fwt && tot < 0.10 * c->fwt_total) return CF_PENALTY;
        double itot = 1.0/tot, bc = 0.0;
        for (int i = 0; i < NB; i++) {
            if (Pf[i] <= 0.0) continue;
            double pfi = Pf[i]*itot;
            for (int j = 0; j < NB; j++) {
                double h = acc[i*NB+j]; if (h <= 0.0 || Pm[j] <= 0.0) continue;
                bc += sqrt((h*itot) * pfi * (Pm[j]*itot));
            }
        }
        return bc;   /* minimize BC == maximize fixed/moving dependence */
    }

    /* CF_COST_CR: directional correlation ratio, fixed bins define the iso-sets.
       Parallelized with one preallocated accumulator block per worker (ctx->cr_acc),
       each block = bn[K], bs[K], bq[K] histograms + {N, Sy, Syy, nin} scalars. Blocks
       are combined in a fixed thread order, so the result is deterministic at a given
       thread count, and the nth==1 path reproduces the serial accumulation exactly
       (byte-identical), preserving the frozen -p 1 behavior. Different thread counts
       differ only by floating-point regrouping (≪ the 0.05 mm agreement floor). */
    int K = c->K;
    int stride = 3*K + 4;
    const int NC = CF_CR_NCHUNK;
    const float *fwt = c->fwt;   /* NULL -> weight 1.0 (default, bit-identical) */
    double *acc = c->cr_acc;
    memset(acc, 0, (size_t)NC * stride * sizeof(double));
#ifdef _OPENMP
    int par = (!c->serial_cost && ns > CF_CR_PAR_MIN);
#endif
    /* Each chunk cix owns sample range [cix*ns/NC, (cix+1)*ns/NC), summed in order
       into its own block. Chunk boundaries and per-chunk order are independent of the
       thread count, so the combined result below is bit-identical for any -p N. */
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) if(par)
#endif
    for (int cix = 0; cix < NC; cix++) {
        int lo = (int)((int64_t)cix * ns / NC);
        int hi = (int)((int64_t)(cix + 1) * ns / NC);
        double *B = acc + (size_t)cix * stride;
        double *bn = B, *bs = B + K, *bq = B + 2*K;
        double lN = 0, lSy = 0, lSyy = 0; long lnin = 0;
        for (int s = lo; s < hi; s++) {
            double ix, iy, iz; cf_apply(&gam, c->sx[s], c->sy[s], c->sz[s], &ix, &iy, &iz);
            int ok; float v = cf_trilerp(Lm, ix, iy, iz, &ok);
            if (!ok) continue;
            double mv = v;
            lnin++;   /* raw in-FOV count for the overlap floor (weight-independent) */
            int k = c->fbin[s];
            double w = fwt ? fwt[s] : 1.0, wmv = w * mv;
            bn[k] += w; bs[k] += wmv; bq[k] += wmv*mv;
            lN += w; lSy += wmv; lSyy += wmv*mv;
        }
        B[3*K] = lN; B[3*K+1] = lSy; B[3*K+2] = lSyy; B[3*K+3] = (double)lnin;
    }
    /* combine chunk blocks into block 0 in fixed chunk order (deterministic) */
    for (int cix = 1; cix < NC; cix++) {
        const double *Bt = acc + (size_t)cix * stride;
        for (int i = 0; i < stride; i++) acc[i] += Bt[i];
    }
    double *bn = acc, *bs = acc + K, *bq = acc + 2*K;
    double N = acc[3*K], Sy = acc[3*K+1], Syy = acc[3*K+2]; long nin = (long)acc[3*K+3];
    if (nin < (long)(0.10*ns) || nin < 16) return CF_PENALTY;
    /* Weighted ROI support floor (see the HEL path): require the in-FOV weighted mass N to be
       >=10% of the total sample weight so a pose can't push the high-weight ROI out of FOV and
       pass the raw-count floor on zero-weight voxels. Also guards N==0 (all-weight-zero -> NaN
       below). No effect unweighted (fwt NULL). */
    if (c->fwt && N < 0.10 * c->fwt_total) return CF_PENALTY;
    if (!(N > 0.0)) return CF_PENALTY;
    double varTot = Syy/N - (Sy/N)*(Sy/N);
    if (!(varTot > 1e-12)) return CF_PENALTY;   /* '!(x>eps)' also rejects NaN (weighted N->0) */
    double within = 0.0;
    for (int k = 0; k < K; k++) {
        /* Drop only truly-empty iso-set bins. `<= 0.0` (not the old `< 1.0`) makes the CR cost
           SCALE-INVARIANT for weighted bn[k] while staying bit-identical unweighted, where bn[k]
           is an exact integer count (0 <-> empty either way). A positive-but-tiny weighted bin has
           a well-defined vk; its bn[k]*vk contribution scales with the weight and cancels in the
           final ratio, so multiplying the weight image by a constant leaves the fit unchanged. */
        if (bn[k] <= 0.0) continue;
        double vk = bq[k]/bn[k] - (bs[k]/bn[k])*(bs[k]/bn[k]);
        if (vk < 0) vk = 0;
        within += bn[k]*vk;
    }
    return within / (N * varTot);   /* = 1 - eta^2, minimize */
}

/* Fraction of the fixed foreground sample cloud covered by the moving FOV. For
 * comparing two translations of the same moving FOV this ranks identically to
 * intersection/moving-FOV coverage; the fixed denominator matches the cost samples. */
static double cf_sample_coverage(const double p[12]) {
    cf_ctx *c = g_cf;
    mat44 F2M = cf_affine(c, p);
    mat44 gam = nifti_mat44_mul(c->Lm->w2v, nifti_mat44_mul(F2M, c->Sf_lvl));
    long nin = 0;
    for (int s = 0; s < c->ns; s++) {
        double ix, iy, iz;
        cf_apply(&gam, c->sx[s], c->sy[s], c->sz[s], &ix, &iy, &iz);
        int ok;
        (void)cf_trilerp(c->Lm, ix, iy, iz, &ok);
        if (ok) nin++;
    }
    return c->ns ? (double)nin / c->ns : 0.0;
}

/* All initialization/candidate choices use this SAME overlap-aware dependence
 * score. Keep the arithmetic order stable: raw minimized cost times overlap would
 * reward low overlap, while (1-cost)*overlap rewards useful statistical support. */
static double cf_dependence_overlap_score(double cost, double overlap) {
    double dep = 1.0 - cost;
    if (!cf_finite(dep) || dep < 0.0 || cost >= CF_PENALTY) dep = 0.0;
    if (dep > 1.0) dep = 1.0;
    return dep * overlap;
}

/* NEWUOA callback (cf_cost_eval already counts the evaluation). */
static double cf_ufunc(int n, double *x) {
    (void)n;
    double p[12]; cf_expand(g_cf, x, p);
    return cf_cost_eval(p);
}

/*==========================================================================*/
/* sample set + level stage setup                                            */
/*==========================================================================*/

/* Build the fixed sample set (foreground voxels) for the current fixed level and
 * precompute per-cost statistics. Frees any previous set. */
static void cf_free_samples(cf_ctx *c) {
    free(c->sx); free(c->sy); free(c->sz); free(c->fval); free(c->fbin); free(c->cr_acc); free(c->fwt);
    c->sx=NULL; c->sy=NULL; c->sz=NULL; c->fval=NULL; c->fbin=NULL; c->cr_acc=NULL; c->fwt=NULL;
    c->ns=0; c->K=0; c->fwt_total=0.0; c->fmean=0.0; c->fvar=0.0;
}

static int cf_make_samples(cf_ctx *c, int K) {
    const cf_level *Lf = c->Lf;
    size_t nv = (size_t)Lf->nx * Lf->ny * Lf->nz;
#ifdef COREG_FAST_TEST_ALLOC
    int inject_partial_oom =
        (++g_cf_test_sample_calls == g_cf_test_fail_sample_call);
    if (inject_partial_oom) g_cf_test_fail_sample_call = 0; /* fail once */
#endif
    /* A failed rebuild must leave an EMPTY, internally consistent sample state.
       This is load-bearing for multi-start: a later strategy must never polish a
       prior winner with stale ns against partial replacement arrays. */
    cf_free_samples(c);
    /* foreground threshold: 5% of dynamic range above min */
    double mn = DBL_MAX, mx = -DBL_MAX;
    for (size_t i = 0; i < nv; i++) { float v = Lf->data[i]; if (v<mn) mn=v; if (v>mx) mx=v; }
    if (!(mx > mn)) return 1;
    double thr = mn + c->thr_frac*(mx - mn);   /* thr_frac validated once at the boundary */
    /* count */
    int ns = 0;
    for (size_t i = 0; i < nv; i++) if (Lf->data[i] > thr) ns++;
    if (ns < 16) return 1;
    c->sx = malloc(sizeof(int32_t)*ns);
#ifdef COREG_FAST_TEST_ALLOC
    c->sy = inject_partial_oom ? NULL : malloc(sizeof(int32_t)*ns);
#else
    c->sy = malloc(sizeof(int32_t)*ns);
#endif
    c->sz = malloc(sizeof(int32_t)*ns); c->fval = malloc(sizeof(float)*ns);
    c->fbin = malloc(sizeof(int)*ns);
    if (!c->sx||!c->sy||!c->sz||!c->fval||!c->fbin) {
        cf_free_samples(c);
        return 1;
    }
    /* Optional per-sample weight: the weight level shares Lf's grid, so the sample voxel
       index below reads it directly. Absent (Lw==NULL) leaves fwt NULL -> unweighted. */
    if (c->Lw) {
        c->fwt = malloc(sizeof(float)*ns);
        if (!c->fwt) { cf_free_samples(c); return 1; }
    }
    double inv = (mx > mn) ? (K - 1) / (mx - mn) : 0.0;
    int idx = 0; double sum = 0;
    for (int z = 0; z < Lf->nz; z++)
      for (int y = 0; y < Lf->ny; y++)
        for (int x = 0; x < Lf->nx; x++) {
            size_t vi = x + (size_t)y*Lf->nx + (size_t)z*Lf->nx*Lf->ny;
            float v = Lf->data[vi];
            if (v <= thr) continue;
            c->sx[idx]=x; c->sy[idx]=y; c->sz[idx]=z; c->fval[idx]=v;
            int k = (int)((v - mn)*inv + 0.5); if (k<0) k=0; if (k>K-1) k=K-1;
            c->fbin[idx]=k;
            if (c->fwt) { float w = c->Lw->data[vi]; c->fwt[idx] = (w > 0.0f) ? w : 0.0f; }
            sum += v; idx++;
        }
    c->fwt_total = 0.0;
    double fwt_max = 0.0;
    if (c->fwt) for (int s = 0; s < ns; s++) { c->fwt_total += c->fwt[s];
                                               if (c->fwt[s] > fwt_max) fwt_max = c->fwt[s]; }
    /* Reject a weight with no positive support over the fixed foreground ONCE here (AFNI errors
       if -weight is never positive; we scope it to the sampled foreground since only foreground
       voxels enter the cost). With the [0,1]-normalized graded weight a real ROI reads up to ~1.0
       over the foreground, so require at least one foreground sample above a small margin — this
       rejects a mask whose positive voxels lie only in the fixed BACKGROUND (foreground reads ~0)
       and subsumes the all-zero-weight case (`wmax==0` skips the remap, so `fwt_max==0`). The
       per-pose weighted-coverage floor (N >= 10% of fwt_total) is retained for the accepted case. */
    if (c->fwt && !(fwt_max > 1e-3)) {
        fprintf(stderr, "coreg fast: -weight is never positive over the fixed foreground\n");
        cf_free_samples(c);
        return 1;
    }
    c->ns = ns; c->K = K;
    double mean = sum/ns; c->fmean = mean;
    double var = 0; for (int s=0;s<ns;s++){ double d=c->fval[s]-mean; var+=d*d; } var/=ns;
    c->fvar = var;
    /* CR/HEL accumulator scratch: CF_CR_NCHUNK fixed blocks, reused across this level's evals
       (thread-count-independent → deterministic reduction). Sized for the larger of the CR block
       (3K+4) and the Hellinger joint-histogram block (NB*NB+1). */
    size_t blk = 3*(size_t)K + 4;
    size_t hblk = (size_t)CF_HEL_NBIN*CF_HEL_NBIN + 1;
    if (hblk > blk) blk = hblk;
    c->cr_acc = (double *)malloc((size_t)CF_CR_NCHUNK * blk * sizeof(double));
    if (!c->cr_acc) { cf_free_samples(c); return 1; }
    return 0;
}

/* Set the active free-dof subset for a stage. dof: 6/7/9/12. */
static void cf_set_dof(cf_ctx *c, int dof, const double base[12]) {
    memcpy(c->base, base, sizeof(double)*12);
    c->global_scale = 0; c->nfree = 0;
    int idx[12], n = 0;
    idx[n++]=0; idx[n++]=1; idx[n++]=2;          /* translation */
    if (dof >= 6) { idx[n++]=3; idx[n++]=4; idx[n++]=5; }  /* rotation */
    if (dof == 7) { idx[n++]=6; c->global_scale = 1; }     /* global scale */
    if (dof >= 9) { idx[n++]=6; idx[n++]=7; idx[n++]=8; }  /* per-axis scale */
    if (dof >= 12){ idx[n++]=9; idx[n++]=10; idx[n++]=11; }/* shear */
    for (int i=0;i<n;i++) c->free_idx[i]=idx[i];
    c->nfree = n;
}

/* Run one NEWUOA refinement from base params, writing the result back into base.
 * powell_newuoa() returns the number of function calls, or a NEGATIVE code on a bad
 * argument or a workspace allocation failure (-7). Any negative return is propagated
 * via c->opt_err so the driver can fail the whole estimate rather than silently
 * returning an unrefined (or partially optimized) fit as success. */
static double cf_refine(cf_ctx *c, int dof, double base[12], double rstart, double rend, int maxcall) {
    cf_set_dof(c, dof, base);
    double x[12];
    for (int i = 0; i < c->nfree; i++) x[i] = base[c->free_idx[i]];
    int nc = powell_newuoa(c->nfree, x, rstart, rend, maxcall, cf_ufunc);
    if (nc < 0) c->opt_err = 1;
    if (dof > c->dof_run) c->dof_run = dof;   /* record the real highest DOF executed */
    double p[12]; cf_expand(c, x, p);
    memcpy(base, p, sizeof(double)*12);
    return cf_cost_eval(p);
}

/* One independent fastx coarse fit. The job owns its mutable cf_ctx fields and histogram
 * scratch; fixed samples and pyramid levels are shared read-only. Cost evaluation is kept
 * serial inside the job so the outer four-way OpenMP loop owns the available parallelism and
 * cannot create nested teams whose worker TLS lacks this job's callback context. */
typedef struct {
    cf_ctx c;
    double seed[12];
    double result[12];
    double result_cost;
    double rstart, rend;
    int cdof, coarse_search;
    float mfac, afac;
} cf_coarse_job;

static void cf_run_coarse_job(cf_coarse_job *j) {
    const int NTOP = 3;
    cf_ctx *c = &j->c;
    c->serial_cost = 1;
    g_cf = c;
    float old_mfac, old_afac;
    powell_get_mfac(&old_mfac, &old_afac);
    powell_set_mfac(j->mfac, j->afac);

    double bp[12]; memcpy(bp, j->seed, sizeof bp);
    cf_nms_candidate top[3]; memset(top, 0, sizeof top);
    for (int i=0; i<NTOP; i++) top[i].c = CF_PENALTY;
    if (j->coarse_search) {
        for (int a=0; a<5; a++) for (int b=0; b<5; b++) for (int cc=0; cc<5; cc++) {
            double p[12]; memcpy(p, bp, sizeof p);
            p[3]=bp[3]+CF_ANG[a]*CF_DEG2RAD/CF_PS_ROT;
            p[4]=bp[4]+CF_ANG[b]*CF_DEG2RAD/CF_PS_ROT;
            p[5]=bp[5]+CF_ANG[cc]*CF_DEG2RAD/CF_PS_ROT;
            p[6]=p[7]=p[8]=0.0;
            double cv = cf_cost_eval(p);
            cf_insert_sorted_candidate(top, NTOP, p, cv, (size_t)((a*5+b)*5+cc));
        }
    } else {
        memcpy(top[0].p, bp, sizeof bp); top[0].c=cf_cost_eval(bp); top[0].order=0;
    }

    double seedbest = CF_PENALTY;
    memcpy(j->result, bp, sizeof j->result);
    for (int t=0; t<NTOP; t++) {
        if (top[t].c >= CF_PENALTY) continue;
        double p[12]; memcpy(p, top[t].p, sizeof p);
        double cv = cf_refine(c, j->cdof, p, j->rstart, j->rend, 150);
        if (cv < seedbest) { seedbest=cv; memcpy(j->result, p, sizeof j->result); }
    }
    j->result_cost = seedbest;
    powell_newuoa_free_threadlocal();
    powell_set_mfac(old_mfac, old_afac);
    g_cf = NULL;
}

/*==========================================================================*/
/* driver                                                                    */
/*==========================================================================*/

/* intensity-weighted centroid (positive voxels) -> world coords */
static int cf_com_world(const float *d, int nx, int ny, int nz, const mat44 *S,
                        double *wx, double *wy, double *wz) {
    double sw=0, cx=0, cy=0, cz=0;
    for (int z=0;z<nz;z++) for (int y=0;y<ny;y++) for (int x=0;x<nx;x++) {
        double v = d[x + (size_t)y*nx + (size_t)z*nx*ny];
        if (v > 0 && v <= FLT_MAX) { sw+=v; cx+=v*x; cy+=v*y; cz+=v*z; }
    }
    if (sw <= 0) return 1;
    cf_apply(S, cx/sw, cy/sw, cz/sw, wx, wy, wz);
    return 0;
}

/* Point the cost context at a moving-pyramid level and cache that level's intensity
 * range (m_bg/m_top drive the HEL moving-axis binning and the LS background fill).
 * Used at each pyramid level. */
static void cf_use_moving_level(cf_ctx *c, const cf_level *Lmlv) {
    c->Lm = Lmlv;
    double mn = DBL_MAX, mx = -DBL_MAX;
    size_t n = (size_t)Lmlv->nx * Lmlv->ny * Lmlv->nz;
    for (size_t i = 0; i < n; i++) { double v = Lmlv->data[i]; if (v<mn) mn=v; if (v>mx) mx=v; }
    c->m_bg = (mn < DBL_MAX) ? mn : 0.0;
    c->m_top = (mx > mn) ? mx : (c->m_bg + 1.0);
}

/* Build a whole moving pyramid (levels 2->1->0) from a source image plus its brightness-centroid
 * seed, holding only ONE full-resolution copy at a time (extract -> build finest -> free full-res
 * -> quadrature-add coarser levels). Returns 0 on success; on failure frees any partially-built
 * levels (leaving them NULL-safe for the caller's cleanup) and returns 1. This is deterministic
 * construction OFF the numerical cost path, so it does not affect -p1==-pN. */
static int cf_build_moving_pyramid(const nifti_image *src, int nx, int ny, int nz,
                                   const mat44 *Sm, int use_cmass,
                                   const double SEP[3], const double FWHM[3],
                                   double com[3], int *have_com, cf_level L[3]) {
    float *full = cf_extract_float(src, NULL);
    if (!full) return 1;
    *have_com = use_cmass && !cf_com_world(full, nx, ny, nz, Sm, &com[0], &com[1], &com[2]);
    int ok = (cf_build_level(full, nx, ny, nz, Sm, SEP[2], FWHM[2], &L[2]) == 0);
    free(full);
    for (int lv = 1; lv >= 0 && ok; lv--) {
        double add = sqrt(FWHM[lv]*FWHM[lv] - FWHM[lv+1]*FWHM[lv+1]);
        if (cf_build_level(L[lv+1].data, L[lv+1].nx, L[lv+1].ny, L[lv+1].nz, &L[lv+1].v2w,
                           SEP[lv], add, &L[lv])) ok = 0;
    }
    if (!ok) { for (int i = 0; i < 3; i++) cf_free_level(&L[i]); return 1; }
    return 0;
}

/* Has the base been masked/skull-stripped — i.e. is its background hard-set to one value?
 * Measured as the fraction of voxels sitting exactly AT the minimum: a masked volume has a
 * huge exact-valued background delta, while a real acquisition or a whole-head template has a
 * noise/interpolation floor spread over many values. Gates the HEL joint-foreground drop (see
 * CF_HEL_MDARK_FRAC). Computed once on the FULL-RES base, off the cost path and independent of
 * pose, thread count and pyramid level, so -p1==-pN stays bit-identical. */
static int cf_base_is_hard_zeroed(const float *d, size_t n) {
    if (!d || n == 0) return 0;
    float mn = d[0];
    for (size_t i = 1; i < n; i++) if (d[i] < mn) mn = d[i];
    size_t nmin = 0;
    for (size_t i = 0; i < n; i++) if (d[i] == mn) nmin++;
    return (double)nmin > CF_HEL_MDARK_BASEFRAC * (double)n;
}

int coreg_fast_estimate(const nifti_image *moving, const nifti_image *fixed,
                        const coreg_fast_opts *opts, coreg_fast_result *result) {
    coreg_fast_opts O = opts ? *opts : coreg_fast_opts_default();
    /* The estimator is configured solely by `opts`. Ambient environment variables must not
     * override an embedding application's explicit request or make a saved affine impossible
     * to reproduce. Validate/normalize the public options once at the boundary. */
    if (O.cost != CF_COST_LS && O.cost != CF_COST_CR && O.cost != CF_COST_HEL &&
        O.cost != CF_COST_HEL_CR)
        O.cost = CF_COST_HEL_CR;
    /* Snap max_dof to the documented schedule set {6,7,9,12} (8/10/11 have no stage). */
    O.max_dof = (O.max_dof < 7) ? 6 : (O.max_dof < 9) ? 7 : (O.max_dof < 12) ? 9 : 12;
#ifdef AL_PROFILE
    g_prof_blur = g_prof_resample = 0.0;   /* per-call reset (accumulators are file-global) */
#endif
    double t0 = cf_wtime();
    if (!moving || !fixed || !result) return 1;
    if (cf_dims_ok(fixed, "fixed") || cf_dims_ok(moving, "moving")) return 1;
    /* The optional weight is a fixed(base)-space image: it must be a valid 3D volume with
       dimensions identical to `fixed` so it shares the fixed pyramid grids voxel-for-voxel. */
    if (O.weight && (cf_dims_ok(O.weight, "weight") ||
                     O.weight->nx != fixed->nx || O.weight->ny != fixed->ny ||
                     O.weight->nz != fixed->nz)) {
        fprintf(stderr, "coreg fast: -weight image dims must match the fixed image\n");
        return 1;
    }
    /* The LS cost does not read the per-sample weight (it is a diagnostic, not CLI-reachable).
       Reject weight+LS at the boundary rather than silently ignoring the weight. */
    if (O.weight && O.cost == CF_COST_LS) {
        fprintf(stderr, "coreg fast: -weight is not supported with the LS cost\n");
        return 1;
    }
    int have_weight = (O.weight != NULL);
    mat44 Sf, Sm;
    /* Unified no-form policy (shared with al_register/-sym/-com/apply): coded sform/qform
       when usable, else a pixdim-centered frame — so both-codes-zero inputs still register. */
    al_image_xform_or_pixdim(fixed, &Sf, NULL);
    al_image_xform_or_pixdim(moving, &Sm, NULL);

    int fnx=fixed->nx, fny=fixed->ny, fnz=fixed->nz;
    int mnx=moving->nx, mny=moving->ny, mnz=moving->nz;

    /* The weight shares the fixed grid VOXEL-FOR-VOXEL (its pyramid is built with Sf and it is
       sampled by the fixed index, never resampled), so the two grids must be the SAME grid, not
       merely close: a fractional-voxel origin offset would shift the physical ROI by up to a
       voxel while still "covering" the FOV. Map the 8 fixed-index corners through each frame and
       reject if any world position diverges by more than a small fraction of a voxel — a header-
       rounding tolerance (a real stationary-space mask carries the fixed sform/qform and agrees
       to FP noise; a resampled/reoriented/shifted one does not). */
    if (have_weight) {
        mat44 Sw; al_image_xform_or_pixdim(O.weight, &Sw, NULL);
        double vx = cf_colnorm(&Sf,0), vy = cf_colnorm(&Sf,1), vz = cf_colnorm(&Sf,2);
        double vmin = vx; if (vy < vmin) vmin = vy; if (vz < vmin) vmin = vz;
        double tol = 0.05 * vmin; if (!(tol > 0.0)) tol = 0.05;   /* ~1/20 voxel: FP noise, not sub-voxel shift */
        int bad = 0;
        for (int c = 0; c < 8 && !bad; c++) {
            double ix = (c&1)?fnx-1:0, iy = (c&2)?fny-1:0, iz = (c&4)?fnz-1:0;
            double fx,fy,fz, wx,wy,wz;
            cf_apply(&Sf, ix,iy,iz, &fx,&fy,&fz);
            cf_apply(&Sw, ix,iy,iz, &wx,&wy,&wz);
            double dx=fx-wx, dy=fy-wy, dz=fz-wz;
            if (sqrt(dx*dx+dy*dy+dz*dz) > tol) bad = 1;
        }
        if (bad) {
            fprintf(stderr, "coreg fast: -weight image is in a different world frame than the "
                            "fixed image (matching dims but sform/qform disagree)\n");
            g_cf = NULL;
            return 1;
        }
    }

    cf_ctx ctx; memset(&ctx, 0, sizeof ctx);
    /* HEL_CR is an orchestration mode, not a cost evaluator: it uses HEL everywhere except
       while explicitly fitting its CR coarse candidates below. */
    ctx.cost = (O.cost == CF_COST_HEL_CR) ? CF_COST_HEL : O.cost;
    ctx.thr_frac = 0.05;  /* fixed foreground threshold (fraction of dynamic range) */
    cf_apply(&Sf, (fnx-1)*0.5, (fny-1)*0.5, (fnz-1)*0.5, &ctx.cx, &ctx.cy, &ctx.cz);
    g_cf = &ctx;

    /* Pyramid schedule: 8 -> 4 -> 2 mm. FWHM[] is the target effective smoothing
       (full-resolution-equivalent) at each level. The pyramid is built HIERARCHICALLY:
       the finest level (2 mm) is blurred+decimated from full resolution, and each
       coarser level is an additional blur + decimation of the next-finer level, so the
       expensive coarse smoothing runs on already-small data instead of re-blurring the
       full-resolution volume three times. Per-step added blur combines in quadrature. */
    const double SEP[3]  = { 8.0, 4.0, 2.0 };
    const double FWHM[3] = { 8.0, 5.0, 3.0 };
    const int nlv = 3;
    int levels_done = 0;
    double best[12];
    double best_all[12]; double cost_all = CF_PENALTY; int have_all = 0;  /* winner across starts */
    double cands[6][12]; int ncand = 0; memset(best, 0, sizeof best);
    double final_cost = CF_PENALTY;   /* local; written to *result only on full success */

    /* Build each image's whole pyramid and FREE its full-res copy BEFORE extracting the
       other image, so only ONE full-resolution volume is ever held at a time (plus the
       small pyramids and one blur scratch). The moving COM is taken from its full-res
       buffer just before that buffer is freed. */
    cf_level Lf[3] = {{0}}, Lm[3] = {{0}}, Lwt[3] = {{0}};  /* Lwt: fine-stage weight pyramid */
    double comm[3] = {0};            /* moving brightness centroid (world mm), for the -com seed */
    int build_ok = 1, have_com_m = 0;
#ifdef AL_PROFILE
    double _tb0 = cf_wtime();
#endif
    /* fixed pyramid */
    float *ffull = cf_extract_float(fixed, NULL);
    if (!ffull) { g_cf = NULL; fprintf(stderr, "coreg fast: fixed unsupported datatype / OOM\n"); return 1; }
    /* Decide the joint-foreground gate here, while the full-res base is still in hand. */
    ctx.hel_mdark = cf_base_is_hard_zeroed(ffull, (size_t)fnx*fny*fnz);
    if (O.verbose)
        fprintf(stderr, "[coreg fast] base %s -> HEL joint-foreground drop %s\n",
                ctx.hel_mdark ? "masked/skull-stripped" : "whole-head",
                ctx.hel_mdark ? "ON" : "OFF");
    if (cf_build_level(ffull, fnx,fny,fnz, &Sf, SEP[2], FWHM[2], &Lf[2])) build_ok = 0;
    free(ffull); ffull = NULL;
    for (int lv = 1; lv >= 0 && build_ok; lv--) {
        double add = sqrt(FWHM[lv]*FWHM[lv] - FWHM[lv+1]*FWHM[lv+1]);
        if (cf_build_level(Lf[lv+1].data, Lf[lv+1].nx,Lf[lv+1].ny,Lf[lv+1].nz, &Lf[lv+1].v2w, SEP[lv], add, &Lf[lv])) build_ok = 0;
    }
    /* moving pyramid — extracted only after the fixed full-res copy is freed */
    if (build_ok && cf_build_moving_pyramid(moving, mnx,mny,mnz, &Sm, O.use_cmass, SEP, FWHM,
                                            comm, &have_com_m, Lm))
        build_ok = 0;
    /* -nocmass (O.use_cmass==0) already zeroed have_com_m above by short-circuiting
       cf_com_world(), so the supplied-affine start is all that remains. */
    /* Weight pyramid: only the FINEST level (2 mm, Lwt[2]) is ever consumed — the cost applies the
       weight solely at lv==2 (see the `ctx.Lw = have_weight && lv>=2` gate below), so the 8 mm
       coarse capture, the 4 mm global-scale bracket, and seed selection stay whole-head. The weight
       shares `fixed`'s dims and frame, so building Lwt[2] from the full-res weight with the same
       (Sf, SEP[2], FWHM[2]) chain yields a grid identical to Lf[2] voxel-for-voxel (blurring softens
       a binary mask, which is desirable). Lwt[0]/Lwt[1] are never built (NULL-safe frees below). */
    if (build_ok && have_weight) {
        float *wfull = cf_extract_float(O.weight, NULL);
        if (!wfull) { build_ok = 0; fprintf(stderr, "coreg fast: weight unsupported datatype / OOM\n"); }
        else {
            /* Clamp negatives and NaN to 0 at FULL RESOLUTION, before the pyramid blur — a
               negative sample would otherwise spread through the Gaussian and depress neighbouring
               positive support (a weight is conceptually non-negative). `!(w > 0)` also maps NaN->0. */
            size_t wnv = (size_t)fnx * fny * fnz;
            float wmax = 0.0f;
            for (size_t i = 0; i < wnv; i++) {
                /* Clamp negatives/NaN AND +Inf to 0 (magnitude guard, since this TU is
                   -ffast-math and a +Inf would poison the remap: wmax=+Inf -> inv=0 ->
                   Inf*0=NaN for that voxel and every finite ROI value collapses to the floor). */
                if (!(wfull[i] > 0.0f && wfull[i] <= FLT_MAX)) wfull[i] = 0.0f;
                else if (wfull[i] > wmax) wmax = wfull[i];
            }
            /* AFNI-style whole-head anchor (the wisdom of 3dAllineate's -autobox default): do NOT
               zero the cost outside the ROI. Restricting the cost to the mask alone leaves GLOBAL
               (isotropic) extent underdetermined — the outer head boundary that fixes absolute size
               has been masked away — so the masked cost slides monotonically toward shrinking the
               source until its bright scalp is dragged into the ROI (the FLIRT Fig-1 degeneracy).
               The failure is ill-conditioned, hence compiler/FP-sensitive: it surfaced as a large
               clang-vs-gcc divergence. AFNI's -weight avoids it not with a floor but by supplying
               a whole-head weight that keeps the scalp ATTENUATED (nonzero), so scale stays
               anchored; we mirror AFNI exactly here — normalize the supplied weight to [0,1] over
               the whole fixed grid (divide by its max) and use it graded. Normalizing by the max
               keeps the fit invariant to a global weight scaling (HEL/CR are scale-invariant). */
            if (wmax > 0.0f) {
                /* Divide directly rather than multiply by 1/wmax: a subnormal/tiny wmax makes
                   1.0f/wmax overflow to +Inf, and then 0*Inf=NaN (background) / tiny*Inf=+Inf
                   (foreground) would poison the remap (+Inf is not caught downstream). Because
                   every wfull[i] <= wmax, the ratio is in [0,1] and CANNOT overflow. */
                for (size_t i = 0; i < wnv; i++)
                    wfull[i] = wfull[i] / wmax;
            }
            if (cf_build_level(wfull, fnx,fny,fnz, &Sf, SEP[2], FWHM[2], &Lwt[2])) build_ok = 0;
            free(wfull);
        }
    }
    if (!build_ok) {
        for (int i=0;i<3;i++){ cf_free_level(&Lf[i]); cf_free_level(&Lm[i]); cf_free_level(&Lwt[i]); }
        g_cf = NULL;
        fprintf(stderr, "coreg fast: pyramid construction failed\n");
        return 1;
    }
#ifdef AL_PROFILE
    fprintf(stderr, " [coreg profile] pyramid build: %.3f s (blur %.3f, resample %.3f)\n",
            (cf_wtime()-_tb0), g_prof_blur, g_prof_resample);
#endif

    /* ---- optimize coarse -> fine over the prebuilt pyramid, once per coarse start ---- */
    /* Multi-start ONLY pays on a hard-zeroed base. Measured over 35 pairs: on a whole-head
       base the extra starts gain a mean 0.00014 (max 0.00078) in final cost -- pure noise --
       whereas on a hard-zeroed base they gain a mean 0.00857 (max 0.04586). Reuse the base
       detector that already gates the joint-foreground skip, so a whole-head base keeps the
       single-pass runtime and a stripped/defaced base gets the full search. */
    /* The mixed default uses all three trajectories. Explicit fasthel/fastcr must remain
       genuinely single-cost selectors, so they use only rigid + scale-bracketed strategies
       and never inject the other cost through strategy 2. */
    int nstrat = ctx.hel_mdark
                   ? ((O.cost == CF_COST_HEL_CR) ? CF_NSTRAT :
                      (O.cost == CF_COST_HEL || O.cost == CF_COST_CR) ? 2 : 1)
                   : 1;
    int sample_err = 0;
#ifdef COREG_FAST_TEST_ALLOC
    g_cf_test_sample_calls = 0;
#endif
    for (int strat = 0; strat < nstrat; strat++) {
      ncand = 0;
      levels_done = 0;
      final_cost = CF_PENALTY;
      for (int lv = 0; lv < nlv; lv++) {
        /* Strategy 2 seeds the pyramid with a CORRELATION-RATIO coarse stage and then runs the
           ordinary HEL fine stages from it. CR and HEL fail on different pairs, so this is a
           cost-space multi-start to complement the DOF-space one; the winner is still chosen on
           the common HEL fine-level cost, so the comparison stays apples-to-apples. */
        ctx.cost = (strat == 2 && lv == 0 &&
                    (O.cost == CF_COST_HEL || O.cost == CF_COST_HEL_CR))
                       ? CF_COST_CR
                       : (O.cost == CF_COST_HEL_CR ? CF_COST_HEL : O.cost);
        ctx.Lf=&Lf[lv]; ctx.Sf_lvl=Lf[lv].v2w;
        /* Weight ONLY the finest level (lv==2, 2 mm); the 8 mm rigid coarse (lv==0) AND the
           4 mm global-scale bracket + its refine (lv==1) stay UNWEIGHTED so the graded ROI weight
           cannot hijack global-scale SELECTION. This extends "rigid coarse, scale later": the
           scale-later stage must also see the full head, or a brain-concentrated weight (e.g. an
           @SSwarper T1w weight applied cross-modal) picks the scale that shrinks the moving into
           the ROI (the FLIRT Fig-1 basin) — the 4 mm bracket scored a spurious ~1.13x isotropic
           shrink on T2w/T1w_MICCAI -> SSW, det 1.45, collapsing the fit. Scoring the bracket on the
           full head recovers the correct expand scale (AFNI-parity: T2w masked HEL 0.153 vs 0.095,
           det 0.69 vs 1.45), and the weighted 2 mm descent then refines WITHIN that basin without
           collapsing. Neutral (<0.001) on the non-collapsing weighted cases; unweighted fits are
           unaffected (have_weight==0 leaves ctx.Lw NULL). The weight's role is fine-stage steering,
           matching its doc: it must not re-select DOF the coarse pyramid already resolved. */
        ctx.Lw = (have_weight && lv >= 2) ? &Lwt[lv] : NULL;
        /* The moving level (ctx.Lm + m_bg/m_top) is set here for the fine levels; the coarse
           level (lv==0) sets it itself below. cf_make_samples reads only ctx.Lf. */
        if (lv != 0) cf_use_moving_level(&ctx, &Lm[lv]);
#ifdef AL_PROFILE
        int _ev0 = ctx.evals; double _tsamp0 = cf_wtime();
#endif
        int K = (int)lround(256.0 / SEP[lv]); if (K<32) K=32; if (K>256) K=256;
        if (cf_make_samples(&ctx, K)) { sample_err = 1; break; }
        double rstart = SEP[lv], rend = SEP[lv]*0.03;
#ifdef AL_PROFILE
        double _tsamp = cf_wtime() - _tsamp0, _topt0 = cf_wtime();
#endif

        if (lv == 0) {
          cf_use_moving_level(&ctx, &Lm[0]);
          /* Choose the supplied-affine frame or the exact `-com` recentered frame
             before the expensive orientation search. In the original world frame,
             applying `-com` is exactly a translation by the moving brightness
             centroid, so both starts can be scored against the SAME pyramid.

             HEL and CR are minimized unexplained-dependence coefficients
             (1 == independent), hence the joint selection score is:

                 (1 - cost) * fraction_of_fixed_foreground_inside_moving_FOV

             Maximizing this rewards both statistical alignment and useful overlap.
             Multiplying the minimized cost itself by coverage would perversely favor
             low overlap. In the single-cost paths and hard-zero multi-start, only the
             selected seed receives the coarse grid and local descent. Whole-head fastx
             deliberately keeps both frames for its HEL/CR competition. Forced `-com`
             arrives with a recentered header and
             O.use_cmass==0; `-nocmass` likewise keeps only its supplied frame. */
          /* The cheap mixed selector is for whole-head bases. On a hard-zeroed base its
             8 mm winner is NOT a safe 2 mm/12-DOF arbiter: adding that early-selected
             trajectory to the imported multi-start made it win at 7 DOF but lose badly
             after polish (T1w1mm->MNIstrip NCC 0.477 -> 0.145; MICCAI->avgstrip
             0.099 -> -0.037). Keep the three full-depth hard-zero strategies independent:
             historical HEL, scale-bracketed HEL, and CR-seeded HEL. */
          int compete = (O.cost == CF_COST_HEL_CR && !ctx.hel_mdark && strat == 0);
          double seeds[2][12]; int nseed = 1;
          memset(seeds[0], 0, sizeof seeds[0]);
          /* NB: the seed CHOICE below must be deterministic across thread counts (the byte-
             parity contract). It is, because HEL/CR use fixed-order chunked reductions — but
             cf_cost_eval's LS branch uses reduction(+:) and is NOT thread-stable, so do not
             wire a future LS caller to use_cmass expecting -p1==-pN seed selection. */
          if (O.use_cmass && have_com_m) {
              double hp[12] = {0}, cp[12] = {0};
              cp[0] = comm[0] / CF_PS_TRANS;
              cp[1] = comm[1] / CF_PS_TRANS;
              cp[2] = comm[2] / CF_PS_TRANS;
              if (compete) {
                  /* fastx keeps BOTH frames: each is independently fitted with HEL and CR,
                     yielding the four coarse candidates requested by the mixed mode. */
                  memcpy(seeds[1], cp, sizeof cp);
                  nseed = 2;
              } else {
                  double hc = cf_cost_eval(hp), cc = cf_cost_eval(cp);
                  double hov = cf_sample_coverage(hp), cov = cf_sample_coverage(cp);
                  double hs = cf_dependence_overlap_score(hc, hov);
                  double cs = cf_dependence_overlap_score(cc, cov);
                  if (cs > hs) memcpy(seeds[0], cp, sizeof cp);
                  if (O.verbose)
                      fprintf(stderr, "[coreg fast] initial affine cost=%.5f overlap=%.5f score=%.5f; "
                                      "COM cost=%.5f overlap=%.5f score=%.5f -> %s\n",
                              hc, hov, hs, cc, cov, cs, (cs > hs) ? "COM" : "affine");
              }
          }

          /* Search and refine the selected start. Whole-head fastx independently searches both
             available frames under both coarse costs; the hard-zero path instead gives HEL,
             scale-bracketed HEL, and CR their own full-depth strategies. */
          typedef cf_nms_candidate cand;
          const int NTOP = (strat == 0) ? 3 : 6;   /* pass 0 == historical */
          /* RIGID coarse (8 mm): lock global scale at identity and search only
             orientation+translation. Freeing scale this coarse lets the correlation-ratio
             cost commit to a spurious isotropic-shrink basin on wide-FOV / short-axis
             cross-modal inputs (observed: T2w->avg152T1 collapsed to ~1.2x). Global scale
             is introduced at 4 mm from this good rigid start; aniso/shear at 2 mm. This
             mirrors FLIRT (coarse rotation/translation search before scaling). Honors a
             lower max_dof request. See AGENTS.md "rigid coarse". */
          int cdof = (O.max_dof < 6) ? O.max_dof : 6;
          double bestc = CF_PENALTY;
          if (!compete) {
            /* Historical fast/fastcr path: keep its ordering and arithmetic unchanged. */
            for (int sd=0; sd<nseed; sd++) {
              double bp[12]; memcpy(bp,seeds[sd],sizeof bp);
              cand top[6]; memset(top, 0, sizeof top);
              for (int i=0;i<NTOP;i++) top[i].c=CF_PENALTY;
              cand pool[5*5*5*CF_CSC_N*CF_CZR_N]; size_t npool = 0;
              if (O.coarse_search) {
                  /* strat 0: rigid (scale locked at identity) -- the historical trajectory.
                     strat 1: also bracket an isotropic scale x z-ratio, searched JOINTLY with
                     orientation, because on a skull-stripped base the true orientation is only
                     reachable WITH scale (see CF_NSTRAT). Discrete and bounded, so it cannot
                     run away the way a free continuous coarse scale did. */
                  static const double CSC[CF_CSC_N] = {0.85, 0.95, 1.05};
                  static const double CZR[CF_CZR_N] = {0.85, 1.0};
                  /* A scale-bearing start is itself a scale fit: never inject it when the
                     public max_dof contract stops at rigid. Strategy 2 still contributes a
                     distinct rigid CR seed in that configuration. */
                  int explore_scale = (strat != 0 && O.max_dof >= 7);
                  int nsc = explore_scale ? CF_CSC_N : 1;
                  int nzr = explore_scale ? CF_CZR_N : 1;
                  for (int a=0;a<5;a++) for (int b=0;b<5;b++) for (int cc=0;cc<5;cc++)
                  for (int sI=0; sI<nsc; sI++) for (int zI=0; zI<nzr; zI++) {
                      double p[12]; memcpy(p,bp,sizeof p);
                      p[3]=bp[3]+CF_ANG[a]*CF_DEG2RAD/CF_PS_ROT;
                      p[4]=bp[4]+CF_ANG[b]*CF_DEG2RAD/CF_PS_ROT;
                      p[5]=bp[5]+CF_ANG[cc]*CF_DEG2RAD/CF_PS_ROT;
                      if (!explore_scale) { p[6]=p[7]=p[8]=0.0; }
                      else { p[6]=p[7]=(CSC[sI]-1.0)/CF_PS_SCALE;
                             p[8]=(CSC[sI]*CZR[zI]-1.0)/CF_PS_SCALE; }
                      double cv = cf_cost_eval(p);
                      if (strat != 0) {
                          memcpy(pool[npool].p, p, sizeof p);
                          pool[npool].c = cv;
                          pool[npool].order = npool;
                          npool++;
                      } else {
                          cf_insert_sorted_candidate(top, NTOP, p, cv, npool);
                      }
                  }
                  /* The exploratory grid's raw top-N tends to contain neighbours of one
                     minimum. Select globally by cost before greedily retaining distinct
                     basins; in-place replacement can break ordering and separation. */
                  if (strat != 0) {
                      int nsel = cf_select_angular_nms(
                          pool, npool, top, NTOP, 20.0,
                          CF_PS_ROT/CF_DEG2RAD, CF_PENALTY);
                      for (int t=nsel; t<NTOP; t++) top[t].c=CF_PENALTY;
                  }
              } else {
                  memcpy(top[0].p,bp,sizeof bp); top[0].c=cf_cost_eval(bp); top[0].order=0;
              }
              double seedbest = CF_PENALTY; double seedp[12]; memcpy(seedp,bp,sizeof seedp);
              int kmax_cand = (strat == 0) ? 1 : 3;
              for (int t=0;t<NTOP;t++) {
                  if (top[t].c >= CF_PENALTY) continue;
                  double p[12]; memcpy(p,top[t].p,sizeof p);
                  double cv = cf_refine(&ctx, cdof, p, rstart, rend, 150);
                  if (ncand < kmax_cand) memcpy(cands[ncand++], p, sizeof p);
                  if (cv < seedbest) { seedbest=cv; memcpy(seedp,p,sizeof seedp); }
              }
              if (O.verbose) fprintf(stderr, "[coreg fast]  seed %d refined best %.4f rot(%.1f,%.1f,%.1f)deg\n",
                  sd, seedbest, seedp[3]*CF_PS_ROT/CF_DEG2RAD, seedp[4]*CF_PS_ROT/CF_DEG2RAD, seedp[5]*CF_PS_ROT/CF_DEG2RAD);
              if (seedbest < bestc) {
                  bestc=seedbest; memcpy(best,seedp,sizeof best);
              }
            }
          } else {
              /* Four jobs in the normal auto mode:
                    HEL×affine, HEL×COM, CR×affine, CR×COM.
                 Each owns its callback context and accumulator. Their inner HEL/CR loops are
                 serial; OpenMP distributes whole optimizers across up to four workers, avoiding
                 nested-team oversubscription. Forced -com/-nocmass supplies one frame -> two jobs. */
              cf_coarse_job jobs[4]; memset(jobs, 0, sizeof jobs);
              float pmfac, pafac; powell_get_mfac(&pmfac, &pafac);
              size_t acc_blk = 3*(size_t)ctx.K + 4;
              size_t hel_blk = (size_t)CF_HEL_NBIN*CF_HEL_NBIN + 1;
              if (hel_blk > acc_blk) acc_blk = hel_blk;
              size_t acc_n = (size_t)CF_CR_NCHUNK * acc_blk;
              int njob = 0, alloc_failed = 0;
              for (int ci=0; ci<2; ci++) for (int sd=0; sd<nseed; sd++) {
                  cf_coarse_job *j = &jobs[njob++];
                  j->c = ctx;
                  j->c.cost = (ci == 0) ? CF_COST_HEL : CF_COST_CR;
                  j->c.evals = 0; j->c.opt_err = 0;
                  j->c.cr_acc = (double *)malloc(acc_n * sizeof(double));
                  if (!j->c.cr_acc) alloc_failed = 1;
                  memcpy(j->seed, seeds[sd], sizeof j->seed);
                  j->rstart=rstart; j->rend=rend; j->cdof=cdof;
                  j->coarse_search=O.coarse_search; j->mfac=pmfac; j->afac=pafac;
              }
              if (!alloc_failed) {
#ifdef _OPENMP
                  int job_threads = njob, max_threads = omp_get_max_threads();
                  if (job_threads > max_threads) job_threads = max_threads;
                  #pragma omp parallel for schedule(static) num_threads(job_threads)
#endif
                  for (int j=0; j<njob; j++) cf_run_coarse_job(&jobs[j]);
              } else {
                  ctx.opt_err = 1;
              }
              /* A worker that used the main thread leaves its TLS callback NULL. Restore the
                 driver's context before deterministic fixed-order aggregation and HEL scoring. */
              g_cf = &ctx; ctx.cost = CF_COST_HEL;
              double coarse_cand[4][12]; int ncoarse = 0;
              for (int j=0; j<njob; j++) {
                  if (!alloc_failed) {
                      ctx.evals += jobs[j].c.evals;
                      if (jobs[j].c.opt_err) ctx.opt_err = 1;
                      if (jobs[j].c.dof_run > ctx.dof_run) ctx.dof_run = jobs[j].c.dof_run;
                      if (jobs[j].result_cost < CF_PENALTY)
                          memcpy(coarse_cand[ncoarse++], jobs[j].result,
                                 sizeof coarse_cand[0]);
                      if (O.verbose)
                          fprintf(stderr, "[coreg fast]  fastx %s seed %d refined %.5f\n",
                                  jobs[j].c.cost == CF_COST_HEL ? "HEL" : "CR",
                                  j % nseed, jobs[j].result_cost);
                  }
                  free(jobs[j].c.cr_acc); jobs[j].c.cr_acc = NULL;
              }

              /* Judge EVERY candidate with the same HEL objective. Retain the established
                 overlap term: raw HEL alone can reward a partial-FOV pose that maps fixed
                 foreground outside the moving acquisition. Strict comparison makes ties choose
                 the earlier HEL candidate deterministically. */
              ctx.cost = CF_COST_HEL;
              double best_score = -1.0;
              for (int i=0; i<ncoarse; i++) {
                  double hc = cf_cost_eval(coarse_cand[i]);
                  double ov = cf_sample_coverage(coarse_cand[i]);
                  double score = cf_dependence_overlap_score(hc, ov);
                  if (O.verbose)
                      fprintf(stderr, "[coreg fast]  fastx candidate %d HEL=%.5f overlap=%.5f "
                                      "score=%.5f\n", i, hc, ov, score);
                  if (score > best_score) {
                      best_score = score; bestc = hc;
                      memcpy(best, coarse_cand[i], sizeof best);
                  }
              }
          }
          if (bestc >= CF_PENALTY) memcpy(best,seeds[0],sizeof best);
          /* The overall coarse winner MUST be candidate 0: the 4 mm stage starts from
             cands[0], and pass 0 carries exactly one candidate, so this is what makes pass 0
             reproduce the historical single-candidate trajectory. Extra (distinct) candidates
             follow it, for the exploratory pass only. */
          { double tmp[6][12]; int nt = 0, kmax = (strat == 0) ? 1 : 3;
            memcpy(tmp[nt++], best, sizeof best);
            for (int i = 0; i < ncand && nt < kmax; i++) {
                double d0=(cands[i][3]-best[3])*CF_PS_ROT/CF_DEG2RAD;
                double d1=(cands[i][4]-best[4])*CF_PS_ROT/CF_DEG2RAD;
                double d2=(cands[i][5]-best[5])*CF_PS_ROT/CF_DEG2RAD;
                if (sqrt(d0*d0+d1*d1+d2*d2) < 1.0) continue;   /* duplicate of best */
                memcpy(tmp[nt++], cands[i], sizeof cands[i]);
            }
            for (int i = 0; i < nt; i++) memcpy(cands[i], tmp[i], sizeof tmp[i]);
            ncand = nt; }
          if (O.verbose) fprintf(stderr, "[coreg fast]  coarse seeds=%d -> best %.4f rot(%.1f,%.1f,%.1f)deg\n",
              nseed, bestc,
              best[3]*CF_PS_ROT/CF_DEG2RAD, best[4]*CF_PS_ROT/CF_DEG2RAD, best[5]*CF_PS_ROT/CF_DEG2RAD);
        } else if (lv == 1) {
            /* Global scale is introduced HERE (the 8 mm coarse level was rigid). A
               continuous refine from identity cannot reach a large TRUE scale at the
               envelope edge — a 25% downscale has a weak gradient (the moving brain is
               smaller than the fixed sample cloud), so gscale_dn(0.75) stalls near 1.0.
               First bracket the global scale with a discrete line search AT the good,
               locked-in rigid pose: this seeds a real large scale difference without
               re-opening the coarse spurious-shrink basin (which needed the coarse
               orientation search to be attractive). Then warm-started 7-DOF refine. */
            if (ncand < 1) { memcpy(cands[0], best, sizeof best); ncand = 1; }
            double selbest = -1.0, selp[12];
            memcpy(selp, cands[0], sizeof selp);
            for (int ci = 0; ci < ncand; ci++) {
                double cp[12]; memcpy(cp, cands[ci], sizeof cp);
                if (O.max_dof >= 7) {
                    const double SCB[5] = {0.78, 0.9, 1.0, 1.12, 1.3};
                    double sbest = cp[6], cbest = CF_PENALTY;
                    for (int i = 0; i < 5; i++) {
                        double p[12]; memcpy(p, cp, sizeof p);
                        p[6] = p[7] = p[8] = (SCB[i]-1.0)/CF_PS_SCALE;
                        double cv = cf_cost_eval(p);
                        if (cv < cbest) { cbest = cv; sbest = (SCB[i]-1.0)/CF_PS_SCALE; }
                    }
                    cp[6] = cp[7] = cp[8] = sbest;
                }
                /* warm-started polish: half-level radius so the fit refines rather than
                   re-searches (a full-level radius lets the extra DOF wander off-brain). */
                double cv = cf_refine(&ctx, (O.max_dof<7)?O.max_dof:7, cp, SEP[lv]*0.5, SEP[lv]*0.02, 250);
                /* Same overlap-aware score the seed selector uses: raw cost alone can prefer a
                   pose that slides the base foreground out of the moving FOV. */
                double ov = cf_sample_coverage(cp);
                double sc = cf_dependence_overlap_score(cv, ov);
                if (sc > selbest) { selbest = sc; memcpy(selp, cp, sizeof selp); }
            }
            memcpy(best, selp, sizeof best);
        } else { /* lv == 2: only the 7-DOF refine here -- see the polish note below */
            cf_refine(&ctx, (O.max_dof>=7)?7:6, best, SEP[lv]*0.5,  SEP[lv]*0.02, 200);
        }
        /* Any optimizer allocation failure is fatal and atomic; do not spend more
           work on finer levels or later strategies when success is impossible. */
        if (ctx.opt_err) break;
#ifdef AL_PROFILE
        double _topt = cf_wtime() - _topt0;
        fprintf(stderr, " [coreg profile] level %g mm (%dx%dx%d): samples %.3f s (%d pts), "
                "optimize %.3f s (%d evals)\n",
                SEP[lv], Lf[lv].nx, Lf[lv].ny, Lf[lv].nz,
                _tsamp, ctx.ns, _topt, ctx.evals-_ev0);
#endif
        final_cost = cf_cost_eval(best);
        levels_done++;
        if (O.verbose) fprintf(stderr, "[coreg fast] start %d: level %g mm done, cost=%.5f evals=%d\n",
                               strat, SEP[lv], final_cost, ctx.evals);
      }
      if (sample_err || ctx.opt_err) break;
      /* Keep the better start by FINE-level cost (strictly, so pass 0 wins ties and an
         unchanged pair keeps its historical trajectory byte-for-byte). */
      if (levels_done == nlv && final_cost < CF_PENALTY && (!have_all || final_cost < cost_all)) {
          memcpy(best_all, best, sizeof best_all); cost_all = final_cost; have_all = 1;
      }
      if (O.verbose) fprintf(stderr, "[coreg fast] start %d final cost=%.5f (best so far %.5f)\n",
                             strat, final_cost, cost_all);
    }
    if (!sample_err && !ctx.opt_err && have_all) {
        memcpy(best, best_all, sizeof best); final_cost = cost_all; levels_done = nlv;
    }
    /* Final 9- and 12-DOF polish, for the WINNING start only. The multi-start needs every
       start carried to 2 mm (the 4 mm cost is NOT a safe arbiter: it disagrees with the 2 mm
       winner on 25/40 cases and would cost up to +0.064, losing M2017->SSW and both
       T1w1mm->stripped fits). The 2 mm 7-DOF cost, however, IS a safe arbiter -- measured over
       the same 40 cases, every disagreement with the fully-polished winner is worth at most
       +0.0026 in final cost, and all four cases that break 4 mm pruning agree. Since the
       9/12-DOF stages are ~76% of the 2 mm work and 2 mm is ~57% of runtime, ranking on the
       7-DOF result and polishing once recovers roughly half the multi-start's cost.
       ctx is still configured for the finest level (samples depend only on the FIXED level, so
       they are identical for every start) -- do not reorder this above the strategy loop. */
    if (!sample_err && !ctx.opt_err && levels_done == nlv && final_cost < CF_PENALTY) {
        if (O.max_dof >= 9)  cf_refine(&ctx, 9,  best, SEP[nlv-1]*0.35, SEP[nlv-1]*0.015, 250);
        if (O.max_dof >= 12) cf_refine(&ctx, 12, best, SEP[nlv-1]*0.25, SEP[nlv-1]*0.01,  350);
        final_cost = cf_cost_eval(best);
        if (O.verbose) fprintf(stderr, "[coreg fast] polished winner -> cost=%.5f (%d evals)\n",
                               final_cost, ctx.evals);
    }
    int opt_err = ctx.opt_err, total_evals = ctx.evals;
    cf_free_samples(&ctx);
    for (int i=0;i<3;i++){ cf_free_level(&Lf[i]); cf_free_level(&Lm[i]); cf_free_level(&Lwt[i]); }
    g_cf = NULL;
    powell_newuoa_free_threadlocal();  /* release the grow-only NEWUOA workspace */

    /* Failure is atomic: *result is written only after ALL requested levels completed,
       the optimizer reported no error, the fit has a valid (finite, in-overlap) cost,
       and the final affine passes validity. Any earlier failure leaves *result
       unchanged (the documented contract). A shortened DEBUG schedule (nlv<3) still
       requires all nlv levels to complete — it is distinct from a construction failure. */
    if (sample_err) {
        fprintf(stderr, "coreg fast: sample construction failed\n"); return 1; }
    if (levels_done < nlv) {
        fprintf(stderr, "coreg fast: level %d/%d failed (sample build)\n", levels_done, nlv); return 1; }
    if (opt_err) { fprintf(stderr, "coreg fast: optimizer error (allocation failure?)\n"); return 1; }
    if (!(final_cost < CF_PENALTY)) {   /* also rejects NaN (NaN < X is false) */
        fprintf(stderr, "coreg fast: no valid overlap / degenerate fit (cost=%.3g)\n", final_cost); return 1; }

    /* final world affine + validity */
    mat44 F2M = cf_affine(&ctx, best);
    if (!cf_affine_finite(F2M)) {
        fprintf(stderr, "coreg fast: non-finite affine\n"); return 1; }
    double det = (double)F2M.m[0][0]*(F2M.m[1][1]*F2M.m[2][2]-F2M.m[1][2]*F2M.m[2][1])
               - (double)F2M.m[0][1]*(F2M.m[1][0]*F2M.m[2][2]-F2M.m[1][2]*F2M.m[2][0])
               + (double)F2M.m[0][2]*(F2M.m[1][0]*F2M.m[2][1]-F2M.m[1][1]*F2M.m[2][0]);
    if (!(det > 0.0)) { fprintf(stderr, "coreg fast: non-positive determinant\n"); return 1; }
    /* guardrail check on final params */
    for (int i=0;i<12;i++) if (best[i] <= CF_LO[i]+1e-6 || best[i] >= CF_HI[i]-1e-6) {
        fprintf(stderr, "coreg fast: fit reached a parameter guardrail (dof %d)\n", i); return 1; }

    result->fixed_to_moving = F2M;
    result->final_cost = final_cost;
    result->evaluations = total_evals;
    result->levels_completed = levels_done;
    result->resolved_cost = O.cost;      /* validated config actually used */
    result->resolved_dof = ctx.dof_run;  /* highest DOF actually fitted (not the requested cap) */
    result->registration_ms = (cf_wtime() - t0) * 1000.0;
#ifdef AL_PROFILE
    fprintf(stderr, " [coreg profile] TOTAL estimate %.3f s: blur %.3f s, resample %.3f s, "
            "%d evals over %d level(s)  (final warp is applied separately by the caller)\n",
            result->registration_ms/1000.0, g_prof_blur, g_prof_resample,
            result->evaluations, result->levels_completed);
#endif
    return 0;
}
