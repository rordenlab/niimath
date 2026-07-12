/* Affine image registration engine - adapted from AFNI 3dAllineate (public domain)
   Original by RW Cox, National Institutes of Health.
   Ported to standalone NIfTI by the niimath project.

   This file implements affine (12 DOF) registration with four cost functions:
   lpc (signed local Pearson, cross-modal), lpa (absolute local Pearson),
   Hellinger (histogram-based, default), and global Pearson (ls, within-modality).
   Compile with -DAL_LPC_MICHO for lpc+ZZ/lpa+ZZ combined cost variant.
   Features twopass coarse-to-fine optimization, autoweight, source automask with
   noise fill (-source_automask), TOHD bloks, and configurable output interpolation.
   The Powell/NEWUOA optimizer is in a separate file. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <stdint.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif


/* Wall-clock timer: returns seconds as double */
static inline double al_wtime(void) {
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

static int al_mul_size(size_t a, size_t b, size_t *out)
{
    if (out == NULL) return 1;
    if (a != 0 && b > SIZE_MAX / a) return 1;
    *out = a * b;
    return 0;
}

/* Profiling support: compile with -DAL_PROFILE to enable timing output */
#ifdef AL_PROFILE
#define PROFILE_START(name) double _prof_##name = al_wtime()
#define PROFILE_END(name, label) fprintf(stderr, " [PROFILE] %-30s %.3f s\n", label, al_wtime() - _prof_##name)
#else
#define PROFILE_START(name) ((void)0)
#define PROFILE_END(name, label) ((void)0)
#endif

#include "nifti_io.h"
#include "allineate.h"
#include "coreg_fast.h"   /* fast SPM/FLIRT-inspired engine, shared by -allineate and -deface */

/* Thread-local storage qualifier for OpenMP safety */
#ifdef _OPENMP
#define AL_TLOCAL __thread
#else
#define AL_TLOCAL
#endif

/* Convert nifti_dmat44 (double) to mat44 (float) */
static mat44 dmat44_to_mat44(nifti_dmat44 d) {
    mat44 f;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            f.m[i][j] = (float)d.m[i][j];
    return f;
}

/* Convert mat44 (float) to nifti_dmat44 (double) */
static nifti_dmat44 mat44_to_dmat44(mat44 f) {
    nifti_dmat44 d;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            d.m[i][j] = (double)f.m[i][j];
    return d;
}

/*==========================================================================*/
/*============================== CONSTANTS ================================*/
/*==========================================================================*/

#define AL_BIGVAL   1.e+38f
#define AL_SMAGIC   208921148

/* Sparse sampling: fraction of task voxels for fine-pass matching (AFNI default) */
#define AL_SPARSE_SAMPLE_FRAC  0.47
/* Coarse-pass matching-point floor. MUST match AFNI's nmatch_setup (98765). An earlier
   distillation used 9999 (~10x fewer), which under-samples the coarse cost surface: on a
   cross-modal / offset case (T2w->avg152T1, header centroids ~25 mm apart) the noisy coarse
   Hellinger landed in a wrong rotation/translation basin from the header start, needing
   -cmass to recover (0.760 -> 0.840 L-R symmetry). Restoring AFNI's 98765 makes the coarse
   surface clean enough to find the right basin WITHOUT -cmass (0.833, matching AFNI's 0.826). */
#define AL_NPT_MATCH_MIN      98765

/* Cost function method codes (subset of AFNI's full list) */
#define GA_MATCH_PEARSON_SCALAR        1
#define GA_MATCH_HELLINGER_SCALAR      7
#define GA_MATCH_PEARSON_LOCALS       11  /* pure lpc (signed local Pearson) */
#define GA_MATCH_PEARSON_LOCALA       12  /* pure lpa (absolute local Pearson) */
#ifdef AL_LPC_MICHO
#define GA_MATCH_LPC_MICHO_SCALAR     13  /* lpc+ZZ (lpc + helper costs) */
#define GA_MATCH_LPA_MICHO_SCALAR     14  /* lpa+ZZ (lpa + helper costs) */
#endif

/* Test if method uses BLOK-based local correlation */
#ifdef AL_LPC_MICHO
#define METH_USES_BLOKS(m) ((m) == GA_MATCH_PEARSON_LOCALS    || \
                            (m) == GA_MATCH_PEARSON_LOCALA     || \
                            (m) == GA_MATCH_LPC_MICHO_SCALAR   || \
                            (m) == GA_MATCH_LPA_MICHO_SCALAR)
#define METH_IS_LPA(m) ((m) == GA_MATCH_LPA_MICHO_SCALAR || \
                         (m) == GA_MATCH_PEARSON_LOCALA)
#else
#define METH_USES_BLOKS(m) ((m) == GA_MATCH_PEARSON_LOCALS || \
                            (m) == GA_MATCH_PEARSON_LOCALA)
#define METH_IS_LPA(m) ((m) == GA_MATCH_PEARSON_LOCALA)
#endif

#define DEFAULT_TBEST_LPA 17

/* Interpolation codes */
/* Interpolation constants are defined in allineate.h (AL_INTERP_*) */

/* Smooth codes */
#define GA_SMOOTH_GAUSSIAN 1

/* BLOK types */
#define GA_BLOK_RHDD  1
#define GA_BLOK_TOHD  2

/* Matrix ordering */
#define MATORDER_SDU  1
#define MATORDER_SUD  2
#define MATORDER_DSU  3
#define MATORDER_DUS  4
#define MATORDER_USD  5
#define MATORDER_UDS  6
#define SMAT_LOWER    1
#define SMAT_UPPER    2
#define DELTA_BEFORE  1
#define DELTA_AFTER   2

#define FWHM_TO_SIGMA(f) ((f)/2.3548f)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define D2R ((float)(M_PI/180.0))

#define CMAX 0.9999f /* max correlation magnitude to allow */

#define NPER 262144  /* max points to warp at once */

#define PARAM_MAXTRIAL 29 /* AFNI maximum number of trial parameter sets */

/*==========================================================================*/
/*========================== TYPE DEFINITIONS ==============================*/
/*==========================================================================*/

/* BLOK set for local Pearson correlation */
typedef struct {
    int    num;       /* number of bloks */
    int   *nelm;      /* elements per blok */
    int  **elm;       /* element index lists */
    int    nx, ny, nz;
    float  dx, dy, dz;
    float  ppow;
} GA_BLOK_set;

/* Per-parameter descriptor */
typedef struct {
    float min, max, siz;
    float delta, toler;
    float ident;
    float val_init, val_pinit, val_fixed, val_out;
    float val_trial[PARAM_MAXTRIAL+2];
    int   idx_trial[PARAM_MAXTRIAL+2];
    int   fixed;       /* 0=free, 1=temp-fixed, 2=perm-fixed */
    char  name[32];
} GA_param;

/* Warp function pointer */
typedef void (*GA_warpfunc)(int npar, float *wpar,
                            int npt, float *xi, float *yi, float *zi,
                                     float *xo, float *yo, float *zo);

/* Master alignment setup structure */
typedef struct {
    int setup;
    int match_code;
    int interp_code;
    int smooth_code;
    float smooth_radius_base, smooth_radius_targ;

    /* Base image */
    float *bsim;           /* base image data */
    float *bsims;          /* smoothed base */
    int bnx, bny, bnz;    /* base dimensions */
    float bdx, bdy, bdz;  /* base voxel sizes */
    float bsbot, bstop, bsclip;
    float bs_topclip;      /* AFNI mri_topclip: histogram-MEMBERSHIP top (drop above); distinct
                              from bsclip = clipate top used for edge-binning */
    unsigned char *bmask;  /* base mask */
    float *bwght;          /* base weights */
    int nmask;

    /* Source/target image */
    float *ajim;           /* source image data */
    float *ajims;          /* smoothed source */
    int anx, any, anz;     /* source dimensions */
    float adx, ady, adz;   /* source voxel sizes */
    float ajbot, ajtop, ajclip;
    float aj_topclip;      /* AFNI mri_topclip: histogram-MEMBERSHIP top (drop above); distinct
                              from ajclip = clipate top used for edge-binning */
    float ajmin, ajmax;    /* original data range (for cubic overshoot clamping) */
    /* AFNI set_2Dhist_xyclip: edge-bin clips from the SAMPLE distribution, recomputed per
       stage. hxc=source(x), hyc=base(y). need_hist_setup is set at each stage's setup and
       consumed (compute + clear) by the stage's FIRST, sequential cost eval — never in the
       parallel region — so the result is deterministic (mirrors AFNI GA_scalar_fitter). */
    float hxc_bot, hxc_top, hyc_bot, hyc_top;
    int   need_hist_setup;
    unsigned char *ajmask;  /* source mask */
    int najmask;

    /* Coordinate transforms */
    mat44 base_cmat, base_imat;  /* base index->xyz and xyz->index */
    mat44 targ_cmat, targ_imat;  /* source index->xyz and xyz->index */
    float base_di, base_dj, base_dk;
    float targ_di, targ_dj, targ_dk;

    /* Control points */
    int npt_match;
    float *im_ar, *jm_ar, *km_ar;  /* base space control point indices */
    float *bvm;                      /* base values at control points */
    float *wvm;                      /* weight values at control points */

    /* Warp parameters */
    GA_warpfunc wfunc;
    int wfunc_numpar, wfunc_numfree;
    GA_param *wfunc_param;
    int wfunc_ntrial;

    /* BLOK stuff */
    GA_BLOK_set *blokset;
    int bloktype, blokmin;
    float blokrad;

    /* Best cost */
    float vbest;

#ifdef AL_LPC_MICHO
    /* micho (lpc+ZZ) parameters */
    double micho_mi, micho_nmi, micho_crA, micho_hel, micho_ov;
    int micho_zfinal;
#endif

    /* source automask noise fill (AFNI's -source_automask) */
    int ajmask_ranfill;       /* if nonzero, fill outside ajmask with noise */
    float aj_ubot, aj_usiz;   /* noise range: ubot + usiz*(u1+u2) */
    float *ajim_orig;          /* backup of original source data before noise fill */

    /* dark (image-minimum) automask (-dark_automask): drop a matched pair from the
       cost when the base OR warped-source value sits at that image's darkest value
       (background / zero-pad / NaN-filled-to-zero). dark_base/dark_targ are those
       per-image minima; applied per-eval in GA_scalar_fitter. */
    int do_dark_automask;
    float dark_base, dark_targ;
} GA_setup;

/* For qsort_floatint */
typedef struct { float a; int b; } float_int;

/*==========================================================================*/
/*====================== EXTERNAL OPTIMIZER DECLARATIONS ==================*/
/*==========================================================================*/

extern int powell_newuoa(int ndim, double *x, double rstart, double rend,
                         int maxcall, double (*ufunc)(int, double *));
extern int powell_newuoa_con(int ndim, double *x, double *xbot, double *xtop,
                             int nrand, double rstart, double rend,
                             int maxcall, double (*ufunc)(int, double *));
extern void powell_set_mfac(float mm, float aa);
extern void powell_get_mfac(float *mm, float *aa);
extern void powell_newuoa_free_threadlocal(void);

/*==========================================================================*/
/*========================== GLOBAL STATE =================================*/
/*==========================================================================*/

static GA_setup *gstup = NULL;  /* current setup for optimizer callback */

static int aff_use_before = 0, aff_use_after = 0;
static mat44 aff_before, aff_after;

/* World-space FIXED(base)->MOVING(source) affine from the most recent successful
   nii_allineate() fit, exposed via nii_last_affine() so the CLI can save it
   (-savemat) without changing nii_allineate's signature or making it dereference a
   new al_opts field (which would break the niimath drop-in, where nii_allineate must
   ignore demo-only fields). Serial-only, like the other engine globals. */
static mat44 g_last_affine;
static int g_last_affine_valid = 0;

/* Wall-clock ms spent in the coarse and fine passes of the most recent al_register()
   (always measured, unlike the -DAL_PROFILE breakdown). nii_allineate() reports them
   in its "Registration completed" line. Serial-only, like the other engine globals. */
static double g_last_coarse_ms = 0.0, g_last_fine_ms = 0.0;

/* When >0, an upper bound (in mm) on each |x/y/z-shift| parameter range in
   al_register. Set transiently by nii_symmetry() around a mirror registration to
   keep the recovered translation near isocenter, then reset to 0. Single-threaded
   at the point of use, so touching this file-scope global is safe (mirrors the
   aff_use_before/after pattern). */
static float al_shift_max_override = 0.0f;

/* When nonzero, GA_setup_affine ties the y- and z-scale to the x-scale (parvec[6])
   so a single free scale parameter drives a GLOBAL ISOTROPIC zoom (no shape
   distortion), and al_register widens param[6]'s range. Set transiently by
   nii_sagseed() around the -zoom seed fit + its transform extraction, then reset to
   0 — so the isotropic *tie* never affects the main nii_allineate() fit. (The main
   fit's scale *range* is still widened separately, via al_register's `relax_scale`
   arg, which -zoom also sets; all other regularization is unchanged.) Read during
   the parallel cost evals but never written there (mirrors aff_use_*). */
static int al_zoom_isotropic = 0;

/* Thread-local workspace for cost function evaluation (avoids per-call malloc) */
static AL_TLOCAL float *tl_avm = NULL;
static AL_TLOCAL int    tl_avm_len = 0;
static AL_TLOCAL float *tl_weff = NULL;   /* -dark_automask per-eval effective weights */
static AL_TLOCAL int    tl_weff_len = 0;
static AL_TLOCAL float *tl_wpar = NULL;
static AL_TLOCAL int    tl_wpar_len = 0;
static AL_TLOCAL float *tl_wbuf = NULL;  /* holds imw,jmw,kmw contiguously */
static AL_TLOCAL int    tl_wbuf_len = 0;

/* These are fixed for our minimal affine-only use case */
#define AL_MATORDER MATORDER_SDU
#define AL_SMAT     SMAT_LOWER
#define AL_DCODE    DELTA_AFTER
#define AL_OUTVAL   0.0f

/*==========================================================================*/
/*========================= PRNG (from AFNI) ==============================*/
/*==========================================================================*/

static unsigned long long MYa = 62003, MYb = 15485863, MYx = 15482917;

static float myunif(void)
{
    MYx = MYa * MYx + MYb;
    return ((unsigned int)MYx) / 4294967296.0f;
}

static void myunif_reset(unsigned long long x) { MYx = x; }

/*==========================================================================*/
/*================ PRED01: periodic reduction to [0,1] ====================*/
/*==========================================================================*/

#define PRED01(x) fabsf((float)(x) - 2.0f*floorf(0.5f*((float)(x)+1.0f)))

/*==========================================================================*/
/*====================== SECTION 3: UTILITY FUNCTIONS =====================*/
/*==========================================================================*/

static int floatint_compare(const void *a, const void *b)
{
    float fa = ((const float_int *)a)->a;
    float fb = ((const float_int *)b)->a;
    if (fa < fb) return -1;
    if (fa > fb) return  1;
    return 0;
}

static void qsort_floatint(int n, float *a, int *b)
{
    float_int *fi;
    int i;
    if (n < 1 || a == NULL) return;
    fi = (float_int *)malloc(sizeof(float_int) * n);
    if (fi == NULL) return;
    for (i = 0; i < n; i++) { fi[i].a = a[i]; fi[i].b = (b != NULL) ? b[i] : i; }
    qsort(fi, n, sizeof(float_int), floatint_compare);
    for (i = 0; i < n; i++) { a[i] = fi[i].a; if (b != NULL) b[i] = fi[i].b; }
    free(fi);
}

/* In-place selection of the k-th smallest element (0-based) via 3-way quickselect.
 * O(n) average; the 3-way (Dutch-flag) partition stays O(n) even when many values
 * are equal (image background). Reorders `a`. Returns a[k]. This replaces the
 * libc qsort+comparator the percentile helpers used to call — a comparator is a
 * per-comparison indirect call (call_indirect), a severe penalty under WASM. */
static float al_select_rank(float *a, int n, int k)
{
    if (n < 1) return 0.0f;
    if (k < 0) k = 0; else if (k >= n) k = n - 1;
    int lo = 0, hi = n - 1;
    while (lo < hi) {
        float pivot = a[lo + ((hi - lo) >> 1)];
        int lt = lo, gt = hi, i = lo;
        while (i <= gt) {
            if (a[i] < pivot)      { float t = a[lt]; a[lt++] = a[i]; a[i++] = t; }
            else if (a[i] > pivot) { float t = a[gt]; a[gt--] = a[i]; a[i]   = t; }
            else i++;
        }
        if (k < lt) hi = lt - 1;
        else if (k > gt) lo = gt + 1;
        else break;   /* k in the ==pivot band -> a[k] == pivot */
    }
    return a[k];
}

/* Column norm of mat44 */
static float mat44_colnorm(mat44 m, int col)
{
    float s = 0.0f;
    int i;
    for (i = 0; i < 3; i++) s += m.m[i][col] * m.m[i][col];
    return sqrtf(s);
}

/* Apply mat44 to a vector */
static inline void mat44_vec(mat44 m, float x, float y, float z,
                             float *ox, float *oy, float *oz)
{
    *ox = m.m[0][0]*x + m.m[0][1]*y + m.m[0][2]*z + m.m[0][3];
    *oy = m.m[1][0]*x + m.m[1][1]*y + m.m[1][2]*z + m.m[1][3];
    *oz = m.m[2][0]*x + m.m[2][1]*y + m.m[2][2]*z + m.m[2][3];
}

/* Macro version for tight loops — avoids pass-by-value overhead.
   Affine: row 3 is always [0 0 0 1], so only 3 rows computed. */
#define MAT44_VEC(A,x,y,z,a,b,c) \
  ( (a) = (A).m[0][0]*(x) + (A).m[0][1]*(y) + (A).m[0][2]*(z) + (A).m[0][3] , \
    (b) = (A).m[1][0]*(x) + (A).m[1][1]*(y) + (A).m[1][2]*(z) + (A).m[1][3] , \
    (c) = (A).m[2][0]*(x) + (A).m[2][1]*(y) + (A).m[2][2]*(z) + (A).m[2][3]  )

/* Load a diagonal mat44 */
static mat44 mat44_diag(float a, float b, float c)
{
    mat44 m;
    memset(&m, 0, sizeof(mat44));
    m.m[0][0] = a; m.m[1][1] = b; m.m[2][2] = c; m.m[3][3] = 1.0f;
    return m;
}

/* Index->world (mm) matrix for a pixdim-centered frame (no sform/qform): a
 * diagonal scale by |pixdim| with the volume centered on the origin. Shared by
 * al_register's sform-absent fallback and nii_symmetry (single source of truth
 * for the centered-frame convention). Non-positive pixdims clamp to 1.0. */
static mat44 al_pixdim_frame(int nx, int ny, int nz, float dx, float dy, float dz)
{
    if (dx <= 0.0f) dx = 1.0f;
    if (dy <= 0.0f) dy = 1.0f;
    if (dz <= 0.0f) dz = 1.0f;
    mat44 m = mat44_diag(dx, dy, dz);
    m.m[0][3] = -(nx - 1) * 0.5f * dx;
    m.m[1][3] = -(ny - 1) * 0.5f * dy;
    m.m[2][3] = -(nz - 1) * 0.5f * dz;
    return m;
}

/* Finite (non-NaN, non-Inf) test via a magnitude guard: this TU is built with
 * -ffast-math, under which isfinite()/isnan() are unreliable, but the ordered
 * comparison `-FLT_MAX <= v <= FLT_MAX` still excludes NaN (compares false) and
 * ±Inf. Single source of truth for the several places that filter non-finite. */
static inline int al_finitef(float v) { return v >= -FLT_MAX && v <= FLT_MAX; }

/* Return nonzero if m is usable as an index->world transform: every entry is
 * finite (see al_finitef — isfinite() is unreliable under -ffast-math) and the
 * upper-left 3x3 is non-singular. The singularity test is
 * scale-invariant: |det| divided by the three column norms is the product of the
 * sines of the inter-axis angles, which collapses to ~0 for an all-zero or
 * coplanar ("bogus") matrix regardless of voxel size. */
static int al_mat44_usable(mat44 m)
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 4; j++)
            if (!al_finitef(m.m[i][j])) return 0;
    double c[3];
    for (int a = 0; a < 3; a++) {
        c[a] = sqrt((double)m.m[0][a]*m.m[0][a] +
                    (double)m.m[1][a]*m.m[1][a] +
                    (double)m.m[2][a]*m.m[2][a]);
        if (!(c[a] > 0.0)) return 0;
    }
    double det =
        (double)m.m[0][0] * ((double)m.m[1][1]*m.m[2][2] - (double)m.m[1][2]*m.m[2][1])
      - (double)m.m[0][1] * ((double)m.m[1][0]*m.m[2][2] - (double)m.m[1][2]*m.m[2][0])
      + (double)m.m[0][2] * ((double)m.m[1][0]*m.m[2][1] - (double)m.m[1][1]*m.m[2][0]);
    return fabs(det) / (c[0]*c[1]*c[2]) > 1e-4;
}

/* Choose an image's index->world (mm) transform with NIfTI precedence: prefer
 * the sform when its code is >= the qform's, otherwise the qform; if the
 * preferred form is unset or degenerate ("bogus"), fall back to whichever form
 * is usable. Writes *out and returns 0 on success; returns 1 when neither the
 * sform nor the qform yields a usable transform (caller should error out). */
int al_image_xform(const nifti_image *nim, mat44 *out)
{
    mat44 s, q;
    int shave = (nim->sform_code > 0);
    int qhave = (nim->qform_code > 0);
    if (shave) { s = dmat44_to_mat44(nim->sto_xyz); shave = al_mat44_usable(s); }
    if (qhave) { q = dmat44_to_mat44(nim->qto_xyz); qhave = al_mat44_usable(q); }
    if (shave && (!qhave || nim->sform_code >= nim->qform_code)) { *out = s; return 0; }
    if (qhave) { *out = q; return 0; }
    return 1;
}

/* Index->world (mm) transform with the SINGLE no-form fallback policy: al_image_xform's
   coded sform/qform selection when a usable form exists, else a pixdim-centered frame
   (al_pixdim_frame). This is the one policy used by every registration/geometry entry
   point — al_register (base/source), the fast engine (coreg_fast), -sym, -com, -sagseed,
   and apply (nii_apply_affine) — so a valid NIfTI/ANALYZE with both form codes 0 always
   yields a usable frame instead of some callers erroring and others falling back. Never
   fails (a centered frame is always constructible). Pass `who` (non-NULL) to log the
   fallback; NULL for silent. Exposed non-static so coreg_fast.c shares the exact policy. */
void al_image_xform_or_pixdim(const nifti_image *nim, mat44 *out, const char *who)
{
    if (al_image_xform(nim, out)) {
        *out = al_pixdim_frame(nim->nx, nim->ny, nim->nz,
                               (float)fabs(nim->dx), (float)fabs(nim->dy), (float)fabs(nim->dz));
        if (who)
            fprintf(stderr, " + [%s] no valid sform/qform: using pixdim-centered frame\n", who);
    }
}

/* GCD */
static int ga_gcd(int m, int n)
{
    while (m > 0) {
        if (n > m) { int t = m; m = n; n = t; }
        m -= n;
    }
    return n;
}

/* Find number relatively prime to n */
static int ga_find_relprime_fixed(int n)
{
    int dj, n5 = n / 5;
    if (n5 < 2) return 1;
    for (dj = n5; ga_gcd(n, dj) > 1; dj++) ;
    return dj;
}

/*==========================================================================*/
/*==================== SECTION 4: GAUSSIAN BLUR ===========================*/
/*==========================================================================*/

/* Separable 3D Gaussian blur in-place.
   sigma is in the same units as dx/dy/dz (i.e., in voxel-size units). */
static int gaussian_blur_3d(float *data, int nx, int ny, int nz,
                            float dx, float dy, float dz, float sigma)
{
    int krad, nxy, ii, jj, kk, rr;
    float *kernel, *buf;
    float sum, sigma_v, ksum;

    if (data == NULL || sigma <= 0.0f) return 0;
    nxy = nx * ny;

    /* Blur along X */
    sigma_v = sigma / dx;
    if (sigma_v > 0.25f) {
        krad = (int)(4.0f * sigma_v + 0.5f);
        if (krad < 1) krad = 1;
        kernel = (float *)calloc(2 * krad + 1, sizeof(float));
        buf    = (float *)calloc(nx, sizeof(float));
        if (!kernel || !buf) { free(kernel); free(buf); return 1; }
        ksum = 0.0f;
        for (rr = -krad; rr <= krad; rr++) {
            float v = expf(-0.5f * (rr * rr) / (sigma_v * sigma_v));
            kernel[rr + krad] = v; ksum += v;
        }
        for (rr = 0; rr < 2 * krad + 1; rr++) kernel[rr] /= ksum;

        for (kk = 0; kk < nz; kk++) {
            for (jj = 0; jj < ny; jj++) {
                float *row = data + jj * nx + kk * nxy;
                memcpy(buf, row, sizeof(float) * nx);
                for (ii = 0; ii < nx; ii++) {
                    sum = 0.0f;
                    for (rr = -krad; rr <= krad; rr++) {
                        int idx = ii + rr;
                        if (idx < 0) idx = 0; else if (idx >= nx) idx = nx - 1;
                        sum += buf[idx] * kernel[rr + krad];
                    }
                    row[ii] = sum;
                }
            }
        }
        free(kernel); free(buf);
    }

    /* Blur along Y */
    sigma_v = sigma / dy;
    if (sigma_v > 0.25f) {
        krad = (int)(4.0f * sigma_v + 0.5f);
        if (krad < 1) krad = 1;
        kernel = (float *)calloc(2 * krad + 1, sizeof(float));
        buf    = (float *)calloc(ny, sizeof(float));
        if (!kernel || !buf) { free(kernel); free(buf); return 1; }
        ksum = 0.0f;
        for (rr = -krad; rr <= krad; rr++) {
            float v = expf(-0.5f * (rr * rr) / (sigma_v * sigma_v));
            kernel[rr + krad] = v; ksum += v;
        }
        for (rr = 0; rr < 2 * krad + 1; rr++) kernel[rr] /= ksum;

        for (kk = 0; kk < nz; kk++) {
            for (ii = 0; ii < nx; ii++) {
                for (jj = 0; jj < ny; jj++)
                    buf[jj] = data[ii + jj * nx + kk * nxy];
                for (jj = 0; jj < ny; jj++) {
                    sum = 0.0f;
                    for (rr = -krad; rr <= krad; rr++) {
                        int idx = jj + rr;
                        if (idx < 0) idx = 0; else if (idx >= ny) idx = ny - 1;
                        sum += buf[idx] * kernel[rr + krad];
                    }
                    data[ii + jj * nx + kk * nxy] = sum;
                }
            }
        }
        free(kernel); free(buf);
    }

    /* Blur along Z */
    sigma_v = sigma / dz;
    if (sigma_v > 0.25f && nz > 1) {
        krad = (int)(4.0f * sigma_v + 0.5f);
        if (krad < 1) krad = 1;
        kernel = (float *)calloc(2 * krad + 1, sizeof(float));
        buf    = (float *)calloc(nz, sizeof(float));
        if (!kernel || !buf) { free(kernel); free(buf); return 1; }
        ksum = 0.0f;
        for (rr = -krad; rr <= krad; rr++) {
            float v = expf(-0.5f * (rr * rr) / (sigma_v * sigma_v));
            kernel[rr + krad] = v; ksum += v;
        }
        for (rr = 0; rr < 2 * krad + 1; rr++) kernel[rr] /= ksum;

        for (jj = 0; jj < ny; jj++) {
            for (ii = 0; ii < nx; ii++) {
                for (kk = 0; kk < nz; kk++)
                    buf[kk] = data[ii + jj * nx + kk * nxy];
                for (kk = 0; kk < nz; kk++) {
                    sum = 0.0f;
                    for (rr = -krad; rr <= krad; rr++) {
                        int idx = kk + rr;
                        if (idx < 0) idx = 0; else if (idx >= nz) idx = nz - 1;
                        sum += buf[idx] * kernel[rr + krad];
                    }
                    data[ii + jj * nx + kk * nxy] = sum;
                }
            }
        }
        free(kernel); free(buf);
    }
    return 0;
}

/* Smooth a float image (FWHM-based radius like AFNI's GA_smooth).
   Allocates and returns a smoothed copy. Caller must free. */
static float *al_smooth(float *im, int nx, int ny, int nz,
                        float dx, float dy, float dz, float fwhm)
{
    float *om;
    float sigma;
    int nvox = nx * ny * nz;

    if (im == NULL || fwhm <= 0.0f) return NULL;
    om = (float *)malloc(sizeof(float) * nvox);
    if (om == NULL) return NULL;
    memcpy(om, im, sizeof(float) * nvox);
    sigma = FWHM_TO_SIGMA(fwhm);
    if (gaussian_blur_3d(om, nx, ny, nz, dx, dy, dz, sigma)) {
        free(om);
        return NULL;
    }
    return om;
}

/*==========================================================================*/
/*==================== SECTION 5: INTERPOLATION FUNCTIONS =================*/
/*==========================================================================*/

#undef  FAR
#define FAR(i,j,k)  far[(i)+(j)*nx+(k)*nxy]

#undef  CLIP
#define CLIP(mm,nn) if(mm < 0)mm=0; else if(mm > nn)mm=nn

#undef  ISTINY
#define ISTINY(a) (fabsf(a) < 0.0001f)

/* Nearest-neighbor interpolation */
static void GA_interp_NN(float *far, int nx, int ny, int nz,
                         int npp, float *ip, float *jp, float *kp, float *vv)
{
    int nxy = nx * ny, pp, ii, jj, kk;
    float nxh = nx - 0.501f, nyh = ny - 0.501f, nzh = nz - 0.501f;
    float xx, yy, zz;

    for (pp = 0; pp < npp; pp++) {
        xx = ip[pp]; if (xx < -0.499f || xx > nxh) { vv[pp] = AL_OUTVAL; continue; }
        yy = jp[pp]; if (yy < -0.499f || yy > nyh) { vv[pp] = AL_OUTVAL; continue; }
        zz = kp[pp]; if (zz < -0.499f || zz > nzh) { vv[pp] = AL_OUTVAL; continue; }
        ii = (int)(xx + 0.5f); jj = (int)(yy + 0.5f); kk = (int)(zz + 0.5f);
        vv[pp] = FAR(ii, jj, kk);
    }
}

/* Linear interpolation */
static void GA_interp_linear(float *far, int nx, int ny, int nz,
                             int npp, float *ip, float *jp, float *kp, float *vv)
{
    int nxy = nx * ny, pp;
    float nxh = nx - 0.501f, nyh = ny - 0.501f, nzh = nz - 0.501f;
    float xx, yy, zz, fx, fy, fz;
    int nx1 = nx - 1, ny1 = ny - 1, nz1 = nz - 1;
    int ix_00, ix_p1, jy_00, jy_p1, kz_00, kz_p1;
    float wt_00, wt_p1;
    float f_j00_k00, f_jp1_k00, f_j00_kp1, f_jp1_kp1, f_k00, f_kp1;
    float ix, jy, kz;

    for (pp = 0; pp < npp; pp++) {
        xx = ip[pp]; if (xx < -0.499f || xx > nxh) { vv[pp] = AL_OUTVAL; continue; }
        yy = jp[pp]; if (yy < -0.499f || yy > nyh) { vv[pp] = AL_OUTVAL; continue; }
        zz = kp[pp]; if (zz < -0.499f || zz > nzh) { vv[pp] = AL_OUTVAL; continue; }

        ix = floorf(xx); fx = xx - ix;
        jy = floorf(yy); fy = yy - jy;
        kz = floorf(zz); fz = zz - kz;

        ix_00 = (int)ix; ix_p1 = ix_00 + 1; CLIP(ix_00, nx1); CLIP(ix_p1, nx1);
        jy_00 = (int)jy; jy_p1 = jy_00 + 1; CLIP(jy_00, ny1); CLIP(jy_p1, ny1);
        kz_00 = (int)kz; kz_p1 = kz_00 + 1; CLIP(kz_00, nz1); CLIP(kz_p1, nz1);

        wt_00 = 1.0f - fx; wt_p1 = fx;

#undef  XINT
#define XINT(j,k) wt_00*FAR(ix_00,j,k)+wt_p1*FAR(ix_p1,j,k)

        f_j00_k00 = XINT(jy_00, kz_00); f_jp1_k00 = XINT(jy_p1, kz_00);
        f_j00_kp1 = XINT(jy_00, kz_p1); f_jp1_kp1 = XINT(jy_p1, kz_p1);

        wt_00 = 1.0f - fy; wt_p1 = fy;
        f_k00 = wt_00 * f_j00_k00 + wt_p1 * f_jp1_k00;
        f_kp1 = wt_00 * f_j00_kp1 + wt_p1 * f_jp1_kp1;

        vv[pp] = (1.0f - fz) * f_k00 + fz * f_kp1;
    }
}

/* Lagrange cubic interpolation polynomials */
#undef  P_M1
#undef  P_00
#undef  P_P1
#undef  P_P2
#undef  P_FACTOR
#define P_M1(x)  (-(x)*((x)-1)*((x)-2))
#define P_00(x)  (3*((x)+1)*((x)-1)*((x)-2))
#define P_P1(x)  (-3*(x)*((x)+1)*((x)-2))
#define P_P2(x)  ((x)*((x)+1)*((x)-1))
#define P_FACTOR 4.62962963e-3f   /* 1/216 */

/* Cubic interpolation */
static void GA_interp_cubic(float *far, int nx, int ny, int nz,
                            int npp, float *ip, float *jp, float *kp, float *vv)
{
    int nxy = nx * ny, pp;
    float nxh = nx - 0.501f, nyh = ny - 0.501f, nzh = nz - 0.501f;
    float xx, yy, zz, fx, fy, fz;
    int nx1 = nx - 1, ny1 = ny - 1, nz1 = nz - 1;
    int ix, jy, kz;
    int ix_m1, ix_00, ix_p1, ix_p2;
    int jy_m1, jy_00, jy_p1, jy_p2;
    int kz_m1, kz_00, kz_p1, kz_p2;
    float wt_m1, wt_00, wt_p1, wt_p2;
    float f_jm1_km1, f_j00_km1, f_jp1_km1, f_jp2_km1,
          f_jm1_k00, f_j00_k00, f_jp1_k00, f_jp2_k00,
          f_jm1_kp1, f_j00_kp1, f_jp1_kp1, f_jp2_kp1,
          f_jm1_kp2, f_j00_kp2, f_jp1_kp2, f_jp2_kp2,
          f_km1,     f_k00,     f_kp1,     f_kp2;

    for (pp = 0; pp < npp; pp++) {
        xx = ip[pp]; if (xx < -0.499f || xx > nxh) { vv[pp] = AL_OUTVAL; continue; }
        yy = jp[pp]; if (yy < -0.499f || yy > nyh) { vv[pp] = AL_OUTVAL; continue; }
        zz = kp[pp]; if (zz < -0.499f || zz > nzh) { vv[pp] = AL_OUTVAL; continue; }

        ix = (int)floorf(xx); fx = xx - ix;
        jy = (int)floorf(yy); fy = yy - jy;
        kz = (int)floorf(zz); fz = zz - kz;

        if (ISTINY(fx) && ISTINY(fy) && ISTINY(fz)) {
            CLIP(ix, nx1); CLIP(jy, ny1); CLIP(kz, nz1);
            vv[pp] = FAR(ix, jy, kz); continue;
        }

        ix_m1 = ix - 1;    ix_00 = ix;        ix_p1 = ix + 1;    ix_p2 = ix + 2;
        CLIP(ix_m1, nx1); CLIP(ix_00, nx1); CLIP(ix_p1, nx1); CLIP(ix_p2, nx1);

        jy_m1 = jy - 1;    jy_00 = jy;        jy_p1 = jy + 1;    jy_p2 = jy + 2;
        CLIP(jy_m1, ny1); CLIP(jy_00, ny1); CLIP(jy_p1, ny1); CLIP(jy_p2, ny1);

        kz_m1 = kz - 1;    kz_00 = kz;        kz_p1 = kz + 1;    kz_p2 = kz + 2;
        CLIP(kz_m1, nz1); CLIP(kz_00, nz1); CLIP(kz_p1, nz1); CLIP(kz_p2, nz1);

        wt_m1 = P_M1(fx); wt_00 = P_00(fx); wt_p1 = P_P1(fx); wt_p2 = P_P2(fx);

#undef  XINT
#define XINT(j,k) wt_m1*FAR(ix_m1,j,k)+wt_00*FAR(ix_00,j,k)\
                 +wt_p1*FAR(ix_p1,j,k)+wt_p2*FAR(ix_p2,j,k)

        f_jm1_km1 = XINT(jy_m1, kz_m1); f_j00_km1 = XINT(jy_00, kz_m1);
        f_jp1_km1 = XINT(jy_p1, kz_m1); f_jp2_km1 = XINT(jy_p2, kz_m1);
        f_jm1_k00 = XINT(jy_m1, kz_00); f_j00_k00 = XINT(jy_00, kz_00);
        f_jp1_k00 = XINT(jy_p1, kz_00); f_jp2_k00 = XINT(jy_p2, kz_00);
        f_jm1_kp1 = XINT(jy_m1, kz_p1); f_j00_kp1 = XINT(jy_00, kz_p1);
        f_jp1_kp1 = XINT(jy_p1, kz_p1); f_jp2_kp1 = XINT(jy_p2, kz_p1);
        f_jm1_kp2 = XINT(jy_m1, kz_p2); f_j00_kp2 = XINT(jy_00, kz_p2);
        f_jp1_kp2 = XINT(jy_p1, kz_p2); f_jp2_kp2 = XINT(jy_p2, kz_p2);

        wt_m1 = P_M1(fy); wt_00 = P_00(fy); wt_p1 = P_P1(fy); wt_p2 = P_P2(fy);

        f_km1 = wt_m1 * f_jm1_km1 + wt_00 * f_j00_km1
               + wt_p1 * f_jp1_km1 + wt_p2 * f_jp2_km1;
        f_k00 = wt_m1 * f_jm1_k00 + wt_00 * f_j00_k00
               + wt_p1 * f_jp1_k00 + wt_p2 * f_jp2_k00;
        f_kp1 = wt_m1 * f_jm1_kp1 + wt_00 * f_j00_kp1
               + wt_p1 * f_jp1_kp1 + wt_p2 * f_jp2_kp1;
        f_kp2 = wt_m1 * f_jm1_kp2 + wt_00 * f_j00_kp2
               + wt_p1 * f_jp1_kp2 + wt_p2 * f_jp2_kp2;

        wt_m1 = P_M1(fz); wt_00 = P_00(fz); wt_p1 = P_P1(fz); wt_p2 = P_P2(fz);

        vv[pp] = P_FACTOR * (wt_m1 * f_km1 + wt_00 * f_k00
                            + wt_p1 * f_kp1 + wt_p2 * f_kp2);
    }
}

/*==========================================================================*/
/*==================== SECTION 6: BLOK SET CREATION =======================*/
/*==========================================================================*/

/* Inside-RHDD test: |a+b|<=siz AND |a-b|<=siz AND |a+c|<=siz
   AND |a-c|<=siz AND |b+c|<=siz AND |b-c|<=siz */
#define FAS(a,s)  ((a) <= (s) && (a) >= -(s))

static int blok_inside_rhdd(float a, float b, float c, float siz)
{
    return FAS(a + b, siz) && FAS(a - b, siz) &&
           FAS(a + c, siz) && FAS(a - c, siz) &&
           FAS(b + c, siz) && FAS(b - c, siz);
}

static int blok_inside_tohd(float a, float b, float c, float siz)
{
    return FAS(a, siz) && FAS(b, siz) && FAS(c, siz) &&
           FAS(a + b + c, 1.5f * siz) && FAS(a - b + c, 1.5f * siz) &&
           FAS(a + b - c, 1.5f * siz) && FAS(a - b - c, 1.5f * siz);
}

static int blok_inside(int bt, float a, float b, float c, float siz)
{
    if (bt == GA_BLOK_RHDD) return blok_inside_rhdd(a, b, c, siz);
    if (bt == GA_BLOK_TOHD) return blok_inside_tohd(a, b, c, siz);
    return 0;
}

/* 3x3 matrix operations for BLOK lattice */
typedef struct { float m[3][3]; } mat33f;

static mat33f mat33_inv(mat33f a)
{
    mat33f r;
    float det = a.m[0][0] * (a.m[1][1]*a.m[2][2] - a.m[1][2]*a.m[2][1])
              - a.m[0][1] * (a.m[1][0]*a.m[2][2] - a.m[1][2]*a.m[2][0])
              + a.m[0][2] * (a.m[1][0]*a.m[2][1] - a.m[1][1]*a.m[2][0]);
    float idet;
    if (fabsf(det) < 1.e-30f) det = 1.e-30f;
    idet = 1.0f / det;
    r.m[0][0] =  (a.m[1][1]*a.m[2][2]-a.m[1][2]*a.m[2][1]) * idet;
    r.m[0][1] = -(a.m[0][1]*a.m[2][2]-a.m[0][2]*a.m[2][1]) * idet;
    r.m[0][2] =  (a.m[0][1]*a.m[1][2]-a.m[0][2]*a.m[1][1]) * idet;
    r.m[1][0] = -(a.m[1][0]*a.m[2][2]-a.m[1][2]*a.m[2][0]) * idet;
    r.m[1][1] =  (a.m[0][0]*a.m[2][2]-a.m[0][2]*a.m[2][0]) * idet;
    r.m[1][2] = -(a.m[0][0]*a.m[1][2]-a.m[0][2]*a.m[1][0]) * idet;
    r.m[2][0] =  (a.m[1][0]*a.m[2][1]-a.m[1][1]*a.m[2][0]) * idet;
    r.m[2][1] = -(a.m[0][0]*a.m[2][1]-a.m[0][1]*a.m[2][0]) * idet;
    r.m[2][2] =  (a.m[0][0]*a.m[1][1]-a.m[0][1]*a.m[1][0]) * idet;
    return r;
}

static void mat33_vec(mat33f m, float x, float y, float z,
                      float *ox, float *oy, float *oz)
{
    *ox = m.m[0][0]*x + m.m[0][1]*y + m.m[0][2]*z;
    *oy = m.m[1][0]*x + m.m[1][1]*y + m.m[1][2]*z;
    *oz = m.m[2][0]*x + m.m[2][1]*y + m.m[2][2]*z;
}

/* Create BLOK set - spatial neighborhood lists for local Pearson correlation.
   Ported from AFNI's mri_genalign_util.c create_GA_BLOK_set(). */
static GA_BLOK_set *create_GA_BLOK_set(
    int nx, int ny, int nz,
    float dx, float dy, float dz,
    int npt, float *im, float *jm, float *km,
    int bloktype, float blokrad, int minel,
    float shfac, int verb)
{
    GA_BLOK_set *gbs;
    float dxp, dyp, dzp, dxq, dyq, dzq, dxr, dyr, dzr, siz;
    float xx, yy, zz, uu, vv, ww;
    mat33f latmat, invlatmat;
    float px, py, pz;
    int pb, pt, qb, qt, rb, rt, pp, qq, rr, nblok, ii, nxy;
    int aa, bb, cc, dd, ss, np, nq, nr, npq;
    int *nelm, *nalm, **elm, ntot, nsav, ndup;
    int alloc_failed = 0;

    if (nx < 3 || ny < 3 || nz < 1) return NULL;
    if (dx <= 0.0f) dx = 1.0f;
    if (dy <= 0.0f) dy = 1.0f;
    if (dz <= 0.0f) dz = 1.0f;
    if (shfac < 0.2f || shfac > 5.0f) shfac = 1.0f;

    if (npt <= 0 || im == NULL || jm == NULL || km == NULL) {
        im = jm = km = NULL; npt = 0;
    }

    /* Set lattice vectors based on blok type */
    switch (bloktype) {
        case GA_BLOK_RHDD: {
            float a = blokrad;
            siz = a; a *= shfac;
            dxp = a;    dyp = a;    dzp = 0.0f;
            dxq = 0.0f; dyq = a;    dzq = a;
            dxr = a;    dyr = 0.0f; dzr = a;
        } break;
        case GA_BLOK_TOHD: {
            float a = blokrad;
            siz = a; a *= shfac;
            dxp = -a;  dyp =  a;  dzp =  a;
            dxq =  a;  dyq = -a;  dzq =  a;
            dxr =  a;  dyr =  a;  dzr = -a;
        } break;
        default: return NULL;
    }

    /* Build lattice matrix and its inverse */
    latmat.m[0][0] = dxp; latmat.m[0][1] = dxq; latmat.m[0][2] = dxr;
    latmat.m[1][0] = dyp; latmat.m[1][1] = dyq; latmat.m[1][2] = dyr;
    latmat.m[2][0] = dzp; latmat.m[2][1] = dzq; latmat.m[2][2] = dzr;
    invlatmat = mat33_inv(latmat);

    /* Find range of lattice indices by testing corners */
    float xt = (nx - 1) * dx, yt = (ny - 1) * dy, zt = (nz - 1) * dz;
    pb = pt = qb = qt = rb = rt = 0;

    /* Test 7 nonzero corners */
    float corners[7][3] = {
        {xt, 0, 0}, {xt, yt, 0}, {xt, 0, zt}, {xt, yt, zt},
        {0, yt, 0}, {0, 0, zt}, {0, yt, zt}
    };
    for (ii = 0; ii < 7; ii++) {
        mat33_vec(invlatmat, corners[ii][0], corners[ii][1], corners[ii][2],
                  &px, &py, &pz);
        pp = (int)floorf(px); pb = (pp < pb) ? pp : pb; pp++; pt = (pp > pt) ? pp : pt;
        qq = (int)floorf(py); qb = (qq < qb) ? qq : qb; qq++; qt = (qq > qt) ? qq : qt;
        rr = (int)floorf(pz); rb = (rr < rb) ? rr : rb; rr++; rt = (rr > rt) ? rr : rt;
    }

    np = pt - pb + 1;
    nq = qt - qb + 1; npq = np * nq;
    nr = rt - rb + 1;
    nblok = npq * nr;

    /* Create empty lists */
    nelm = (int *)  calloc(nblok, sizeof(int));
    nalm = (int *)  calloc(nblok, sizeof(int));
    elm  = (int **) calloc(nblok, sizeof(int *));
    if (!nelm || !nalm || !elm) { free(nelm); free(nalm); free(elm); return NULL; }

    nxy = nx * ny;
    if (npt == 0) npt = nxy * nz;

    /* Assign points to bloks */
    for (ndup = ntot = ii = 0; ii < npt; ii++) {
        if (im != NULL) {
            pp = (int)im[ii]; qq = (int)jm[ii]; rr = (int)km[ii];
        } else {
            pp = ii % nx; rr = ii / nxy; qq = (ii - rr * nxy) / nx;
        }
        ss = ii;
        xx = pp * dx; yy = qq * dy; zz = rr * dz;
        mat33_vec(invlatmat, xx, yy, zz, &px, &py, &pz);
        pp = (int)floorf(px + 0.499f);
        qq = (int)floorf(py + 0.499f);
        rr = (int)floorf(pz + 0.499f);

        int nsaved = 0;
        for (cc = rr - 1; cc <= rr + 1; cc++) {
            if (cc < rb || cc > rt) continue;
            for (bb = qq - 1; bb <= qq + 1; bb++) {
                if (bb < qb || bb > qt) continue;
                for (aa = pp - 1; aa <= pp + 1; aa++) {
                    if (aa < pb || aa > pt) continue;
                    float cx, cy, cz;
                    mat33_vec(latmat, (float)aa, (float)bb, (float)cc, &cx, &cy, &cz);
                    uu = xx - cx; vv = yy - cy; ww = zz - cz;
                    if (blok_inside(bloktype, uu, vv, ww, siz)) {
                        dd = (aa - pb) + (bb - qb) * np + (cc - rb) * npq;
                        /* Add to blok: expand array if needed */
                        if (nelm[dd] == nalm[dd]) {
                            int new_nalm = (int)(1.5f * nalm[dd]) + 16;
                            int *new_elm = (int *)realloc(elm[dd], sizeof(int) * new_nalm);
                            if (new_elm == NULL) {
                                alloc_failed = 1;
                                goto blok_alloc_fail;
                            }
                            elm[dd] = new_elm;
                            nalm[dd] = new_nalm;
                        }
                        elm[dd][nelm[dd]++] = ss;
                        ntot++; nsaved++;
                    }
                }
            }
        }
        if (nsaved > 1) ndup++;
    }

blok_alloc_fail:
    if (alloc_failed) {
        fprintf(stderr, "allineate: BLOK set allocation failed\n");
        for (dd = 0; dd < nblok; dd++) free(elm[dd]);
        free(nalm); free(nelm); free(elm);
        return NULL;
    }

    /* Compute minel if not specified */
    if (minel < 9) {
        for (minel = dd = 0; dd < nblok; dd++)
            minel = (nelm[dd] > minel) ? nelm[dd] : minel;
        minel = (int)(0.456 * minel) + 1;
    }

    /* Cull underpopulated bloks */
    for (nsav = dd = 0; dd < nblok; dd++) {
        if (nelm[dd] < minel) {
            if (elm[dd] != NULL) { free(elm[dd]); elm[dd] = NULL; }
            nelm[dd] = 0;
        } else {
            /* Clip array to actual size */
            if (nelm[dd] < nalm[dd] && nelm[dd] > 0) {
                int *new_elm = (int *)realloc(elm[dd], sizeof(int) * nelm[dd]);
                if (new_elm != NULL) elm[dd] = new_elm;
            }
            nsav++;
        }
    }
    free(nalm);

    if (nsav == 0) {
        fprintf(stderr, "allineate: BLOK set has 0 surviving bloks\n");
        for (dd = 0; dd < nblok; dd++) free(elm[dd]);
        free(nelm); free(elm); return NULL;
    }

    /* Build output struct */
    gbs = (GA_BLOK_set *)malloc(sizeof(GA_BLOK_set));
    if (gbs == NULL) {
        for (dd = 0; dd < nblok; dd++) free(elm[dd]);
        free(nelm); free(elm); return NULL;
    }
    gbs->num  = nsav;
    gbs->nelm = (int *)  calloc(nsav, sizeof(int));
    gbs->elm  = (int **) calloc(nsav, sizeof(int *));
    if (gbs->nelm == NULL || gbs->elm == NULL) {
        for (dd = 0; dd < nblok; dd++) free(elm[dd]);
        free(gbs->nelm);
        free(gbs->elm);
        free(gbs);
        free(nelm);
        free(elm);
        return NULL;
    }
    for (ntot = nsav = dd = 0; dd < nblok; dd++) {
        if (nelm[dd] > 0 && elm[dd] != NULL) {
            gbs->nelm[nsav] = nelm[dd]; ntot += nelm[dd];
            gbs->elm[nsav]  = elm[dd];  nsav++;
        }
    }
    free(nelm); free(elm);
    gbs->ppow = 1.0f;
    gbs->nx = nx; gbs->ny = ny; gbs->nz = nz;
    gbs->dx = dx; gbs->dy = dy; gbs->dz = dz;

    if (verb)
        fprintf(stderr, " + %d total points in %d bloks (rad=%.1f, %d dups)\n",
                ntot, gbs->num, blokrad, ndup);
    return gbs;
}

static void free_GA_BLOK_set(GA_BLOK_set *gbs)
{
    int i;
    if (gbs == NULL) return;
    if (gbs->elm) {
        for (i = 0; i < gbs->num; i++)
            if (gbs->elm[i]) free(gbs->elm[i]);
        free(gbs->elm);
    }
    if (gbs->nelm) free(gbs->nelm);
    free(gbs);
}

/*==========================================================================*/
/*========== SECTION 6b: HISTOGRAM CLIP LEVELS (CLEQWD from AFNI) =========*/
/*==========================================================================*/

#undef  WAY_BIG
#define WAY_BIG 1.e+10f

/* Compute background clip level from positive values in an array.
   Equivalent to AFNI's THD_cliplevel(im, mfrac).
   Algorithm: histogram positive values, find iterative median-fraction threshold. */
static float al_cliplevel(int n, const float *ar, float mfrac)
{
    if (n < 222 || ar == NULL) return 0.0f;
    if (mfrac <= 0.0f || mfrac >= 0.99f) mfrac = 0.50f;

    /* Find max of positive values and count them */
    float amax = 0.0f;
    int npos = 0;
    for (int ii = 0; ii < n; ii++) {
        if (ar[ii] > 0.0f) {
            npos++;
            if (ar[ii] > amax) amax = ar[ii];
        }
    }
    if (npos <= 222 || amax <= 0.0f) return 0.0f;

    /* Build histogram of positive values */
    int nhist = 10000;
    float sfac = (float)nhist / amax;
    int *hist = (int *)calloc(nhist + 1, sizeof(int));
    if (!hist) return 0.0f;
    double dsum = 0.0;
    for (int ii = 0; ii < n; ii++) {
        if (ar[ii] > 0.0f) {
            int kk = (int)(sfac * ar[ii] + 0.499f);
            if (kk > nhist) kk = nhist;
            hist[kk]++;
            dsum += (double)kk * (double)kk;
        }
    }

    /* Initial cut: RMS bin index halved */
    int ib = (int)(0.5 * sqrt(dsum / npos) + 0.5);

    /* Walk down from top, capturing ~65% of positive voxels */
    int qq = (int)(0.65f * npos);
    int kk = 0, ii;
    for (ii = nhist; ii >= ib && kk < qq; ii--)
        kk += hist[ii];
    int ncut = ii;

    /* Iterate: find median above cut, set new cut = mfrac * median */
    for (int iter = 0; iter < 66; iter++) {
        int nold = ncut;
        /* Count voxels at or above cut */
        npos = 0;
        for (ii = ncut; ii <= nhist; ii++) npos += hist[ii];
        int nhalf = npos / 2;
        /* Find median bin */
        kk = 0;
        for (ii = ncut; ii <= nhist && kk < nhalf; ii++)
            kk += hist[ii];
        ncut = (int)(mfrac * ii);
        if (ncut == nold) break;
    }

    free(hist);
    return (float)ncut / sfac;
}

/* AFNI zero-pads the fixed image so supra-threshold content has at least
   max(8,n/8) background voxels on every face. Return those six pad widths so
   the compact engine can derive AFNI's expanded translation range. */
static void al_fixed_pad_widths(const float *im, int nx, int ny, int nz,
                                int *pxm, int *pxp, int *pym, int *pyp,
                                int *pzm, int *pzp)
{
    *pxm = *pxp = *pym = *pyp = *pzm = *pzp = 0;
    if (im == NULL || nx < 1 || ny < 1 || nz < 1) return;

    size_t nvox = (size_t)nx * ny * nz;
    if (nvox > INT_MAX) return;
    float cv = 0.33f * al_cliplevel((int)nvox, im, 0.33f);
    int xlo = nx, ylo = ny, zlo = nz, xhi = -1, yhi = -1, zhi = -1;
    for (int z = 0; z < nz; z++) for (int y = 0; y < ny; y++) for (int x = 0; x < nx; x++) {
        float v = im[x + (size_t)y * nx + (size_t)z * nx * ny];
        /* MRI_autobbox sees the thresholded image, so a true zero is not content
           even in the degenerate cv==0 case. */
        if (!al_finitef(v) || v == 0.0f || v < cv) continue;
        if (x < xlo) xlo = x; if (x > xhi) xhi = x;
        if (y < ylo) ylo = y; if (y > yhi) yhi = y;
        if (z < zlo) zlo = z; if (z > zhi) zhi = z;
    }
    if (xhi < xlo || yhi < ylo || zhi < zlo) return;

    int mpad = 8;
    if (nx / 8 > mpad) mpad = nx / 8;
    if (ny / 8 > mpad) mpad = ny / 8;
    if (nz / 8 > mpad) mpad = nz / 8;
    *pxm = mpad - xlo; if (*pxm < 0) *pxm = 0;
    *pym = mpad - ylo; if (*pym < 0) *pym = 0;
    *pzm = mpad - zlo; if (*pzm < 0) *pzm = 0;
    *pxp = mpad - (nx - 1 - xhi); if (*pxp < 0) *pxp = 0;
    *pyp = mpad - (ny - 1 - yhi); if (*pyp < 0) *pyp = 0;
    *pzp = mpad - (nz - 1 - zhi); if (*pzp < 0) *pzp = 0;
    if (nz == 1) *pzm = *pzp = 0;  /* AFNI does not pad through-plane in 2D. */
}

static void al_shift_range_mm(mat44 cmat, int nx, int ny, int nz,
                              float *xmm, float *ymm, float *zmm)
{
    float x = 0.321f * (nx - 1), y = 0.321f * (ny - 1), z = 0.321f * (nz - 1);
    *xmm = *ymm = *zmm = 0.01f;
    for (int ii = -1; ii <= 1; ii += 2)
        for (int jj = -1; jj <= 1; jj += 2)
            for (int kk = -1; kk <= 1; kk += 2) {
                float xp = fabsf(cmat.m[0][0]*(ii*x) + cmat.m[0][1]*(jj*y) + cmat.m[0][2]*(kk*z));
                float yp = fabsf(cmat.m[1][0]*(ii*x) + cmat.m[1][1]*(jj*y) + cmat.m[1][2]*(kk*z));
                float zp = fabsf(cmat.m[2][0]*(ii*x) + cmat.m[2][1]*(jj*y) + cmat.m[2][2]*(kk*z));
                if (xp > *xmm) *xmm = xp;
                if (yp > *ymm) *ymm = yp;
                if (zp > *zmm) *zmm = zp;
            }
}

/* Compute the p-th quantile (0 <= p <= 1) of an array.
   O(n) selection (al_select_rank); allocates a temporary copy. */
static float al_quantile(int n, const float *ar, float p)
{
    if (n <= 0 || ar == NULL) return 0.0f;
    if (p <= 0.0f) p = 0.0f;
    if (p >= 1.0f) p = 1.0f;

    float *tmp = (float *)malloc(n * sizeof(float));
    if (!tmp) return 0.0f;
    int nn = 0;
    for (int ii = 0; ii < n; ii++)
        if (ar[ii] < WAY_BIG) tmp[nn++] = ar[ii];
    if (nn == 0) { free(tmp); return 0.0f; }

    int idx = (int)(p * (nn - 1));
    float val = al_select_rank(tmp, nn, idx);   /* O(n) select, no qsort comparator */
    free(tmp);
    return val;
}

/* Compute clip levels for histogram-based cost functions (CLEQWD mode).
   Equivalent to AFNI's clipate() from thd_correlate.c.
   For positive images: bot = THD_cliplevel(0.345), top = quantile(0.966),
   capped at top = min(top, 4.321 * bot).
   Returns bot in pair.a, top in pair.b. If a >= b, clipping is invalid. */
typedef struct { float a, b; } al_float_pair;

static al_float_pair al_clipate(int n, const float *ar)
{
    al_float_pair rr = { 1.0f, 0.0f }; /* invalid by default */
    if (n < 666 || ar == NULL) return rr;

    /* Compute clip levels from positive values (al_cliplevel handles this).
       Do NOT bail out for images with some negative values — negative values
       are simply ignored by the clipping computation, as in AFNI's clipate. */

    float cbot = al_cliplevel(n, ar, 0.345f);
    float ctop = al_quantile(n, ar, 0.966f);
    if (ctop > 4.321f * cbot) ctop = 4.321f * cbot;

    if (cbot >= ctop) return rr; /* invalid */
    rr.a = cbot;
    rr.b = ctop;
    return rr;
}

/*==========================================================================*/
/*========== SECTION 7: 2D HISTOGRAM AND CORRELATION METRICS ==============*/
/*==========================================================================*/

/* Shannon entropy function */
#define SHANENT(z) (((z) <= 0.0f) ? 0.0f : (z)*logf(z))

/* Static 2D histogram state (thread-local for OpenMP safety) */
static AL_TLOCAL float *al_xc = NULL, *al_yc = NULL, *al_xyc = NULL;
static AL_TLOCAL float al_nww = 0.0f;
static AL_TLOCAL int al_nbin = 0, al_nbp = 0, al_nbm = 0;
static AL_TLOCAL int al_nbp_cap = 0;
static AL_TLOCAL unsigned char *al_good = NULL;  /* build_2Dhist per-eval good[] mask, reused */
static AL_TLOCAL int al_good_cap = 0;
static AL_TLOCAL int al_hist_oom = 0;   /* set ONLY on a histogram-buffer alloc failure in
                                           build_2Dhist, so the cost path can distinguish a real
                                           OOM (reject → AL_BIGVAL) from a legitimately empty
                                           histogram (n<=9 / constant range / no overlap → cost 0). */
static double al_hpow = 0.33333333333;

#undef  XYC
#define XYC(p,q) al_xyc[(p)+(q)*al_nbp]

#undef  WAY_BIG
#define WAY_BIG 1.e+10f
#undef  GOODVAL
#define GOODVAL(x) ((x) < WAY_BIG)
#undef  RANGVAL
#define RANGVAL(x,b,t) ((x) >= (b) && (x) <= (t))
#undef  WW
#define WW(i) ((w==NULL) ? 1.0f : w[i])

static void reset_2Dhist_state(void)
{
    al_nbin = al_nbp = al_nbm = 0; al_nww = 0.0f;
}

static void clear_2Dhist(void)
{
    if (al_xc)  { free(al_xc);  al_xc = NULL; }
    if (al_yc)  { free(al_yc);  al_yc = NULL; }
    if (al_xyc) { free(al_xyc); al_xyc = NULL; }
    if (al_good) { free(al_good); al_good = NULL; }
    al_nbp_cap = 0; al_good_cap = 0;
    reset_2Dhist_state();
}

static int ensure_2Dhist_capacity(int nbp)
{
    if (nbp <= 0) return 1;
    if (al_nbp_cap >= nbp && al_xc != NULL && al_yc != NULL && al_xyc != NULL)
        return 0;

    size_t nb = (size_t)nbp;
    float *new_xc = (float *)malloc(nb * sizeof(float));
    float *new_yc = (float *)malloc(nb * sizeof(float));
    float *new_xyc = (float *)malloc(nb * nb * sizeof(float));
    if (!new_xc || !new_yc || !new_xyc) {
        free(new_xc); free(new_yc); free(new_xyc);
        clear_2Dhist();
        return 1;
    }

    free(al_xc); free(al_yc); free(al_xyc);
    al_xc = new_xc; al_yc = new_yc; al_xyc = new_xyc;
    al_nbp_cap = nbp;
    return 0;
}

/* Build 2D histogram with CLEQWD edge-bin mode (matching AFNI's xyclip path).
   xbot/xtop and ybot/ytop are CLEQWD clip ranges. Values outside the clip
   range are binned into edge bins (0 and nbm) rather than excluded, keeping
   all data points in the histogram while concentrating resolution on the
   informative intensity range. */
static void build_2Dhist(int n, float xbot, float xtop, float *x,
                                float ybot, float ytop, float *y, float *w,
                                float xmemtop, float ymemtop)
{
    int ii, jj, kk, ngood;
    float xi, yi, xx, yy, x1, y1, ww;
    unsigned char *good;
    float xdbot, xdtop, ydbot, ydtop; /* actual data range */

    al_hist_oom = 0;   /* cleared at entry; set below ONLY on a genuine malloc failure */
    if (n <= 9 || x == NULL || y == NULL) return;

    reset_2Dhist_state();

    /* Reuse the thread-local good[] mask across evals (build_2Dhist runs once per Hellinger cost
       evaluation — a per-eval malloc/free in the hot path otherwise). Grow-only; freed by
       clear_2Dhist() alongside the histogram buffers. */
    if (al_good_cap < n) {
        free(al_good);
        al_good = (unsigned char *)malloc((size_t)n);
        al_good_cap = al_good ? n : 0;
    }
    good = al_good;
    if (!good) { al_hist_oom = 1; return; }   /* OOM, not a legitimately empty histogram */
    /* Histogram MEMBERSHIP: drop any pair with a value above its image's mri_topclip
       (xmemtop/ymemtop). This is AFNI's `good[]` over [min, topclip] — it excludes bright
       outliers (e.g. a T2 source's CSF/fat/orbits, non-brain the T1 base lacks) so they can't
       bias the fit toward a shrunk overlap. xmemtop<=0 disables the drop for that axis. */
    for (ii = 0; ii < n; ii++)
        good[ii] = GOODVAL(x[ii]) && GOODVAL(y[ii])
                   && (xmemtop <= 0.0f || x[ii] <= xmemtop)
                   && (ymemtop <= 0.0f || y[ii] <= ymemtop);

    /* Find actual data range */
    xdbot = WAY_BIG; xdtop = -WAY_BIG;
    ydbot = WAY_BIG; ydtop = -WAY_BIG;
    for (ii = 0; ii < n; ii++) {
        if (!good[ii]) continue;
        if (x[ii] > xdtop) xdtop = x[ii];
        if (x[ii] < xdbot) xdbot = x[ii];
        if (y[ii] > ydtop) ydtop = y[ii];
        if (y[ii] < ydbot) ydbot = y[ii];
    }
    if (xdbot >= xdtop || ydbot >= ydtop) { return; }

    /* Count good values in actual data range */
    memset(good, 0, n);
    for (ngood = ii = 0; ii < n; ii++) {
        if (RANGVAL(x[ii], xdbot, xdtop) && RANGVAL(y[ii], ydbot, ydtop) && WW(ii) > 0.0f) {
            good[ii] = 1; ngood++;
        }
    }
    if (ngood == 0) { return; }

    /* Compute number of bins from total n (matches AFNI) */
    al_nbin = (int)pow((double)n, al_hpow);
    if (al_nbin > 255) al_nbin = 255;
    else if (al_nbin < 3) al_nbin = 3;
    al_nbp = al_nbin + 1;
    al_nbm = al_nbin - 1;

    if (ensure_2Dhist_capacity(al_nbp)) { al_hist_oom = 1; return; }   /* histogram-buffer OOM */
    memset(al_xc, 0, (size_t)al_nbp * sizeof(float));
    memset(al_yc, 0, (size_t)al_nbp * sizeof(float));
    memset(al_xyc, 0, (size_t)al_nbp * (size_t)al_nbp * sizeof(float));
    /* al_nww already zeroed by reset_2Dhist_state() at function entry */

    /* Determine if CLEQWD clipping is active and data extends beyond clips */
    int xyclip = (xbot < xtop) &&
                 (xdbot < xbot) && (xbot < xtop) && (xtop < xdtop) &&
                 (ydbot < ybot) && (ybot < ytop) && (ytop < ydtop);

    if (xyclip) {
        /* CLEQWD edge-bin mode (AFNI's xyclip path):
           values below clip go to bin 0, above clip go to bin nbm,
           clip range is divided into bins 1..nbm-1 with linear spread */
        xi = (al_nbin - 2.000001f) / (xtop - xbot);
        yi = (al_nbin - 2.000001f) / (ytop - ybot);
        for (ii = 0; ii < n; ii++) {
            if (!good[ii]) continue;
            /* X binning with edge bins */
            xx = x[ii];
            if (xx < xbot) { jj = 0; xx = 0.0f; }
            else if (xx > xtop) { jj = al_nbm; xx = 1.0f; }
            else { xx = 1.0f + (xx - xbot) * xi; jj = (int)xx; xx = xx - jj; }
            /* Y binning with edge bins */
            yy = y[ii];
            if (yy < ybot) { kk = 0; yy = 0.0f; }
            else if (yy > ytop) { kk = al_nbm; yy = 1.0f; }
            else { yy = 1.0f + (yy - ybot) * yi; kk = (int)yy; yy = yy - kk; }

            x1 = 1.0f - xx; y1 = 1.0f - yy;
            ww = WW(ii); al_nww += ww;

            al_xc[jj] += x1 * ww; al_xc[jj + 1] += xx * ww;
            al_yc[kk] += y1 * ww; al_yc[kk + 1] += yy * ww;

            XYC(jj,     kk    ) += x1 * (y1 * ww);
            XYC(jj + 1, kk    ) += xx * (y1 * ww);
            XYC(jj,     kk + 1) += x1 * (yy * ww);
            XYC(jj + 1, kk + 1) += xx * (yy * ww);
        }
    } else {
        /* Equal-size bins with linear spread (no CLEQWD or data within clip) */
        float xb = xdbot, yb = ydbot;
        xi = al_nbm / (xdtop - xdbot);
        yi = al_nbm / (ydtop - ydbot);
        for (ii = 0; ii < n; ii++) {
            if (!good[ii]) continue;
            xx = (x[ii] - xb) * xi;
            jj = (int)xx; xx = xx - jj; x1 = 1.0f - xx;
            yy = (y[ii] - yb) * yi;
            kk = (int)yy; yy = yy - kk; y1 = 1.0f - yy;
            ww = WW(ii); al_nww += ww;

            al_xc[jj] += x1 * ww; al_xc[jj + 1] += xx * ww;
            al_yc[kk] += y1 * ww; al_yc[kk + 1] += yy * ww;

            XYC(jj,     kk    ) += x1 * (y1 * ww);
            XYC(jj + 1, kk    ) += xx * (y1 * ww);
            XYC(jj,     kk + 1) += x1 * (yy * ww);
            XYC(jj + 1, kk + 1) += xx * (yy * ww);
        }
    }
}

static void normalize_2Dhist(void)
{
    if (al_nww > 0.0f && al_xyc != NULL && al_xc != NULL && al_yc != NULL) {
        float ni = 1.0f / al_nww;
        int nbq = al_nbp * al_nbp, ii;
        for (ii = 0; ii < al_nbp; ii++) { al_xc[ii] *= ni; al_yc[ii] *= ni; }
        for (ii = 0; ii < nbq; ii++) { al_xyc[ii] *= ni; }
    }
}

/* Combined Hellinger, MI, NMI, CRA from built histogram
   Returns: .a=hel, .b=mi, .c=nmi, .d=crA */
typedef struct { float a, b, c, d; } float_quad;

static float_quad al_helmicra(int n, float xbot, float xtop, float *x,
                              float ybot, float ytop, float *y, float *w,
                              float xmemtop, float ymemtop)
{
    int ii, jj;
    float hel, pq, vv, uu;
    float cyvar, uyvar, yrat, xrat;
    float_quad hmc = {0.0f, 0.0f, 0.0f, 0.0f};

    build_2Dhist(n, xbot, xtop, x, ybot, ytop, y, w, xmemtop, ymemtop);
    if (al_nbin <= 0 || al_nww <= 0) return hmc;
    normalize_2Dhist();

    /* Hellinger, MI, NMI */
    hel = vv = uu = 0.0f;
    for (ii = 0; ii < al_nbp; ii++) {
        vv += SHANENT(al_xc[ii]) + SHANENT(al_yc[ii]);
        for (jj = 0; jj < al_nbp; jj++) {
            pq = XYC(ii, jj);
            hel += sqrtf(pq * al_xc[ii] * al_yc[jj]);
            uu  += SHANENT(pq);
        }
    }
    hmc.a = 1.0f - hel;                     /* Hellinger */
    hmc.b = uu - vv;                         /* MI */
    hmc.c = (vv != 0.0f) ? uu / vv : 0.0f;  /* NMI */

    /* CR(y|x) */
    cyvar = 0.0f;
    for (ii = 0; ii < al_nbp; ii++) {
        if (al_xc[ii] > 0.0f) {
            float mm = 0.0f, vvv = 0.0f;
            for (jj = 1; jj < al_nbp; jj++) {
                mm  += jj * XYC(ii, jj);
                vvv += jj * (jj * XYC(ii, jj));
            }
            cyvar += vvv - mm * mm / al_xc[ii];
        }
    }
    vv = 0.0f; uu = 0.0f;
    for (jj = 1; jj < al_nbp; jj++) {
        uu += jj * al_yc[jj]; vv += jj * (jj * al_yc[jj]);
    }
    uyvar = vv - uu * uu;
    yrat = (uyvar > 0.0f) ? cyvar / uyvar : 1.0f;

    /* CR(x|y) */
    cyvar = 0.0f;
    for (jj = 0; jj < al_nbp; jj++) {
        if (al_yc[jj] > 0.0f) {
            float mm = 0.0f, vvv = 0.0f;
            for (ii = 1; ii < al_nbp; ii++) {
                mm  += ii * XYC(ii, jj);
                vvv += ii * (ii * XYC(ii, jj));
            }
            cyvar += vvv - mm * mm / al_yc[jj];
        }
    }
    vv = 0.0f; uu = 0.0f;
    for (ii = 1; ii < al_nbp; ii++) {
        uu += ii * al_xc[ii]; vv += ii * (ii * al_xc[ii]);
    }
    uyvar = vv - uu * uu;
    xrat = (uyvar > 0.0f) ? cyvar / uyvar : 1.0f;

    hmc.d = 1.0f - 0.5f * (xrat + yrat); /* additive CRA */
    return hmc;
}

/*==========================================================================*/
/*================== SECTION 8: CORE ALIGNMENT ENGINE =====================*/
/*==========================================================================*/

/*--- Rotation matrix: Q = R3*R2*R1 where Ri rotates about axis axi by thi ---*/
static mat44 rot_matrix(int ax1, double th1,
                        int ax2, double th2, int ax3, double th3)
{
    mat44 q, p, r;
    double c, s;

    /* R1 */
    memset(&q, 0, sizeof(mat44)); q.m[3][3] = 1.0f;
    c = cos(th1); s = sin(th1);
    switch (ax1) {
        case 0: q.m[0][0]=1; q.m[1][1]=c; q.m[1][2]=-s; q.m[2][1]=s; q.m[2][2]=c; break;
        case 1: q.m[1][1]=1; q.m[0][0]=c; q.m[0][2]=s; q.m[2][0]=-s; q.m[2][2]=c; break;
        case 2: q.m[2][2]=1; q.m[0][0]=c; q.m[0][1]=-s; q.m[1][0]=s; q.m[1][1]=c; break;
    }

    /* R2 */
    memset(&p, 0, sizeof(mat44)); p.m[3][3] = 1.0f;
    c = cos(th2); s = sin(th2);
    switch (ax2) {
        case 0: p.m[0][0]=1; p.m[1][1]=c; p.m[1][2]=-s; p.m[2][1]=s; p.m[2][2]=c; break;
        case 1: p.m[1][1]=1; p.m[0][0]=c; p.m[0][2]=s; p.m[2][0]=-s; p.m[2][2]=c; break;
        case 2: p.m[2][2]=1; p.m[0][0]=c; p.m[0][1]=-s; p.m[1][0]=s; p.m[1][1]=c; break;
    }
    q = nifti_mat44_mul(p, q);

    /* R3 */
    memset(&r, 0, sizeof(mat44)); r.m[3][3] = 1.0f;
    c = cos(th3); s = sin(th3);
    switch (ax3) {
        case 0: r.m[0][0]=1; r.m[1][1]=c; r.m[1][2]=-s; r.m[2][1]=s; r.m[2][2]=c; break;
        case 1: r.m[1][1]=1; r.m[0][0]=c; r.m[0][2]=s; r.m[2][0]=-s; r.m[2][2]=c; break;
        case 2: r.m[2][2]=1; r.m[0][0]=c; r.m[0][1]=-s; r.m[1][0]=s; r.m[1][1]=c; break;
    }
    q = nifti_mat44_mul(r, q);
    return q;
}

/*--- Build affine mat44 from 12 parameters:
      [0..2]=shifts, [3..5]=rotations(deg), [6..8]=scales, [9..11]=shears ---*/
static mat44 GA_setup_affine(int npar, float *parvec)
{
    mat44 ss, dd, uu, aa, bb, gam;
    float a, b, c, p, q, r;

    /* Rotation */
    a = b = c = 0.0f;
    if (npar >= 4) a = D2R * parvec[3];
    if (npar >= 5) b = D2R * parvec[4];
    if (npar >= 6) c = D2R * parvec[5];
    if (a != 0.0f || b != 0.0f || c != 0.0f)
        uu = rot_matrix(2, a, 0, b, 1, c);
    else
        uu = mat44_diag(1.0f, 1.0f, 1.0f);

    /* Scaling */
    a = b = c = 1.0f;
    if (npar >= 7) { a = parvec[6]; if (a <= 0.10f || a >= 10.0f) a = 1.0f; }
    if (al_zoom_isotropic) {
        b = c = a;   /* -zoom: y,z scale follow x-scale -> one global isotropic factor */
    } else {
        if (npar >= 8) { b = parvec[7]; if (b <= 0.10f || b >= 10.0f) b = 1.0f; }
        if (npar >= 9) { c = parvec[8]; if (c <= 0.10f || c >= 10.0f) c = 1.0f; }
    }
    dd = mat44_diag(a, b, c);

    /* Shear (lower triangular by default) */
    a = b = c = 0.0f;
    if (npar >= 10) { a = parvec[9];  if (fabsf(a) > 0.3333f) a = 0.0f; }
    if (npar >= 11) { b = parvec[10]; if (fabsf(b) > 0.3333f) b = 0.0f; }
    if (npar >= 12) { c = parvec[11]; if (fabsf(c) > 0.3333f) c = 0.0f; }
    memset(&ss, 0, sizeof(mat44)); ss.m[3][3] = 1.0f;
    switch (AL_SMAT) {
        default:
        case SMAT_LOWER:
            ss.m[0][0] = 1.0f;
            ss.m[1][0] = a; ss.m[1][1] = 1.0f;
            ss.m[2][0] = b; ss.m[2][1] = c; ss.m[2][2] = 1.0f;
            break;
        case SMAT_UPPER:
            ss.m[0][0] = 1.0f; ss.m[0][1] = a; ss.m[0][2] = b;
            ss.m[1][1] = 1.0f; ss.m[1][2] = c;
            ss.m[2][2] = 1.0f;
            break;
    }

    /* Multiply as ordered: default SDU */
    switch (AL_MATORDER) {
        default:
        case MATORDER_SDU: aa = nifti_mat44_mul(ss, dd); bb = uu; break;
        case MATORDER_SUD: aa = nifti_mat44_mul(ss, uu); bb = dd; break;
        case MATORDER_DSU: aa = nifti_mat44_mul(dd, ss); bb = uu; break;
        case MATORDER_DUS: aa = nifti_mat44_mul(dd, uu); bb = ss; break;
        case MATORDER_USD: aa = nifti_mat44_mul(uu, ss); bb = dd; break;
        case MATORDER_UDS: aa = nifti_mat44_mul(uu, dd); bb = ss; break;
    }
    gam = nifti_mat44_mul(aa, bb);

    /* Shifts */
    a = b = c = 0.0f;
    if (npar >= 1) a = parvec[0];
    if (npar >= 2) b = parvec[1];
    if (npar >= 3) c = parvec[2];
    if (AL_DCODE == DELTA_BEFORE) {
        mat44_vec(gam, a, b, c, &p, &q, &r);
        a = p; b = q; c = r;
    }
    gam.m[0][3] = a; gam.m[1][3] = b; gam.m[2][3] = c;

    /* Before/after coordinate transforms */
    if (aff_use_before) gam = nifti_mat44_mul(gam, aff_before);
    if (aff_use_after)  gam = nifti_mat44_mul(aff_after, gam);

    return gam;
}

/*--- Warp function for affine transforms ---*/
static void al_wfunc_affine(int npar, float *wpar,
                            int npt, float *xi, float *yi, float *zi,
                                     float *xo, float *yo, float *zo)
{
    static AL_TLOCAL mat44 gam;
    int ii;

    if (npar > 0 && wpar != NULL)
        gam = GA_setup_affine(npar, wpar);

    if (npt <= 0 || xi == NULL || xo == NULL) return;

    for (ii = 0; ii < npt; ii++)
        MAT44_VEC(gam, xi[ii], yi[ii], zi[ii], xo[ii], yo[ii], zo[ii]);
}

/*--- Interpolate at a set of points ---*/
static void al_interp(float *far, int nx, int ny, int nz,
                      int interp_code, int npp,
                      float *ip, float *jp, float *kp, float *vv)
{
    switch (interp_code) {
        case AL_INTERP_NN:
            GA_interp_NN(far, nx, ny, nz, npp, ip, jp, kp, vv); break;
        case AL_INTERP_LINEAR:
            GA_interp_linear(far, nx, ny, nz, npp, ip, jp, kp, vv); break;
        case AL_INTERP_CUBIC:
        default:
            GA_interp_cubic(far, nx, ny, nz, npp, ip, jp, kp, vv); break;
    }
}

/*--- Fused incremental warp + interpolation for all-voxels case ---*/
/* When processing all base voxels in order, the affine warp is separable:
   xo = m[0][0]*i + (m[0][1]*j + m[0][2]*k + m[0][3])
   We compute the row base once per (j,k), then step by m[*][0] per voxel.
   This eliminates intermediate coordinate arrays and the chunking overhead.

   For linear interpolation, a per-row safe-range pre-pass identifies the
   interior voxels where all warped coordinates are guaranteed in-bounds.
   The interior loop skips bounds checking and CLIP, replacing floorf() with
   direct (int) casts. Edge voxels still use the full bounds-checked path. */

/* Compute the i-range [i_lo, i_hi) where all three warped coordinates stay
   inside [0.5, dim-1.5], guaranteeing (int)coord in [0, dim-2] and +1 in [1, dim-1].
   Returns 0 if no valid range exists. */
/* Compute range of output voxel indices ii where all three source coordinates
   fall within [coord_lo, coord_hi] (safe for unchecked interpolation).
   coord_lo = 0.5 for linear, 1.0 for cubic. coord_hi = *_safe parameter. */
static inline int al_safe_irange_ex(int bnx,
    float bj_x, float mx0, float anx_safe,
    float bj_y, float my0, float any_safe,
    float bj_z, float mz0, float anz_safe,
    float coord_lo,
    int *i_lo_out, int *i_hi_out)
{
    float lo = 0.0f, hi = (float)(bnx - 1);
    float bases[3] = {bj_x, bj_y, bj_z};
    float steps[3] = {mx0, my0, mz0};
    float highs[3] = {anx_safe, any_safe, anz_safe};

    for (int a = 0; a < 3; a++) {
        float b = bases[a], s = steps[a], h = highs[a];
        if (s > 1e-12f) {
            float t0 = (coord_lo - b) / s, t1 = (h - b) / s;
            if (t0 > lo) lo = t0;
            if (t1 < hi) hi = t1;
        } else if (s < -1e-12f) {
            float t0 = (h - b) / s, t1 = (coord_lo - b) / s;
            if (t0 > lo) lo = t0;
            if (t1 < hi) hi = t1;
        } else {
            if (b < coord_lo || b > h) return 0;
        }
    }
    int ilo = (int)ceilf(lo);
    int ihi = (int)floorf(hi) + 1;
    if (ilo < 0) ilo = 0;
    if (ihi > bnx) ihi = bnx;
    if (ilo >= ihi) return 0;
    *i_lo_out = ilo;
    *i_hi_out = ihi;
    return 1;
}

static inline int al_safe_irange(int bnx,
    float bj_x, float mx0, float anx_safe,
    float bj_y, float my0, float any_safe,
    float bj_z, float mz0, float anz_safe,
    int *i_lo_out, int *i_hi_out)
{
    return al_safe_irange_ex(bnx, bj_x, mx0, anx_safe, bj_y, my0, any_safe,
                              bj_z, mz0, anz_safe, 0.5f, i_lo_out, i_hi_out);
}

/* Linear interpolation for a single voxel with bounds check and CLIP */
static inline float al_interp_linear_checked(
    float xx, float yy, float zz,
    float nxh, float nyh, float nzh,
    int anx1, int any1, int anz1,
    const float * restrict aim, int anx, int nxy_a, int *oob)
{
    if (xx < -0.499f || xx > nxh || yy < -0.499f || yy > nyh ||
        zz < -0.499f || zz > nzh) { *oob = 1; return 0.0f; }
    *oob = 0;
    float ix = floorf(xx), fx = xx - ix;
    float jy = floorf(yy), fy = yy - jy;
    float kz = floorf(zz), fz = zz - kz;
    int ix_00 = (int)ix, ix_p1 = ix_00 + 1;
    int jy_00 = (int)jy, jy_p1 = jy_00 + 1;
    int kz_00 = (int)kz, kz_p1 = kz_00 + 1;
    CLIP(ix_00, anx1); CLIP(ix_p1, anx1);
    CLIP(jy_00, any1); CLIP(jy_p1, any1);
    CLIP(kz_00, anz1); CLIP(kz_p1, anz1);
    float w0 = 1.0f - fx, w1 = fx;
/* outer parens are load-bearing: XLINT is used as a multiplicand, e.g. fy*XLINT(...),
 * so without them the weight would bind only to the first term (the w1 term would escape). */
#define XLINT(j,k) (w0*aim[(ix_00)+(j)*anx+(k)*nxy_a]+w1*aim[(ix_p1)+(j)*anx+(k)*nxy_a])
    float fk0 = (1.0f-fy)*XLINT(jy_00,kz_00) + fy*XLINT(jy_p1,kz_00);
    float fk1 = (1.0f-fy)*XLINT(jy_00,kz_p1) + fy*XLINT(jy_p1,kz_p1);
#undef XLINT
    return (1.0f - fz) * fk0 + fz * fk1;
}

/* Cubic interpolation for a single voxel with bounds check and CLIP */
static inline float al_interp_cubic_checked(
    float xx, float yy, float zz,
    float nxh, float nyh, float nzh,
    int anx1, int any1, int anz1,
    const float * restrict aim, int anx, int nxy_a, int *oob)
{
    if (xx < -0.499f || xx > nxh || yy < -0.499f || yy > nyh ||
        zz < -0.499f || zz > nzh) { *oob = 1; return 0.0f; }
    *oob = 0;
    int cix = (int)floorf(xx); float cfx = xx - cix;
    int cjy = (int)floorf(yy); float cfy = yy - cjy;
    int ckz = (int)floorf(zz); float cfz = zz - ckz;
    if (ISTINY(cfx) && ISTINY(cfy) && ISTINY(cfz)) {
        CLIP(cix, anx1); CLIP(cjy, any1); CLIP(ckz, anz1);
        return aim[cix + cjy * anx + ckz * nxy_a];
    }
    int cix_m1 = cix-1, cix_00 = cix, cix_p1 = cix+1, cix_p2 = cix+2;
    CLIP(cix_m1, anx1); CLIP(cix_00, anx1); CLIP(cix_p1, anx1); CLIP(cix_p2, anx1);
    int cjy_m1 = cjy-1, cjy_00 = cjy, cjy_p1 = cjy+1, cjy_p2 = cjy+2;
    CLIP(cjy_m1, any1); CLIP(cjy_00, any1); CLIP(cjy_p1, any1); CLIP(cjy_p2, any1);
    int ckz_m1 = ckz-1, ckz_00 = ckz, ckz_p1 = ckz+1, ckz_p2 = ckz+2;
    CLIP(ckz_m1, anz1); CLIP(ckz_00, anz1); CLIP(ckz_p1, anz1); CLIP(ckz_p2, anz1);
    float cwx_m1 = P_M1(cfx), cwx_00 = P_00(cfx), cwx_p1 = P_P1(cfx), cwx_p2 = P_P2(cfx);
#define XINTC_(j,k) cwx_m1*aim[(cix_m1)+(j)*anx+(k)*nxy_a]+cwx_00*aim[(cix_00)+(j)*anx+(k)*nxy_a]\
                    +cwx_p1*aim[(cix_p1)+(j)*anx+(k)*nxy_a]+cwx_p2*aim[(cix_p2)+(j)*anx+(k)*nxy_a]
    float f_jm1_km1=XINTC_(cjy_m1,ckz_m1), f_j00_km1=XINTC_(cjy_00,ckz_m1);
    float f_jp1_km1=XINTC_(cjy_p1,ckz_m1), f_jp2_km1=XINTC_(cjy_p2,ckz_m1);
    float f_jm1_k00=XINTC_(cjy_m1,ckz_00), f_j00_k00=XINTC_(cjy_00,ckz_00);
    float f_jp1_k00=XINTC_(cjy_p1,ckz_00), f_jp2_k00=XINTC_(cjy_p2,ckz_00);
    float f_jm1_kp1=XINTC_(cjy_m1,ckz_p1), f_j00_kp1=XINTC_(cjy_00,ckz_p1);
    float f_jp1_kp1=XINTC_(cjy_p1,ckz_p1), f_jp2_kp1=XINTC_(cjy_p2,ckz_p1);
    float f_jm1_kp2=XINTC_(cjy_m1,ckz_p2), f_j00_kp2=XINTC_(cjy_00,ckz_p2);
    float f_jp1_kp2=XINTC_(cjy_p1,ckz_p2), f_jp2_kp2=XINTC_(cjy_p2,ckz_p2);
    float cwy_m1=P_M1(cfy), cwy_00=P_00(cfy), cwy_p1=P_P1(cfy), cwy_p2=P_P2(cfy);
    float f_km1 = cwy_m1*f_jm1_km1+cwy_00*f_j00_km1+cwy_p1*f_jp1_km1+cwy_p2*f_jp2_km1;
    float f_k00 = cwy_m1*f_jm1_k00+cwy_00*f_j00_k00+cwy_p1*f_jp1_k00+cwy_p2*f_jp2_k00;
    float f_kp1 = cwy_m1*f_jm1_kp1+cwy_00*f_j00_kp1+cwy_p1*f_jp1_kp1+cwy_p2*f_jp2_kp1;
    float f_kp2 = cwy_m1*f_jm1_kp2+cwy_00*f_j00_kp2+cwy_p1*f_jp1_kp2+cwy_p2*f_jp2_kp2;
    float cwz_m1=P_M1(cfz), cwz_00=P_00(cfz), cwz_p1=P_P1(cfz), cwz_p2=P_P2(cfz);
    return P_FACTOR * (cwz_m1*f_km1 + cwz_00*f_k00 + cwz_p1*f_kp1 + cwz_p2*f_kp2);
#undef XINTC_
}

static void GA_warp_interp_fused(mat44 gam, const float * restrict aim, int anx, int any, int anz,
                                  int bnx, int bny, int bnz, int interp_code,
                                  float * restrict avm)
{
    int nxy_b = bnx * bny;
    int nxy_a = anx * any;
    float nxh = anx - 0.501f, nyh = any - 0.501f, nzh = anz - 0.501f;
    int anx1 = anx - 1, any1 = any - 1, anz1 = anz - 1;
    /* Safe upper bounds for interior range: coord <= this means (int)coord+1 <= dim-1 */
    float anx_safe = (float)(anx - 2) + 0.999f;
    float any_safe = (float)(any - 2) + 0.999f;
    float anz_safe = (float)(anz - 2) + 0.999f;
    float mx0 = gam.m[0][0], mx1 = gam.m[0][1], mx2 = gam.m[0][2], mx3 = gam.m[0][3];
    float my0 = gam.m[1][0], my1 = gam.m[1][1], my2 = gam.m[1][2], my3 = gam.m[1][3];
    float mz0 = gam.m[2][0], mz1 = gam.m[2][1], mz2 = gam.m[2][2], mz3 = gam.m[2][3];

#ifdef _OPENMP
    int npt_total = bnx * bny * bnz;
    #pragma omp parallel for schedule(static) if(npt_total > 100000)
#endif
    for (int kk = 0; kk < bnz; kk++) {
        float bk_x = mx2 * kk + mx3;
        float bk_y = my2 * kk + my3;
        float bk_z = mz2 * kk + mz3;
        int out_base_k = kk * nxy_b;

        for (int jj = 0; jj < bny; jj++) {
            float bj_x = mx1 * jj + bk_x;
            float bj_y = my1 * jj + bk_y;
            float bj_z = mz1 * jj + bk_z;
            int out_base = out_base_k + jj * bnx;

            if (interp_code == AL_INTERP_LINEAR) {
                int i_lo = 0, i_hi = 0;
                al_safe_irange(bnx, bj_x, mx0, anx_safe,
                               bj_y, my0, any_safe, bj_z, mz0, anz_safe, &i_lo, &i_hi);

                /* Leading edge: bounds-checked */
                for (int ii = 0; ii < i_lo; ii++) {
                    int oob;
                    float v = al_interp_linear_checked(
                        mx0*ii+bj_x, my0*ii+bj_y, mz0*ii+bj_z,
                        nxh, nyh, nzh, anx1, any1, anz1, aim, anx, nxy_a, &oob);
                    avm[out_base+ii] = oob ? AL_OUTVAL : v;
                }

                /* Interior: no bounds check, no CLIP, no floorf */
                for (int ii = i_lo; ii < i_hi; ii++) {
                    float xx = mx0 * ii + bj_x;
                    float yy = my0 * ii + bj_y;
                    float zz = mz0 * ii + bj_z;
                    int ix0 = (int)xx; float fx = xx - ix0;
                    int jy0 = (int)yy; float fy = yy - jy0;
                    int kz0 = (int)zz; float fz = zz - kz0;
                    float w0 = 1.0f - fx, w1 = fx;
                    int a00 = ix0 + jy0*anx, a10 = ix0 + (jy0+1)*anx;
                    int k0off = kz0*nxy_a, k1off = (kz0+1)*nxy_a;
                    float f00 = w0*aim[a00+k0off] + w1*aim[a00+1+k0off];
                    float f10 = w0*aim[a10+k0off] + w1*aim[a10+1+k0off];
                    float f01 = w0*aim[a00+k1off] + w1*aim[a00+1+k1off];
                    float f11 = w0*aim[a10+k1off] + w1*aim[a10+1+k1off];
                    float wy0 = 1.0f - fy;
                    avm[out_base+ii] = (1.0f-fz)*(wy0*f00 + fy*f10)
                                           + fz *(wy0*f01 + fy*f11);
                }

                /* Trailing edge: bounds-checked */
                for (int ii = i_hi; ii < bnx; ii++) {
                    int oob;
                    float v = al_interp_linear_checked(
                        mx0*ii+bj_x, my0*ii+bj_y, mz0*ii+bj_z,
                        nxh, nyh, nzh, anx1, any1, anz1, aim, anx, nxy_a, &oob);
                    avm[out_base+ii] = oob ? AL_OUTVAL : v;
                }
            } else { /* AL_INTERP_CUBIC — split into interior (no bounds check) and border */
                /* For cubic: floor(coord)-1 >= 0 and floor(coord)+2 <= dim-1
                   So we need coord >= 1.0 and coord <= dim-2.0 (conservative) */
                float anx_csafe = (float)(anx - 3) + 0.999f;
                float any_csafe = (float)(any - 3) + 0.999f;
                float anz_csafe = (float)(anz - 3) + 0.999f;
                int c_lo = 0, c_hi = 0;
                int have_interior = (anx >= 4 && any >= 4 && anz >= 4) &&
                    al_safe_irange_ex(bnx, bj_x, mx0, anx_csafe,
                                      bj_y, my0, any_csafe, bj_z, mz0, anz_csafe,
                                      1.0f, &c_lo, &c_hi);

                /* Leading edge: bounds-checked */
                int c_lo_end = have_interior ? c_lo : bnx;
                for (int ii = 0; ii < c_lo_end; ii++) {
                    int oob;
                    float v = al_interp_cubic_checked(
                        mx0*ii+bj_x, my0*ii+bj_y, mz0*ii+bj_z,
                        nxh, nyh, nzh, anx1, any1, anz1, aim, anx, nxy_a, &oob);
                    avm[out_base+ii] = oob ? AL_OUTVAL : v;
                }

                /* Interior: no bounds check, no CLIP, no floorf */
                if (have_interior) {
                    for (int ii = c_lo; ii < c_hi; ii++) {
                        float xx = mx0 * ii + bj_x;
                        float yy = my0 * ii + bj_y;
                        float zz = mz0 * ii + bj_z;
                        int cix = (int)xx; float cfx = xx - cix;
                        int cjy = (int)yy; float cfy = yy - cjy;
                        int ckz = (int)zz; float cfz = zz - ckz;
                        float cwx_m1 = P_M1(cfx), cwx_00 = P_00(cfx), cwx_p1 = P_P1(cfx), cwx_p2 = P_P2(cfx);
                        int cix_m1 = cix-1, cix_p1 = cix+1, cix_p2 = cix+2;
                        int jm1 = (cjy-1)*anx, j00 = cjy*anx, jp1 = (cjy+1)*anx, jp2 = (cjy+2)*anx;
                        int km1 = (ckz-1)*nxy_a, k00 = ckz*nxy_a, kp1 = (ckz+1)*nxy_a, kp2 = (ckz+2)*nxy_a;
#define XINTC_FAST(joff,koff) cwx_m1*aim[cix_m1+joff+koff]+cwx_00*aim[cix+joff+koff]\
                              +cwx_p1*aim[cix_p1+joff+koff]+cwx_p2*aim[cix_p2+joff+koff]
                        float f_jm1_km1=XINTC_FAST(jm1,km1), f_j00_km1=XINTC_FAST(j00,km1);
                        float f_jp1_km1=XINTC_FAST(jp1,km1), f_jp2_km1=XINTC_FAST(jp2,km1);
                        float f_jm1_k00=XINTC_FAST(jm1,k00), f_j00_k00=XINTC_FAST(j00,k00);
                        float f_jp1_k00=XINTC_FAST(jp1,k00), f_jp2_k00=XINTC_FAST(jp2,k00);
                        float f_jm1_kp1=XINTC_FAST(jm1,kp1), f_j00_kp1=XINTC_FAST(j00,kp1);
                        float f_jp1_kp1=XINTC_FAST(jp1,kp1), f_jp2_kp1=XINTC_FAST(jp2,kp1);
                        float f_jm1_kp2=XINTC_FAST(jm1,kp2), f_j00_kp2=XINTC_FAST(j00,kp2);
                        float f_jp1_kp2=XINTC_FAST(jp1,kp2), f_jp2_kp2=XINTC_FAST(jp2,kp2);
                        float cwy_m1=P_M1(cfy), cwy_00=P_00(cfy), cwy_p1=P_P1(cfy), cwy_p2=P_P2(cfy);
                        float f_km1 = cwy_m1*f_jm1_km1+cwy_00*f_j00_km1+cwy_p1*f_jp1_km1+cwy_p2*f_jp2_km1;
                        float f_k00 = cwy_m1*f_jm1_k00+cwy_00*f_j00_k00+cwy_p1*f_jp1_k00+cwy_p2*f_jp2_k00;
                        float f_kp1 = cwy_m1*f_jm1_kp1+cwy_00*f_j00_kp1+cwy_p1*f_jp1_kp1+cwy_p2*f_jp2_kp1;
                        float f_kp2 = cwy_m1*f_jm1_kp2+cwy_00*f_j00_kp2+cwy_p1*f_jp1_kp2+cwy_p2*f_jp2_kp2;
                        float cwz_m1=P_M1(cfz), cwz_00=P_00(cfz), cwz_p1=P_P1(cfz), cwz_p2=P_P2(cfz);
                        avm[out_base+ii] = P_FACTOR * (cwz_m1*f_km1 + cwz_00*f_k00 + cwz_p1*f_kp1 + cwz_p2*f_kp2);
#undef XINTC_FAST
                    }

                    /* Trailing edge: bounds-checked */
                    for (int ii = c_hi; ii < bnx; ii++) {
                        int oob;
                        float v = al_interp_cubic_checked(
                            mx0*ii+bj_x, my0*ii+bj_y, mz0*ii+bj_z,
                            nxh, nyh, nzh, anx1, any1, anz1, aim, anx, nxy_a, &oob);
                        avm[out_base+ii] = oob ? AL_OUTVAL : v;
                    }
                }
            }
        }
    }
}

static inline void al_clip_values(float * restrict v, int n, float lo, float hi)
{
    for (int ii = 0; ii < n; ii++)
        if (v[ii] < lo) v[ii] = lo;
        else if (v[ii] > hi) v[ii] = hi;
}

/* Clamp a warped OUTPUT image to the source range to suppress cubic overshoot. This is the
   FALLBACK used only by al_scalar_warpone's GENERIC (non-affine) path — currently unreached,
   since every al_scalar_warpone caller passes al_wfunc_affine (which takes the fused path,
   clamped precisely in-FOV by al_clip_fused_cubic_infov) and cubic never routes here. Because
   this variant cannot distinguish in-FOV samples from the AL_OUTVAL(0) fill, it extends the
   clip range to include 0: this PRESERVES the fill (never clamps it to a positive source min)
   but, as a known limitation of the whole-array approach, it leaves in-FOV under/overshoot
   between 0 and the nearer source bound unclamped (e.g. a positive-only source's in-FOV
   undershoot in [0,lo) is not pulled up to lo). The fused path does not have this limitation.
   (The cost path clips via al_clip_values directly and is left untouched.) */
static void al_clip_to_source_range(float * restrict v, int n,
                                    const float * restrict src, int nsrc)
{
    if (v == NULL || src == NULL || n <= 0 || nsrc <= 0) return;
    float lo = src[0], hi = src[0];
    for (int ii = 1; ii < nsrc; ii++) {
        if (src[ii] < lo) lo = src[ii];
        else if (src[ii] > hi) hi = src[ii];
    }
    if (lo > AL_OUTVAL) lo = AL_OUTVAL;   /* never clamp the out-of-FOV fill upward */
    if (hi < AL_OUTVAL) hi = AL_OUTVAL;
    al_clip_values(v, n, lo, hi);
}

/*--- Get warped target values at control points (or all points) ---*/
/* Returns 0 on success (avm fully written), 1 if a buffer allocation failed
   (avm left untouched — caller must reject the cost, not score stale values). */
static int GA_get_warped_values(int nmpar, double *mpar, float *avm)
{
    (void)nmpar;
    int npar, ii, jj, kk, qq, pp, npp, mm, nx, ny, nxy, npt, nall, nper;
    float *wpar, v;
    float *imf = NULL, *jmf = NULL, *kmf = NULL;
    float *imw, *jmw, *kmw;
    float *aim;

    npar = gstup->wfunc_numpar;
    /* Reuse thread-local wpar buffer */
    if (tl_wpar_len < npar) {
        free(tl_wpar);
        tl_wpar = (float *)malloc(npar * sizeof(float));
        tl_wpar_len = tl_wpar ? npar : 0;  /* only update length on success */
    }
    wpar = tl_wpar;
    if (!wpar) return 1;
    nper = NPER;

    /* Load warping parameters */
    if (mpar != NULL) {
        for (ii = pp = 0; ii < npar; ii++) {
            if (gstup->wfunc_param[ii].fixed) {
                wpar[ii] = gstup->wfunc_param[ii].val_fixed;
            } else {
                v = (float)mpar[pp++];
                wpar[ii] = gstup->wfunc_param[ii].min
                          + gstup->wfunc_param[ii].siz * PRED01(v);
            }
        }
    } else {
        for (ii = 0; ii < npar; ii++)
            wpar[ii] = gstup->wfunc_param[ii].val_out;
    }

    /* Fast path: all voxels in order — use fused incremental warp+interp */
    if (gstup->im_ar == NULL && mpar != NULL &&
        gstup->wfunc == al_wfunc_affine &&
        (gstup->interp_code == AL_INTERP_LINEAR || gstup->interp_code == AL_INTERP_CUBIC)) {
        /* GA_setup_affine already applies aff_before/aff_after internally */
        mat44 gam = GA_setup_affine(npar, wpar);
        aim = (gstup->ajims != NULL) ? gstup->ajims : gstup->ajim;
        GA_warp_interp_fused(gam, aim, gstup->anx, gstup->any, gstup->anz,
                              gstup->bnx, gstup->bny, gstup->bnz,
                              gstup->interp_code, avm);
        /* Clip cubic overshoot to the original source range, matching the general path. */
        if (gstup->interp_code == AL_INTERP_CUBIC) {
            npt = gstup->bnx * gstup->bny * gstup->bnz;
            al_clip_values(avm, npt, gstup->ajmin, gstup->ajmax);
        }
        return 0;
    }

    /* Space for control points. Snapshot the mode once: allocation, use, and cleanup
       must agree even if this function is later reused in a less strictly serial host. */
    int alloc_ijk = (mpar == NULL || gstup->im_ar == NULL);
    if (alloc_ijk) {
        npt = gstup->bnx * gstup->bny * gstup->bnz;
        nall = (nper < npt) ? nper : npt;
        imf = (float *)calloc(nall, sizeof(float));
        jmf = (float *)calloc(nall, sizeof(float));
        kmf = (float *)calloc(nall, sizeof(float));
    } else {
        npt = gstup->npt_match;
        nall = (nper < npt) ? nper : npt;
    }

    /* Reuse thread-local warp buffer (holds imw,jmw,kmw contiguously) */
    if (tl_wbuf_len < nall * 3) {
        free(tl_wbuf);
        tl_wbuf = (float *)malloc(nall * 3 * sizeof(float));
        tl_wbuf_len = tl_wbuf ? nall * 3 : 0;  /* only update length on success */
    }
    if (!tl_wbuf || (alloc_ijk && (!imf || !jmf || !kmf))) {
        if (alloc_ijk) { free(imf); free(jmf); free(kmf); }
        return 1;
    }
    imw = tl_wbuf;
    jmw = tl_wbuf + nall;
    kmw = tl_wbuf + nall * 2;

    nx = gstup->bnx; ny = gstup->bny; nxy = nx * ny;

    /* Send parameters to warp function for setup */
    gstup->wfunc(npar, wpar, 0, NULL, NULL, NULL, NULL, NULL, NULL);

    /* Choose source image (smoothed if available) */
    aim = (gstup->ajims != NULL && mpar != NULL) ? gstup->ajims : gstup->ajim;

    /* Process in chunks */
    for (pp = 0; pp < npt; pp += nall) {
        npp = nall;
        if (npp > npt - pp) npp = npt - pp;

        if (alloc_ijk) {
            for (qq = 0; qq < npp; qq++) {
                mm = pp + qq;
                ii = mm % nx; kk = mm / nxy; jj = (mm - kk * nxy) / nx;
                imf[qq] = (float)ii; jmf[qq] = (float)jj; kmf[qq] = (float)kk;
            }
        } else {
            imf = gstup->im_ar + pp;
            jmf = gstup->jm_ar + pp;
            kmf = gstup->km_ar + pp;
        }

        /* Warp control points */
        gstup->wfunc(npar, NULL, npp, imf, jmf, kmf, imw, jmw, kmw);

        /* Interpolate source at warped locations */
        al_interp(aim, gstup->anx, gstup->any, gstup->anz,
                  gstup->interp_code, npp, imw, jmw, kmw, avm + pp);
    }

    /* imw/jmw/kmw and wpar are thread-local — not freed per call */
    if (alloc_ijk) {
        free(kmf); free(jmf); free(imf);
    }

    /* Clip interpolated values to original source data range (not CLEQWD range) */
    if (gstup->interp_code == AL_INTERP_CUBIC) {
        al_clip_values(avm, npt, gstup->ajmin, gstup->ajmax);
    }
    return 0;
}

/*--- Ensure blokset is created (avoids lazy-init race in parallel regions) ---*/
static void al_ensure_blokset(GA_setup *stup)
{
    if (stup->blokset != NULL) return;
    float rad = stup->blokrad, mrad;
    if (stup->smooth_code > 0 && stup->smooth_radius_base > 0.0f)
        rad = sqrtf(rad * rad + stup->smooth_radius_base * stup->smooth_radius_base);
    mrad = 1.2345f * (stup->base_di + stup->base_dj + stup->base_dk);
    if (rad < mrad) rad = mrad;
    stup->blokset = create_GA_BLOK_set(
        stup->bnx, stup->bny, stup->bnz,
        stup->base_di, stup->base_dj, stup->base_dk,
        stup->npt_match, stup->im_ar, stup->jm_ar, stup->km_ar,
        stup->bloktype, rad, stup->blokmin, 1.0f, 0);
}

/*--- Local Pearson correlation using BLOKs ---*/
/* Returns signed local correlation (always signed; LPA applies 1-|val| externally) */
static float GA_pearson_local(int npt, float *avm, float *bvm, float *wvm)
{
    (void)npt;
    GA_BLOK_set *gbs;
    int nblok, nelm, *elm, dd, ii, jj;
    float xv, yv, xy, xm, ym, vv, ww, ws, wss, pcor, wt, psum = 0.0f, pabs;

    al_ensure_blokset(gstup);
    if (gstup->blokset == NULL) return AL_BIGVAL;

    gbs   = gstup->blokset;
    nblok = gbs->num;
    if (nblok < 1) return AL_BIGVAL;

    for (wss = 0.0f, dd = 0; dd < nblok; dd++) {
        nelm = gbs->nelm[dd];
        if (nelm < 9) continue;
        elm = gbs->elm[dd];

        if (wvm == NULL) {
            xv = yv = xy = xm = ym = 0.0f; ws = 1.0f;
            for (ii = 0; ii < nelm; ii++) {
                jj = elm[ii];
                xm += avm[jj]; ym += bvm[jj];
            }
            xm /= nelm; ym /= nelm;
            for (ii = 0; ii < nelm; ii++) {
                jj = elm[ii];
                vv = avm[jj] - xm; ww = bvm[jj] - ym;
                xv += vv * vv; yv += ww * ww; xy += vv * ww;
            }
        } else {
            xv = yv = xy = xm = ym = ws = 0.0f;
            for (ii = 0; ii < nelm; ii++) {
                jj = elm[ii];
                wt = wvm[jj]; ws += wt;
                xm += avm[jj] * wt; ym += bvm[jj] * wt;
            }
            if (ws <= 0.0f) continue;
            xm /= ws; ym /= ws;
            for (ii = 0; ii < nelm; ii++) {
                jj = elm[ii];
                wt = wvm[jj]; vv = avm[jj] - xm; ww = bvm[jj] - ym;
                xv += wt * vv * vv; yv += wt * ww * ww; xy += wt * vv * ww;
            }
        }

        if (xv <= 0.0f || yv <= 0.0f) continue;
        pcor = xy / sqrtf(xv * yv);
        if (pcor > CMAX)  pcor = CMAX;
        else if (pcor < -CMAX) pcor = -CMAX;
        pcor = logf((1.0f + pcor) / (1.0f - pcor));  /* 2*arctanh */
        pabs = fabsf(pcor);
        psum += ws * pcor * pabs;  /* signed: emphasize large values */
        wss  += ws;
    }

    if (wss <= 0.0f) return AL_BIGVAL;
    return 0.25f * psum / wss;
}

#ifdef AL_LPC_MICHO
/*--- Overlap fraction between warped base mask and source mask ---*/
static float GA_get_warped_overlap_fraction(void)
{
    int npar, ii, jj, kk, qq, pp, nqq, mm, nx, nxy, nxt, nxyt, npt, nhit;
    float *imf, *jmf, *kmf, *imw, *jmw, *kmw, xx, yy, zz, nxh, nyh, nzh;
    unsigned char *bsar, *tgar;

    if (gstup->bmask == NULL || gstup->ajmask == NULL) return 1.0f;
    bsar = gstup->bmask;
    tgar = gstup->ajmask;

    npar = gstup->wfunc_numpar;
    npt  = gstup->bnx * gstup->bny * gstup->bnz;
    nx   = gstup->bnx; nxy = nx * gstup->bny;

    /* Count base mask voxels */
    for (nqq = pp = 0; pp < npt; pp++) if (bsar[pp]) nqq++;
    if (nqq == 0) return 1.0f;

    nxt  = gstup->anx; nxyt = nxt * gstup->any;
    nxh  = nxt - 0.501f; nyh = gstup->any - 0.501f; nzh = gstup->anz - 0.501f;

    imf = (float *)malloc(sizeof(float) * nqq);
    jmf = (float *)malloc(sizeof(float) * nqq);
    kmf = (float *)malloc(sizeof(float) * nqq);
    if (!imf || !jmf || !kmf) { free(imf); free(jmf); free(kmf); return 1.0f; }

    for (pp = qq = 0; pp < npt; pp++) {
        if (bsar[pp]) {
            ii = pp % nx; kk = pp / nxy; jj = (pp - kk * nxy) / nx;
            imf[qq] = (float)ii; jmf[qq] = (float)jj; kmf[qq] = (float)kk; qq++;
        }
    }

    imw = (float *)malloc(sizeof(float) * nqq);
    jmw = (float *)malloc(sizeof(float) * nqq);
    kmw = (float *)malloc(sizeof(float) * nqq);
    if (!imw || !jmw || !kmw) {
        free(imf); free(jmf); free(kmf);
        free(imw); free(jmw); free(kmw); return 1.0f;
    }

    gstup->wfunc(npar, NULL, nqq, imf, jmf, kmf, imw, jmw, kmw);
    free(kmf); free(jmf); free(imf);

    for (nhit = qq = 0; qq < nqq; qq++) {
        xx = imw[qq]; if (xx < -0.499f || xx > nxh) continue;
        yy = jmw[qq]; if (yy < -0.499f || yy > nyh) continue;
        zz = kmw[qq]; if (zz < -0.499f || zz > nzh) continue;
        ii = (int)(xx + 0.5f); jj = (int)(yy + 0.5f); kk = (int)(zz + 0.5f);
        if (tgar[ii + jj * nxt + kk * nxyt]) nhit++;
    }
    free(kmw); free(jmw); free(imw);

    xx = (float)nqq;
    int najmask = 0;
    for (pp = 0; pp < gstup->anx * gstup->any * gstup->anz; pp++)
        if (tgar[pp]) najmask++;
    yy = (float)najmask * gstup->adx * gstup->ady * gstup->adz /
         (gstup->bdx * gstup->bdy * gstup->bdz);
    float mn = (xx < yy) ? xx : yy;
    if (mn <= 0.0f) return 1.0f;
    return (float)nhit / mn;
}
#endif /* AL_LPC_MICHO */

/*--- Compute cost function ---*/
/* Weighted global Pearson correlation: returns negative correlation (for minimization).
   Two-pass algorithm: compute means, then variance/covariance.
   Auto-vectorized by the compiler with -ffast-math. */
static double GA_pearson_global(int npt, const float * restrict avm,
                                const float * restrict bvm,
                                const float * restrict wvm)
{
    double ws = 0.0, xm = 0.0, ym = 0.0;
    double xv = 0.0, yv = 0.0, xy = 0.0;
    int ii;

    if (wvm != NULL) {
        for (ii = 0; ii < npt; ii++) {
            double w = (double)wvm[ii];
            ws += w; xm += w * avm[ii]; ym += w * bvm[ii];
        }
        if (ws <= 0.0) return (double)AL_BIGVAL;
        xm /= ws; ym /= ws;
        for (ii = 0; ii < npt; ii++) {
            double w = (double)wvm[ii], dx = (double)avm[ii] - xm, dy = (double)bvm[ii] - ym;
            xv += w * dx * dx; yv += w * dy * dy; xy += w * dx * dy;
        }
    } else {
        for (ii = 0; ii < npt; ii++) { xm += avm[ii]; ym += bvm[ii]; }
        ws = npt; xm /= ws; ym /= ws;
        for (ii = 0; ii < npt; ii++) {
            double dx = (double)avm[ii] - xm, dy = (double)bvm[ii] - ym;
            xv += dx * dx; yv += dy * dy; xy += dx * dy;
        }
    }

    if (xv <= 0.0 || yv <= 0.0) return (double)AL_BIGVAL;
    return -(xy / sqrt(xv * yv));
}

static double GA_scalar_costfun(int meth, int npt,
                                float *avm, float *bvm, float *wvm)
{
    double val = 0.0;

    switch (meth) {
        case GA_MATCH_PEARSON_SCALAR: { /* Global Pearson (within-modality, fast) */
            val = GA_pearson_global(npt, avm, bvm, wvm);
        } break;

        case GA_MATCH_HELLINGER_SCALAR: { /* Hellinger only (fast, cross-modal) */
            float_quad hmc;
            hmc = al_helmicra(npt, gstup->hxc_bot, gstup->hxc_top, avm,
                                   gstup->hyc_bot, gstup->hyc_top, bvm, wvm,
                                   gstup->aj_topclip, gstup->bs_topclip);
            if (al_hist_oom) return (double)AL_BIGVAL; /* histogram OOM: reject, not a cost-0 "win" */
            val = -(double)hmc.a; /* aligned → hmc.a > 0 → val negative → good for min */
        } break;

        case GA_MATCH_PEARSON_LOCALS:  /* pure lpc (signed local Pearson) */
            val = (double)GA_pearson_local(npt, avm, bvm, wvm);
        break;

        case GA_MATCH_PEARSON_LOCALA:  /* pure lpa (absolute local Pearson) */
            val = (double)GA_pearson_local(npt, avm, bvm, wvm);
            val = 1.0 - fabs(val);
        break;

#ifdef AL_LPC_MICHO
        case GA_MATCH_LPC_MICHO_SCALAR: /* lpc+ZZ (signed + helper costs) */
        case GA_MATCH_LPA_MICHO_SCALAR: { /* lpa+ZZ: 1 - |lpc| + helpers */
            val = (double)GA_pearson_local(npt, avm, bvm, wvm);
            if (METH_IS_LPA(meth))
                val = 1.0 - fabs(val);  /* LPA: absolute, lower=better */
            if (gstup->micho_hel != 0.0 || gstup->micho_mi  != 0.0 ||
                gstup->micho_nmi != 0.0 || gstup->micho_crA != 0.0) {
                float_quad hmc;
                float ovv;
                hmc = al_helmicra(npt, gstup->hxc_bot, gstup->hxc_top, avm,
                                       gstup->hyc_bot, gstup->hyc_top, bvm, wvm,
                                       gstup->aj_topclip, gstup->bs_topclip);
                if (al_hist_oom) return (double)AL_BIGVAL; /* histogram OOM: reject this eval */
                val += -gstup->micho_hel * hmc.a - gstup->micho_mi * hmc.b
                       + gstup->micho_nmi * hmc.c + gstup->micho_crA * (1.0 - fabs(hmc.d));

                if (gstup->micho_ov != 0.0 && gstup->bmask != NULL && gstup->ajmask != NULL) {
                    ovv = GA_get_warped_overlap_fraction();
                    ovv = (9.95f - 10.0f * ovv);
                    if (ovv < 0.0f) ovv = 0.0f;
                    val += gstup->micho_ov * ovv * ovv;
                }
            }
        } break;
#endif

        default:
            return (double)AL_BIGVAL;
    }

    if (fabs(val) > 1.e+37) val = (double)AL_BIGVAL; /* magnitude guard: safe under -ffast-math */
    return val;
}

/*--- NEWUOA callback: evaluate cost from scaled parameters ---*/
#ifdef AL_PROFILE
static AL_TLOCAL double _fitter_warp_time = 0.0;
static AL_TLOCAL double _fitter_cost_time = 0.0;
static AL_TLOCAL int    _fitter_calls = 0;
#endif

static double GA_scalar_fitter(int npar, double *mpar)
{
    double val;
    float *avm, *bvm, *wvm;
    int npt = gstup->npt_match;

    /* Reuse thread-local buffer instead of per-call malloc/free */
    if (tl_avm_len < npt) {
        free(tl_avm);
        tl_avm = (float *)malloc(npt * sizeof(float));
        tl_avm_len = tl_avm ? npt : 0;  /* only update length on success */
    }
    avm = tl_avm;
    if (!avm) return (double)AL_BIGVAL;
    /* No memset needed: GA_get_warped_values writes every element
       (interpolated value or AL_OUTVAL=0 for out-of-bounds) */

#ifdef AL_PROFILE
    double _t0 = al_wtime();
#endif
    if (GA_get_warped_values(npar, mpar, avm) != 0)
        return (double)AL_BIGVAL;  /* buffer OOM: reject, don't score stale avm */
#ifdef AL_PROFILE
    double _t1 = al_wtime();
#endif

    bvm = gstup->bvm;
    wvm = gstup->wvm;

    /* AFNI set_2Dhist_xyclip (GA_scalar_fitter): on the first cost eval of a stage, compute the
       histogram edge-bin clips from the LIVE sample distribution (warped source avm + base bvm)
       and cache them, then clear the flag. This runs on the stage's first, SEQUENTIAL eval (the
       coarse ransetup center eval / the refinement warm-up eval) — never inside the parallel
       region — so the shared write is race-free and deterministic. Falls back to the image
       clipate defaults on a degenerate clip. Recomputing on the live samples is why AFNI keeps
       the cross-modal histogram resolution on the informative range as alignment improves. */
    if (gstup->need_hist_setup) {
        al_float_pair xc = al_clipate(npt, avm);
        al_float_pair yc = al_clipate(npt, bvm);
        if (xc.a < xc.b) { gstup->hxc_bot = xc.a; gstup->hxc_top = xc.b; }
        if (yc.a < yc.b) { gstup->hyc_bot = yc.a; gstup->hyc_top = yc.b; }
        gstup->need_hist_setup = 0;
        if (getenv("AL_VERB"))
            fprintf(stderr, "[AL_VERB hist] source clip %.4g .. %.4g; base clip %.4g .. %.4g\n",
                    gstup->hxc_bot, gstup->hxc_top, gstup->hyc_bot, gstup->hyc_top);
    }

    /* -dark_automask: zero the weight of any matched pair whose base or warped-source
       value is at that image's darkest value (background/pad). Folds into the existing
       per-point weight, so it costs one pass over the (already-materialized) point
       arrays and no extra work inside the cost functions. Fall back to the unmasked
       weights if too few pairs survive (gross early misalignment), so the optimizer
       still gets a usable gradient toward overlap. */
    if (gstup->do_dark_automask) {
        if (tl_weff_len < npt) {
            free(tl_weff);
            tl_weff = (float *)malloc(npt * sizeof(float));
            tl_weff_len = tl_weff ? npt : 0;
        }
        if (tl_weff) {
            float db = gstup->dark_base, dt = gstup->dark_targ;
            int nsurv = 0;
            for (int ii = 0; ii < npt; ii++) {
                if (avm[ii] > dt && bvm[ii] > db) {
                    tl_weff[ii] = wvm ? wvm[ii] : 1.0f; nsurv++;
                } else {
                    tl_weff[ii] = 0.0f;
                }
            }
            /* Keep the masked weights only if enough pairs survive for a stable cost;
               otherwise fall back to unmasked (gross early misalignment). Scale the
               floor with npt so the mask isn't silently a no-op on tiny point sets
               (npt < 64) while still requiring ~64 for large ones. */
            int min_surv = (npt < 256) ? (npt / 4) : 64;
            if (nsurv >= min_surv) wvm = tl_weff;
        }
    }

    val = GA_scalar_costfun(gstup->match_code, gstup->npt_match, avm, bvm, wvm);
#ifdef AL_PROFILE
    double _t2 = al_wtime();
    _fitter_warp_time += (_t1 - _t0);
    _fitter_cost_time += (_t2 - _t1);
    _fitter_calls++;
#endif
    return val;
}

#ifdef AL_PROFILE
/* Dump fitter profile stats.  Counters are AL_TLOCAL, so during OMP parallel
   regions each thread accumulates independently.  This dump collects from all
   threads (if OpenMP is active) or just the main thread otherwise. */
static void fitter_profile_dump(const char *label) {
    double total_warp = 0.0, total_cost = 0.0;
    int total_calls = 0;
#ifdef _OPENMP
    #pragma omp parallel reduction(+:total_warp,total_cost,total_calls)
    {
        total_warp  += _fitter_warp_time;
        total_cost  += _fitter_cost_time;
        total_calls += _fitter_calls;
        _fitter_warp_time = 0.0;
        _fitter_cost_time = 0.0;
        _fitter_calls = 0;
    }
#else
    total_warp  = _fitter_warp_time;
    total_cost  = _fitter_cost_time;
    total_calls = _fitter_calls;
    _fitter_warp_time = 0.0;
    _fitter_cost_time = 0.0;
    _fitter_calls = 0;
#endif
    fprintf(stderr, " [PROFILE] %s: %d calls, warp=%.3fs, cost=%.3fs (%.1f%% warp)\n",
            label, total_calls, total_warp, total_cost,
            (total_warp + total_cost > 0) ?
            100.0 * total_warp / (total_warp + total_cost) : 0.0);
}
#define FITTER_PROFILE_DUMP(label) fitter_profile_dump(label)
#else
#define FITTER_PROFILE_DUMP(label) ((void)0)
#endif

/*--- Count free params and compute siz ---*/
static void GA_param_setup(GA_setup *stup)
{
    int ii, qq;
    if (stup == NULL || stup->setup != AL_SMAGIC) return;
    for (ii = qq = 0; qq < stup->wfunc_numpar; qq++)
        if (!stup->wfunc_param[qq].fixed) ii++;
    stup->wfunc_numfree = ii;
    for (qq = 0; qq < stup->wfunc_numpar; qq++)
        stup->wfunc_param[qq].siz = stup->wfunc_param[qq].max
                                   - stup->wfunc_param[qq].min;
}

/*--- Setup matching: control points, base/weight values ---*/
static void al_scalar_setup(GA_setup *stup)
{
    int qq, ii, jj, kk, mm, nx, ny, nz, nxy, nmatch;
    float *bsar;
    unsigned char *mask;

    if (stup == NULL) return;
    stup->setup = 0;

    nx = stup->bnx; ny = stup->bny; nz = stup->bnz; nxy = nx * ny;
    mask = stup->bmask;

    /* Smooth base if needed */
    if (stup->smooth_code > 0 && stup->smooth_radius_base > 0.0f) {
        if (stup->bsims) free(stup->bsims);
        stup->bsims = al_smooth(stup->bsim, nx, ny, nz,
                                stup->bdx, stup->bdy, stup->bdz,
                                stup->smooth_radius_base);
        if (stup->bsims == NULL) return;
    } else {
        if (stup->bsims) { free(stup->bsims); stup->bsims = NULL; }
    }

    /* Smooth source if needed */
    if (stup->smooth_code > 0 && stup->smooth_radius_targ > 0.0f) {
        if (stup->ajims) free(stup->ajims);
        stup->ajims = al_smooth(stup->ajim, stup->anx, stup->any, stup->anz,
                                stup->adx, stup->ady, stup->adz,
                                stup->smooth_radius_targ);
        if (stup->ajims == NULL) return;
    } else {
        if (stup->ajims) { free(stup->ajims); stup->ajims = NULL; }
    }

    /* Source automask noise fill: replace non-mask voxels with random noise
       AFTER smoothing, matching AFNI's -source_automask behavior.
       This prevents background from contaminating local correlations. */
    if (stup->ajmask != NULL && stup->ajmask_ranfill && stup->aj_usiz > 0.0f) {
        float *af = (stup->ajims != NULL) ? stup->ajims : stup->ajim;
        unsigned char *mmm = stup->ajmask;
        int nvox = stup->anx * stup->any * stup->anz;
        float ubot = stup->aj_ubot, usiz = stup->aj_usiz;
        myunif_reset(1234567890);  /* deterministic noise */
        for (ii = 0; ii < nvox; ii++) {
            if (!mmm[ii]) {
                float u1 = myunif(), u2 = myunif();
                af[ii] = ubot + usiz * (u1 + u2);
            }
        }
    }

    /* The topclip-membership + live sample-clip (need_hist_setup) machinery is consumed ONLY by
       the histogram costs (Hellinger, and the -DAL_LPC_MICHO helper terms). For ls/lpc/lpa it is
       pure wasted work — a full-image quantile per axis, a per-stage sample clipate, and a
       discarded sequential warm-up eval. Gate it on the cost so those paths skip it. The default
       (Hellinger) is byte-for-byte unchanged (hist_cost == 1). */
    int hist_cost = (stup->match_code == GA_MATCH_HELLINGER_SCALAR);
#ifdef AL_LPC_MICHO
    hist_cost = hist_cost || stup->match_code == GA_MATCH_LPC_MICHO_SCALAR
                          || stup->match_code == GA_MATCH_LPA_MICHO_SCALAR;
#endif

    /* Get min/max and CLEQWD clip levels for source image */
    {
        float *src = (stup->ajims != NULL) ? stup->ajims : stup->ajim;
        int nvox = stup->anx * stup->any * stup->anz;
        stup->ajbot = src[0]; stup->ajtop = src[0];
        for (ii = 1; ii < nvox; ii++) {
            if (src[ii] < stup->ajbot) stup->ajbot = src[ii];
            if (src[ii] > stup->ajtop) stup->ajtop = src[ii];
        }
        stup->ajmin = stup->ajbot;  /* preserve original data range */
        stup->ajmax = stup->ajtop;
        stup->ajclip = stup->ajtop;
        /* AFNI mri_topclip = MIN(3.11*THD_cliplevel(0.511), max) for non-negative images:
           the histogram-MEMBERSHIP top (pairs with a value above this are DROPPED, matching
           AFNI). Negative images (e.g. CT) keep the full range. Distinct from the clipate
           edge-bin clip below — collapsing the two lets bright cross-modal outliers (T2 CSF/
           fat) bias the fit toward shrinking overlap (see AGENTS.md). */
        stup->aj_topclip = stup->ajmax;
        if (hist_cost && stup->ajmin >= 0.0f) {
            float tc = 3.11f * al_cliplevel(nvox, src, 0.511f);
            if (tc < stup->aj_topclip) stup->aj_topclip = tc;
        }
        /* CLEQWD edge-bin clips are consumed only by histogram costs. */
        if (hist_cost) {
            al_float_pair cp = al_clipate(nvox, src);
            if (cp.a < cp.b) { stup->ajbot = cp.a; stup->ajclip = cp.b; }
        }
    }

    /* Get min/max and CLEQWD clip levels for base image */
    {
        float *bas = (stup->bsims != NULL) ? stup->bsims : stup->bsim;
        int nvox = nx * ny * nz;
        stup->bsbot = bas[0]; stup->bstop = bas[0];
        for (ii = 1; ii < nvox; ii++) {
            if (bas[ii] < stup->bsbot) stup->bsbot = bas[ii];
            if (bas[ii] > stup->bstop) stup->bstop = bas[ii];
        }
        stup->bsclip = stup->bstop;
        /* AFNI mri_topclip membership top for the base (see source block above). */
        stup->bs_topclip = stup->bstop;
        if (hist_cost && stup->bsbot >= 0.0f) {
            float tc = 3.11f * al_cliplevel(nvox, bas, 0.511f);
            if (tc < stup->bs_topclip) stup->bs_topclip = tc;
        }
        /* CLEQWD edge-bin clips are consumed only by histogram costs. */
        if (hist_cost) {
            al_float_pair cp = al_clipate(nvox, bas);
            if (cp.a < cp.b) { stup->bsbot = cp.a; stup->bsclip = cp.b; }
        }
    }
    /* Edge-bin clips default to the image clipate (the fallback); the stage's first cost eval
       refreshes them from the live sample distribution (AFNI set_2Dhist_xyclip). */
    stup->hxc_bot = stup->ajbot; stup->hxc_top = stup->ajclip;
    stup->hyc_bot = stup->bsbot; stup->hyc_top = stup->bsclip;
    stup->need_hist_setup = hist_cost;   /* refresh live sample clips only for histogram costs */

    /* Determine number of matching points */
    nmatch = stup->npt_match;
    if (nmatch <= 9 || nmatch > nx * ny * nz)
        nmatch = nx * ny * nz;
    if (stup->nmask > 0 && nmatch > stup->nmask)
        nmatch = stup->nmask;
    stup->npt_match = nmatch;
    if (getenv("AL_VERB"))
        fprintf(stderr, "[AL_VERB setup] smooth=%.2f npt_match=%d nmask=%d | "
                "src range %.4g..%.4g clip %.4g..%.4g topclip %.4g | "
                "base range %.4g..%.4g clip %.4g..%.4g topclip %.4g\n",
                stup->smooth_radius_base, stup->npt_match, stup->nmask,
                stup->ajmin, stup->ajmax, stup->ajbot, stup->ajclip, stup->aj_topclip,
                stup->bsbot, stup->bstop, stup->bsbot, stup->bsclip, stup->bs_topclip);

    /* Free old control point arrays */
    if (stup->im_ar) { free(stup->im_ar); stup->im_ar = NULL; }
    if (stup->jm_ar) { free(stup->jm_ar); stup->jm_ar = NULL; }
    if (stup->km_ar) { free(stup->km_ar); stup->km_ar = NULL; }

    int use_all = 0;
    if (nmatch >= nx * ny * nz || (stup->nmask > 0 && nmatch >= stup->nmask))
        use_all = 1;

    if (use_all && stup->nmask == 0) {
        /* All points, no mask - no index arrays needed, use all voxels in order */
        /* im_ar stays NULL, signaling "use all" */
    } else if (use_all && stup->nmask > 0) {
        /* All mask points */
        stup->im_ar = (float *)malloc(sizeof(float) * nmatch);
        stup->jm_ar = (float *)malloc(sizeof(float) * nmatch);
        stup->km_ar = (float *)malloc(sizeof(float) * nmatch);
        if (!stup->im_ar || !stup->jm_ar || !stup->km_ar) {
            free(stup->im_ar); stup->im_ar = NULL;
            free(stup->jm_ar); stup->jm_ar = NULL;
            free(stup->km_ar); stup->km_ar = NULL;
            return;
        }
        { int nvox_t = nx * ny * nz;
          for (mm = qq = 0; qq < nmatch && mm < nvox_t; mm++) {
            if (mask == NULL || mask[mm]) {
                ii = mm % nx; kk = mm / nxy; jj = (mm - kk * nxy) / nx;
                stup->im_ar[qq] = (float)ii;
                stup->jm_ar[qq] = (float)jj;
                stup->km_ar[qq] = (float)kk;
                qq++;
            }
          }
          if (qq != nmatch) {
              free(stup->im_ar); stup->im_ar = NULL;
              free(stup->jm_ar); stup->jm_ar = NULL;
              free(stup->km_ar); stup->km_ar = NULL;
              return;  /* inconsistent mask count: leave setup invalid and fail closed */
          }
        }
    } else {
        /* Subset of points */
        int nvox = nx * ny * nz;
        int dm = ga_find_relprime_fixed(nvox);
        stup->im_ar = (float *)malloc(sizeof(float) * nmatch);
        stup->jm_ar = (float *)malloc(sizeof(float) * nmatch);
        stup->km_ar = (float *)malloc(sizeof(float) * nmatch);
        if (!stup->im_ar || !stup->jm_ar || !stup->km_ar) {
            free(stup->im_ar); stup->im_ar = NULL;
            free(stup->jm_ar); stup->jm_ar = NULL;
            free(stup->km_ar); stup->km_ar = NULL;
            return;
        }
        mm = (nx / 2) + (ny / 2) * nx + (nz / 2) * nxy;
        for (qq = 0; qq < nmatch; mm = (mm + dm) % nvox)
            if (mask == NULL || mask[mm]) {
                ii = mm % nx; kk = mm / nxy; jj = (mm - kk * nxy) / nx;
                stup->im_ar[qq] = (float)ii;
                stup->jm_ar[qq] = (float)jj;
                stup->km_ar[qq] = (float)kk;
                qq++;
            }
    }

    /* Extract base values at control points */
    if (stup->bvm) { free(stup->bvm); stup->bvm = NULL; }
    if (stup->wvm) { free(stup->wvm); stup->wvm = NULL; }

    bsar = (stup->bsims != NULL) ? stup->bsims : stup->bsim;
    stup->bvm = (float *)calloc(nmatch, sizeof(float));
    if (stup->bwght != NULL)
        stup->wvm = (float *)calloc(nmatch, sizeof(float));
    if (!stup->bvm || (stup->bwght != NULL && !stup->wvm)) {
        free(stup->bvm); stup->bvm = NULL;
        free(stup->wvm); stup->wvm = NULL;
        return;  /* OOM: leave stup->setup = 0 so callers fail closed (checked below) */
    }

    if (stup->im_ar == NULL) {
        memcpy(stup->bvm, bsar, sizeof(float) * nmatch);
        if (stup->wvm && stup->bwght)
            memcpy(stup->wvm, stup->bwght, sizeof(float) * nmatch);
    } else {
        for (qq = 0; qq < nmatch; qq++) {
            int rr = (int)(stup->im_ar[qq] + stup->jm_ar[qq] * nx + stup->km_ar[qq] * nxy);
            stup->bvm[qq] = bsar[rr];
        }
        if (stup->bwght != NULL && stup->wvm != NULL) {
            for (qq = 0; qq < nmatch; qq++) {
                int rr = (int)(stup->im_ar[qq] + stup->jm_ar[qq] * nx + stup->km_ar[qq] * nxy);
                stup->wvm[qq] = stup->bwght[rr];
            }
        }
    }

    /* Destroy old BLOK set (will be recreated when needed) */
    if (stup->blokset) { free_GA_BLOK_set(stup->blokset); stup->blokset = NULL; }

    stup->setup = AL_SMAGIC;
}

/*--- Optimize: run powell_newuoa on the current setup ---*/
static int al_scalar_optim(GA_setup *stup, double rstart, double rend, int nstep)
{
    double *wpar;
    int ii, qq, nfunc;

    if (stup == NULL || stup->setup != AL_SMAGIC) return -1;
    GA_param_setup(stup);
    if (stup->wfunc_numfree <= 0) return -2;

    wpar = (double *)calloc(stup->wfunc_numfree, sizeof(double));
    if (wpar == NULL) return -3;   /* OOM: fail closed rather than deref NULL */
    for (ii = qq = 0; qq < stup->wfunc_numpar; qq++) {
        if (!stup->wfunc_param[qq].fixed) {
            wpar[ii] = (stup->wfunc_param[qq].val_init - stup->wfunc_param[qq].min)
                       / stup->wfunc_param[qq].siz;
            if (wpar[ii] < 0.0 || wpar[ii] > 1.0) wpar[ii] = PRED01(wpar[ii]);
            ii++;
        }
    }

    gstup = stup;
    if (nstep <= 4 * stup->wfunc_numfree + 5) nstep = 6666;
    if (rstart > 0.2)  rstart = 0.2;
    else if (rstart <= 0.0) rstart = 0.1;
    if (rend >= 0.9 * rstart || rend <= 0.0) rend = 0.0666 * rstart;

    nfunc = powell_newuoa(stup->wfunc_numfree, wpar, rstart, rend, nstep, GA_scalar_fitter);
    if (nfunc < 0) { free(wpar); return nfunc; }  /* optimizer OOM: don't accept unoptimized params */
    double final_cost = GA_scalar_fitter(stup->wfunc_numfree, wpar);
    if (!(final_cost < (double)AL_BIGVAL)) {
        free(wpar);
        return -4;  /* invalid/failed final evaluation: do not report an unscored fit */
    }
    stup->vbest = (float)final_cost;

    /* Copy results back */
    for (ii = qq = 0; qq < stup->wfunc_numpar; qq++) {
        if (stup->wfunc_param[qq].fixed) {
            stup->wfunc_param[qq].val_out = stup->wfunc_param[qq].val_fixed;
        } else {
            stup->wfunc_param[qq].val_out = stup->wfunc_param[qq].min
                                           + stup->wfunc_param[qq].siz * PRED01(wpar[ii]);
            ii++;
        }
    }
    free(wpar);
    return nfunc;
}

/*--- Random startup search (coarse pass): returns nonzero on allocation failure ---*/
static int al_scalar_ransetup(GA_setup *stup, int nrand)
{
#define NKEEP (3*PARAM_MAXTRIAL+1)
    double val, vbest, *bpar;
    double *kpar[NKEEP] = {0}, kval[NKEEP];
    double *all_wpar = NULL, *tvals = NULL, *tpars = NULL;
    int *all_isrand = NULL, *tisrand = NULL;
    int ii, qq, nfr, kk, jj, ngrid, ngtot, ngood, twof, maxstep;
    int ival[NKEEP], rval[NKEEP];
    float fval[NKEEP];
    int icod, oom = 0;

    if (stup == NULL || stup->setup != AL_SMAGIC) return 0;
    if (nrand < NKEEP) nrand = NKEEP + 13;

    GA_param_setup(stup); gstup = stup;
    if (stup->wfunc_numfree <= 0) return 0;

    nfr = stup->wfunc_numfree;
    switch (nfr) {
        case 1: ngrid = 9; break;
        case 2: ngrid = 8; break;
        case 3: ngrid = 6; break;
        case 4: ngrid = 4; break;
        case 5: ngrid = 3; break;
        case 6: ngrid = 2; break;
        default: ngrid = 1; break;
    }
    for (ngtot = 1, qq = 0; qq < nfr; qq++) ngtot *= ngrid;
    if (nfr < 4) nrand *= 2;

    /* Save and set interp to linear for coarse pass */
    icod = stup->interp_code;
    stup->interp_code = AL_INTERP_LINEAR;

    for (kk = 0; kk < NKEEP; kk++) {
        kpar[kk] = (double *)calloc(nfr, sizeof(double));
        if (!kpar[kk]) { oom = 1; goto ran_cleanup; }
    }

    /* Evaluate center of parameter space */
    for (qq = 0; qq < nfr; qq++) kpar[0][qq] = 0.5;
    val = GA_scalar_fitter(nfr, kpar[0]);
    kval[0] = val; rval[0] = 0;
    for (kk = 1; kk < NKEEP; kk++) { rval[kk] = 0; kval[kk] = AL_BIGVAL; }

    twof = 1 << nfr;
    myunif_reset(3456789);

    int ntotal = ngtot + nrand;
    fprintf(stderr, " + Coarse search: %d grid + %d random trials\n", ngtot, nrand);

    /* Pre-generate all starting parameter sets (sequential — uses RNG) */
    all_wpar = (double *)malloc((size_t)ntotal * nfr * sizeof(double));
    all_isrand = (int *)malloc((size_t)ntotal * sizeof(int));
    if (!all_wpar || !all_isrand) { oom = 1; goto ran_cleanup; }
    for (ii = 0; ii < ntotal; ii++) {
        double *wp = all_wpar + ii * nfr;
        if (ii < ngtot) {
            double kp;
            double gval = 0.5 / (ngrid + 1.0);
            int s = ii;
            for (qq = 0; qq < nfr; qq++) {
                kk = s % ngrid; s = s / ngrid;
                kp = (kk == 0) ? 0.5 : (kk + 1.0);
                wp[qq] = 0.5 + kp * gval;
            }
            all_isrand[ii] = 0;
        } else {
            for (qq = 0; qq < nfr; qq++) wp[qq] = 0.5 * (1.05 + 0.90 * myunif());
            all_isrand[ii] = 1;
        }
    }

    al_ensure_blokset(stup);

    /* Evaluate all grid+random × reflections in parallel */
    {
        long ntasks = (long)ntotal * twof;
        tvals = (double *)malloc((size_t)ntasks * sizeof(double));
        tpars = (double *)malloc((size_t)ntasks * nfr * sizeof(double));
        tisrand = (int *)malloc((size_t)ntasks * sizeof(int));
        if (!tvals || !tpars || !tisrand) { oom = 1; goto ran_cleanup; }

        /* Pre-generate all reflected parameter sets */
        for (long idx = 0; idx < ntasks; idx++) {
            int ti = (int)(idx / twof);
            int ts = (int)(idx % twof);
            double *wp = all_wpar + ti * nfr;
            double *sp = tpars + idx * nfr;
            for (qq = 0; qq < nfr; qq++)
                sp[qq] = (ts & (1 << qq)) ? 1.0 - wp[qq] : wp[qq];
            tisrand[idx] = all_isrand[ti];
        }

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 16)
#endif
        for (long idx = 0; idx < ntasks; idx++) {
            double *sp = tpars + idx * nfr;
            tvals[idx] = GA_scalar_fitter(nfr, sp);
        }

        /* Insert results into top-K (sequential — fast) */
        for (long idx = 0; idx < ntasks; idx++) {
            val = tvals[idx];
            double *sp = tpars + idx * nfr;
            for (kk = 0; kk < NKEEP; kk++) {
                if (val < kval[kk]) {
                    for (jj = NKEEP - 2; jj >= kk; jj--) {
                        memcpy(kpar[jj + 1], kpar[jj], sizeof(double) * nfr);
                        kval[jj + 1] = kval[jj]; rval[jj + 1] = rval[jj];
                    }
                    memcpy(kpar[kk], sp, sizeof(double) * nfr);
                    kval[kk] = val; rval[kk] = tisrand[idx];
                    break;
                }
            }
        }

        free(tisrand); tisrand = NULL;
        free(tpars); tpars = NULL;
        free(tvals); tvals = NULL;
    }
    free(all_isrand); all_isrand = NULL;
    free(all_wpar); all_wpar = NULL;

    for (ngood = kk = 0; kk < NKEEP && kval[kk] < AL_BIGVAL; kk++, ngood++) ;
    if (ngood < 1) {
        // Fail closed: EVERY coarse candidate (including the near-identity center) scored
        // AL_BIGVAL — from a cost-path OOM or degenerate input (constant/empty base, no
        // surviving bloks). Returning success here let al_register proceed with only the
        // identity transform and print "Registration complete" on an unregistered result —
        // a correctness bug and, for -deface, a privacy failure (identity pulls the mask
        // onto the wrong grid, leaving faces intact). Signal the caller to abort instead.
        fprintf(stderr, "allineate: no good starting locations found\n");
        oom = 1;
        goto ran_cleanup;
    }

    /* Make sure all in 0..1 range */
    for (kk = 0; kk < ngood; kk++)
        for (ii = 0; ii < nfr; ii++) kpar[kk][ii] = PRED01(kpar[kk][ii]);

    /* Little optimization on each set */
    vbest = AL_BIGVAL; jj = 0;
    stup->interp_code = AL_INTERP_LINEAR;
    maxstep = 11 * nfr + 17;
    fprintf(stderr, " + Refining %d candidate parameter sets\n", ngood);

    al_ensure_blokset(stup);

    /* mfac/afac are thread-local: the caller's powell_set_mfac() ran on the main
       thread only, so capture it and re-apply per worker (else workers use the
       default/stale factors → thread-count-dependent npt). See al_register. */
    float rs_mfac, rs_afac; powell_get_mfac(&rs_mfac, &rs_afac);
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (kk = 0; kk < ngood; kk++) {
        if (kval[kk] >= AL_BIGVAL) continue;
        powell_set_mfac(rs_mfac, rs_afac);
        int prc = powell_newuoa(nfr, kpar[kk], 0.05, 0.001, maxstep, GA_scalar_fitter);
        kval[kk] = (prc < 0) ? AL_BIGVAL : GA_scalar_fitter(nfr, kpar[kk]);  /* OOM -> not selectable */
    }
    for (kk = 0; kk < ngood; kk++) {
        if (kval[kk] < vbest) { vbest = kval[kk]; jj = kk; }
    }
    stup->vbest = (float)vbest;

    for (kk = 0; kk < ngood; kk++)
        for (ii = 0; ii < nfr; ii++) kpar[kk][ii] = PRED01(kpar[kk][ii]);

    /* Save best into val_init/val_out */
    bpar = kpar[jj];
    for (ii = qq = 0; qq < stup->wfunc_numpar; qq++) {
        if (!stup->wfunc_param[qq].fixed) {
            stup->wfunc_param[qq].val_init = stup->wfunc_param[qq].min
                                            + stup->wfunc_param[qq].siz * bpar[ii];
            ii++;
        }
        stup->wfunc_param[qq].val_out = stup->wfunc_param[qq].val_init;
    }

    /* Save trial parameter sets (sorted by cost, excluding those too close) */
    for (ii = 0; ii < ngood; ii++) { fval[ii] = (float)kval[ii]; ival[ii] = ii; }
    qsort_floatint(ngood, fval, ival);

    /* Save best into trial #0 */
    for (qq = 0; qq < stup->wfunc_numpar; qq++) {
        if (!stup->wfunc_param[qq].fixed)
            stup->wfunc_param[qq].val_trial[0] = stup->wfunc_param[qq].val_init;
        else
            stup->wfunc_param[qq].val_trial[0] = stup->wfunc_param[qq].val_fixed;
    }

    int nt = 1;
    for (jj = 1; jj < ngood && nt < PARAM_MAXTRIAL; jj++) {
        double *qpar = kpar[ival[jj]];
        int skip = 0;
        for (kk = 0; kk < jj && !skip; kk++) {
            double *cpar = kpar[ival[kk]];
            double dist = 0.0;
            for (ii = 0; ii < nfr; ii++) {
                double d = fabs(qpar[ii] - cpar[ii]);
                if (d > dist) dist = d;
            }
            if (dist < 0.05) skip = 1;
        }
        if (skip) continue;
        for (ii = qq = 0; qq < stup->wfunc_numpar; qq++) {
            if (!stup->wfunc_param[qq].fixed) {
                stup->wfunc_param[qq].val_trial[nt] = stup->wfunc_param[qq].min
                    + stup->wfunc_param[qq].siz * qpar[ii];
                ii++;
            } else {
                stup->wfunc_param[qq].val_trial[nt] = stup->wfunc_param[qq].val_fixed;
            }
        }
        nt++;
    }
    stup->wfunc_ntrial = nt;

ran_cleanup:
    free(tisrand);
    free(tpars);
    free(tvals);
    free(all_isrand);
    free(all_wpar);
    for (kk = 0; kk < NKEEP; kk++) free(kpar[kk]);
    stup->interp_code = icod;
    if (oom)
        fprintf(stderr, "allineate: coarse startup search failed\n"); // OOM or no usable start
    return oom;
#undef NKEEP
}

/* In-FOV cubic clamp for the FUSED final warp (al_scalar_warpone). GA_warp_interp_fused
   wrote AL_OUTVAL for out-of-FOV output voxels and interpolated values for in-FOV ones;
   cubic can overshoot the source range (ringing). This clamps ONLY the in-FOV voxels to
   the full source range [lo,hi] — recomputing the exact FOV predicate the kernel used from
   `gam` — and leaves the AL_OUTVAL fill untouched. So (unlike al_clip_to_source_range's
   range-extension fallback) in-FOV under/overshoot toward zero IS clamped while the fill is
   preserved. Touches neither GA_warp_interp_fused nor the cost path (whose clip is separate).
   `gam` is the base-index -> source-index affine passed to the kernel. */
static void al_clip_fused_cubic_infov(float * restrict war, mat44 gam,
                                      const float * restrict src, int anx, int any, int anz,
                                      int bnx, int bny, int bnz)
{
    size_t nsrc = (size_t)anx * any * anz;
    if (!war || !src || nsrc < 1) return;
    float lo = src[0], hi = src[0];
    for (size_t i = 1; i < nsrc; i++) { if (src[i] < lo) lo = src[i]; else if (src[i] > hi) hi = src[i]; }
    float nxh = anx - 0.501f, nyh = any - 0.501f, nzh = anz - 0.501f;   /* == kernel FOV bounds */
    float mx0=gam.m[0][0], mx1=gam.m[0][1], mx2=gam.m[0][2], mx3=gam.m[0][3];
    float my0=gam.m[1][0], my1=gam.m[1][1], my2=gam.m[1][2], my3=gam.m[1][3];
    float mz0=gam.m[2][0], mz1=gam.m[2][1], mz2=gam.m[2][2], mz3=gam.m[2][3];
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) if((size_t)bnx*bny*bnz > 100000)
#endif
    for (int kk = 0; kk < bnz; kk++) {
        for (int jj = 0; jj < bny; jj++) {
            size_t ob = ((size_t)kk * bny + jj) * bnx;
            float bx = mx1*jj + mx2*kk + mx3;
            float by = my1*jj + my2*kk + my3;
            float bz = mz1*jj + mz2*kk + mz3;
            for (int ii = 0; ii < bnx; ii++) {
                float xx = mx0*ii + bx, yy = my0*ii + by, zz = mz0*ii + bz;
                if (xx < -0.499f || xx > nxh || yy < -0.499f || yy > nyh ||
                    zz < -0.499f || zz > nzh) continue;   /* out-of-FOV: keep AL_OUTVAL fill */
                float v = war[ob+ii];
                if (v < lo) war[ob+ii] = lo; else if (v > hi) war[ob+ii] = hi;
            }
        }
    }
}

/*--- Warp an entire image to base grid ---*/
static float *al_scalar_warpone(int npar, float *wpar, GA_warpfunc wfunc,
                                float *imtarg, int tnx, int tny, int tnz,
                                int nnx, int nny, int nnz, int icode)
{
    int ii, jj, kk, qq, pp, npp, mm, nx, ny, nxy, nz, npt, nall, nper;
    float *imf, *jmf, *kmf, *imw, *jmw, *kmw;
    float *war;

    nper = NPER;

    nx = nnx; ny = nny; nz = nnz; nxy = nx * ny; npt = nxy * nz;
    war = (float *)calloc(npt, sizeof(float));
    if (!war) return NULL;

    if (wfunc == al_wfunc_affine &&
        (icode == AL_INTERP_LINEAR || icode == AL_INTERP_CUBIC)) {
        mat44 gam = GA_setup_affine(npar, wpar);
        GA_warp_interp_fused(gam, imtarg, tnx, tny, tnz,
                              nnx, nny, nnz, icode, war);
        if (icode == AL_INTERP_CUBIC)
            al_clip_fused_cubic_infov(war, gam, imtarg, tnx, tny, tnz, nnx, nny, nnz);
        return war;
    }

    nall = (nper < npt) ? nper : npt;
    imf = (float *)calloc(nall, sizeof(float));
    jmf = (float *)calloc(nall, sizeof(float));
    kmf = (float *)calloc(nall, sizeof(float));
    imw = (float *)calloc(nall, sizeof(float));
    jmw = (float *)calloc(nall, sizeof(float));
    kmw = (float *)calloc(nall, sizeof(float));
    if (!imf || !jmf || !kmf || !imw || !jmw || !kmw) {
        free(imf); free(jmf); free(kmf); free(imw); free(jmw); free(kmw);
        free(war); return NULL;
    }

    /* outval is always 0.0f (AL_OUTVAL) */
    /* Send parameters to warp function for generic chunked warps */
    wfunc(npar, wpar, 0, NULL, NULL, NULL, NULL, NULL, NULL);

    for (pp = 0; pp < npt; pp += nall) {
        npp = nall;
        if (npp > npt - pp) npp = npt - pp;

        for (qq = 0; qq < npp; qq++) {
            mm = pp + qq;
            ii = mm % nx; kk = mm / nxy; jj = (mm - kk * nxy) / nx;
            imf[qq] = (float)ii; jmf[qq] = (float)jj; kmf[qq] = (float)kk;
        }

        wfunc(npar, NULL, npp, imf, jmf, kmf, imw, jmw, kmw);
        al_interp(imtarg, tnx, tny, tnz, icode, npp, imw, jmw, kmw, war + pp);
    }

    /* (outval restore removed - AL_OUTVAL is a constant) */

    /* Clip to source range */
    if (icode == AL_INTERP_CUBIC)
        al_clip_to_source_range(war, npt, imtarg, tnx * tny * tnz);

    free(kmw); free(jmw); free(imw);
    free(kmf); free(jmf); free(imf);
    return war;
}

/*==========================================================================*/
/*============== SECTION 9: AUTOWEIGHT AND AUTOMASK =======================*/
/*==========================================================================*/

/* Compute autoweight from base image */
static float *al_autoweight(float *basim, int nx, int ny, int nz,
                            float dx, float dy, float dz)
{
    if (!basim || nx < 1 || ny < 1 || nz < 1) return NULL;
    int nxyz = nx * ny * nz, ii, jj, kk, ff;
    int xfade, yfade, zfade;
    float *wf, clip, mx;
    float sigma;

    wf = (float *)malloc(sizeof(float) * nxyz);
    if (!wf) return NULL;
    for (ii = 0; ii < nxyz; ii++) wf[ii] = fabsf(basim[ii]);

    /* Zero out edges */
#undef  WV
#define WV(i,j,k) wf[(i)+(j)*nx+(k)*(nx*ny)]
    xfade = (int)(0.05 * nx + 3.0); if (5 * xfade >= nx) xfade = (nx - 1) / 5;
    yfade = (int)(0.05 * ny + 3.0); if (5 * yfade >= ny) yfade = (ny - 1) / 5;
    zfade = (int)(0.05 * nz + 3.0); if (5 * zfade >= nz) zfade = (nz - 1) / 5;
    for (jj = 0; jj < ny; jj++)
        for (ii = 0; ii < nx; ii++)
            for (ff = 0; ff < zfade; ff++) WV(ii, jj, ff) = WV(ii, jj, nz - 1 - ff) = 0.0f;
    for (kk = 0; kk < nz; kk++)
        for (jj = 0; jj < ny; jj++)
            for (ff = 0; ff < xfade; ff++) WV(ff, jj, kk) = WV(nx - 1 - ff, jj, kk) = 0.0f;
    for (kk = 0; kk < nz; kk++)
        for (ii = 0; ii < nx; ii++)
            for (ff = 0; ff < yfade; ff++) WV(ii, ff, kk) = WV(ii, ny - 1 - ff, kk) = 0.0f;

    /* Squash large values: clip at 3× median of nonzero voxels (matching AFNI) */
    {
        int npos = 0;
        for (ii = 0; ii < nxyz; ii++) if (wf[ii] > 0.0f) npos++;
        if (npos > 0) {
            float *tmp = (float *)malloc(sizeof(float) * npos);
            if (tmp) {
                int tt = 0;
                for (ii = 0; ii < nxyz; ii++) if (wf[ii] > 0.0f) tmp[tt++] = wf[ii];
                float median = al_select_rank(tmp, npos, npos / 2);  /* O(n) select */
                free(tmp);
                clip = 3.0f * median;
                for (ii = 0; ii < nxyz; ii++) if (wf[ii] > clip) wf[ii] = clip;
            }
        }
    }

    /* Gaussian blur */
    sigma = 2.25f * (dx + dy + dz) / 3.0f;
    if (gaussian_blur_3d(wf, nx, ny, nz, dx, dy, dz, sigma)) {
        free(wf);
        return NULL;
    }

    /* Threshold */
    mx = 0.0f;
    for (ii = 0; ii < nxyz; ii++) if (wf[ii] > mx) mx = wf[ii];
    clip = 0.05f * mx;
    for (ii = 0; ii < nxyz; ii++) if (wf[ii] < clip) wf[ii] = 0.0f;

    /* Normalize to 0..1 */
    mx = 0.0f;
    for (ii = 0; ii < nxyz; ii++) if (wf[ii] > mx) mx = wf[ii];
    if (mx > 0.0f) {
        float inv_mx = 1.0f / mx;
        for (ii = 0; ii < nxyz; ii++) wf[ii] *= inv_mx;
    }

    return wf;
}

/* --- Source-intensity percentile for the automask --------------------------
 * al_automask needs only the ~98th intensity percentile to set a background
 * threshold. The old code sorted all voxels via libc qsort just to read one
 * value. Default behavior below keeps that exact threshold while using O(n)
 * quickselect instead of O(n log n) qsort. The robust histogram estimator is
 * available only as an explicit behavior change via -DAL_AUTOMASK_ROBUST_RANGE. */

#ifndef AL_AUTOMASK_ROBUST_RANGE
/* Exact p-th percentile (p in 0..1) via 3-way quickselect on a scratch copy.
 * O(n) average; the 3-way (Dutch-flag) partition stays O(n) even when millions of
 * voxels share a value (image background) — the case that degrades 2-way schemes
 * and that a quicksort hits as worst-case. Returns the same value the old
 * sort-then-index did (sorted[(int)(p*nvox)]). */
static int al_pct98_quickselect(const float *src, int nvox, float p, float *out)
{
    if (nvox < 1) { *out = 0.0f; return 0; }
    float *a = (float *)malloc(sizeof(float) * nvox);
    if (!a) return 1;   /* OOM: signal failure so the caller can fail closed */
    /* Copy only finite values: al_select_rank's total-order partition is broken by
     * NaN (both `<` and `>` compare false), which would make the percentile — and
     * thus the mask threshold — arbitrary. NaN voxels never pass `> thresh` anyway. */
    int m = 0;
    for (int i = 0; i < nvox; i++) if (al_finitef(src[i])) a[m++] = src[i];
    if (m < 1) { free(a); *out = 0.0f; return 0; }   /* no finite data: empty mask */
    int k = (int)(p * m); if (k >= m) k = m - 1; if (k < 0) k = 0;
    *out = al_select_rank(a, m, k);   /* percentile over the finite values */
    free(a);
    return 0;
}
#endif

#ifdef AL_AUTOMASK_ROBUST_RANGE
/* FSL-style robust "98th percentile" (port of nifti_robust_range for a float
 * array): 1001-bin histogram with zero handling + the binary-image expansion
 * trick. O(n), no sort/comparator — ideal for (effectively 16-bit) scan data and
 * immune to the emscripten qsort pathology. This intentionally changes the mask
 * threshold on some inputs, so it is opt-in. */
static float al_pct98_robust(const float *src, int nvox)
{
    const int ignoreZero = 0;
    if (nvox < 1) return 1.0f;
    float mn = INFINITY, mx = -INFINITY;
    size_t nZero = 0, nSkip = 0;
    for (int i = 0; i < nvox; i++) {
        float v = src[i];
        if (!al_finitef(v)) { nSkip++; continue; }  /* NaN/inf */
        if (v == 0.0f) { nZero++; if (ignoreZero) continue; }
        if (v < mn) mn = v;
        if (v > mx) mx = v;
    }
    if ((nZero > 0) && (mn > 0.0f) && !ignoreZero) mn = 0.0f;
    if (mn > mx) return 1.0f;        /* everything skipped */
    if (mn == mx) return mx;
    size_t excl = ignoreZero ? (nZero + nSkip) : nSkip;
    size_t n2pct = (size_t)((double)(nvox - excl) * 0.02 + 0.5);
    if (n2pct < 1 || (size_t)nvox - excl < 100) return mx; /* tiny volume */
    enum { NB = 1001 };
    float scl = (NB - 1) / (mx - mn);
    int hist[NB]; for (int i = 0; i < NB; i++) hist[i] = 0;
    for (int i = 0; i < nvox; i++) {
        float v = src[i];
        if (!al_finitef(v)) continue;
        if (ignoreZero && v == 0.0f) continue;
        int b = (int)((v - mn) * scl + 0.5f);
        if (b < 0) b = 0; else if (b >= NB) b = NB - 1;
        hist[b]++;
    }
    size_t n = 0; int lo = 0;
    while (n < n2pct && lo < NB) { n += hist[lo]; lo++; }
    lo--;
    n = 0; int hi = NB;
    while (n < n2pct && hi > 0) { hi--; n += hist[hi]; }
    if (lo == hi) {  /* majority share one bin: expand to nearest non-empty bins */
        int ok = -1;
        while (ok != 0) {
            if (lo > 0) { lo--; if (hist[lo] > 0) ok = 0; }
            if (ok != 0 && hi < NB - 1) { hi++; if (hist[hi] > 0) ok = 0; }
            if (lo == 0 && hi == NB - 1) ok = 0;
        }
    }
    return (float)hi / scl + mn;     /* pct98 */
}
#endif

/* Compute simple source automask (binary): keep voxels above 10% of the ~98th
 * intensity percentile. See the percentile helpers above for the O(n) methods. */
static unsigned char *al_automask(float *srcim, int nvox)
{
    unsigned char *mask = (unsigned char *)calloc(nvox, sizeof(unsigned char));
    if (!mask) return NULL;
#ifdef AL_AUTOMASK_ROBUST_RANGE
    float pct98 = al_pct98_robust(srcim, nvox);
#else
    float pct98;
    if (al_pct98_quickselect(srcim, nvox, 0.98f, &pct98)) {
        /* Percentile scratch OOM: fail closed rather than returning a broad mask
         * built from a bogus zero threshold (an observable, not silent, failure). */
        free(mask);
        return NULL;
    }
#endif
    float thresh = 0.10f * pct98;
    for (int ii = 0; ii < nvox; ii++)
        mask[ii] = (srcim[ii] > thresh) ? 1 : 0;
    return mask;
}

/*==========================================================================*/
/*============= SECTION 10: BEFAFTER COORDINATE SETUP =====================*/
/*==========================================================================*/

static void al_setup_befafter(mat44 base_cmat, mat44 targ_imat)
{
    /* The affine parameters operate in xyz space:
       base index -> base xyz [before] -> transformed source xyz -> source index [after]. */
    aff_before = base_cmat;
    aff_after  = targ_imat;
    aff_use_before = 1;
    aff_use_after  = 1;
}

/* Clean up internal macros before the public API section */
#undef FAR
#undef CLIP
#undef ISTINY
#undef P_M1
#undef P_00
#undef P_P1
#undef P_P2
#undef P_FACTOR
#undef FAS
#undef SHANENT
#undef XYC
#undef WAY_BIG
#undef GOODVAL
#undef RANGVAL
#undef WV

/*==========================================================================*/
/*================ SECTION 11: PUBLIC FUNCTION nii_allineate() ============*/
/*==========================================================================*/

/* Map AL_COST_* to internal GA_MATCH_* code and name string.
   By default, -cost lpc uses pure local Pearson correlation (fast, matching AFNI).
   Compile with -DAL_LPC_MICHO to use the lpc+ZZ combined cost instead. */
static void al_resolve_cost(int cost, int *match_out, const char **name_out)
{
    switch (cost) {
#ifdef AL_LPC_MICHO
        case AL_COST_LPC:       *match_out = GA_MATCH_LPC_MICHO_SCALAR; *name_out = "lpc+ZZ"; return;
        case AL_COST_LPA:       *match_out = GA_MATCH_LPA_MICHO_SCALAR; *name_out = "lpa+ZZ"; return;
#else
        case AL_COST_LPC:       *match_out = GA_MATCH_PEARSON_LOCALS;   *name_out = "lpc"; return;
        case AL_COST_LPA:       *match_out = GA_MATCH_PEARSON_LOCALA;   *name_out = "lpa"; return;
#endif
        case AL_COST_PEARSON:   *match_out = GA_MATCH_PEARSON_SCALAR;   *name_out = "Pearson"; return;
        default:
        case AL_COST_HELLINGER: *match_out = GA_MATCH_HELLINGER_SCALAR; *name_out = "Hellinger"; return;
    }
}

/* Compute intensity-weighted center of mass in voxel indices */
static void al_center_of_mass(const float *data, int nx, int ny, int nz,
                               double *cx, double *cy, double *cz)
{
    double ws = 0.0, xc = 0.0, yc = 0.0, zc = 0.0;
    for (int kk = 0; kk < nz; kk++)
        for (int jj = 0; jj < ny; jj++)
            for (int ii = 0; ii < nx; ii++) {
                float fv = data[ii + jj * nx + kk * nx * ny];
                /* Skip non-finite voxels (NaN/±Inf reach the engine unvalidated and,
                 * under -ffast-math, `NaN <= 0` is false, so a single NaN would poison
                 * ws/xc and make nii_center_of_mass fold NaN into the sform/qform). The
                 * magnitude guard (al_finitef, not isfinite) is required by -ffast-math. */
                if (!al_finitef(fv) || fv <= 0.0f) continue;
                double v = (double)fv;
                ws += v; xc += v * ii; yc += v * jj; zc += v * kk;
            }
    if (ws > 0.0) { *cx = xc / ws; *cy = yc / ws; *cz = zc / ws; }
    else { *cx = (nx - 1) * 0.5; *cy = (ny - 1) * 0.5; *cz = (nz - 1) * 0.5; }
}



/* Overflow-safe 3D voxel count (nx*ny*nz), kept within int range so the int loop
 * counters/indices used throughout allineate stay valid. NIfTI dims are int64; this
 * uses division-based checked multiplication so the product is never formed unless it
 * is provably <= INT_MAX. Returns 0 and sets *out on success, 1 (message) on failure. */
static int al_safe_nvox(const nifti_image *n, const char *who, size_t *out)
{
    long long nx = n->nx, ny = n->ny, nz = n->nz;
    if (nx <= 0 || ny <= 0 || nz <= 0) {
        fprintf(stderr, "%s: non-positive dimensions\n", who); return 1; }
    if (nx > (long long)INT_MAX / ny || nx * ny > (long long)INT_MAX / nz) {
        fprintf(stderr, "%s: voxel count out of range\n", who); return 1; }
    *out = (size_t)(nx * ny * nz);   /* proven <= INT_MAX */
    return 0;
}

/* Extract float data from NIfTI image. Returns malloc'd float array (NULL on bad
 * dims/datatype) — the dimension guard is centralized here so every caller (allineate,
 * deface, reslice) is protected without each needing to pre-validate. */
static float *nii_to_float(nifti_image *nim)
{
    size_t ii;
    float *fdata;
    if (nim == NULL || nim->data == NULL) return NULL;
    size_t nvox;
    if (al_safe_nvox(nim, "nii_to_float", &nvox)) return NULL;
    fdata = (float *)calloc(nvox, sizeof(float));
    if (!fdata) return NULL;

    switch (nim->datatype) {
        case DT_UINT8: {
            unsigned char *p = (unsigned char *)nim->data;
            for (ii = 0; ii < nvox; ii++) fdata[ii] = (float)p[ii];
        } break;
        case DT_INT16: {
            short *p = (short *)nim->data;
            for (ii = 0; ii < nvox; ii++) fdata[ii] = (float)p[ii];
        } break;
        case DT_INT32: {
            int *p = (int *)nim->data;
            for (ii = 0; ii < nvox; ii++) fdata[ii] = (float)p[ii];
        } break;
        case DT_FLOAT32: {
            float *p = (float *)nim->data;
            for (ii = 0; ii < nvox; ii++) fdata[ii] = p[ii];
        } break;
        case DT_FLOAT64: {
            double *p = (double *)nim->data;
            for (ii = 0; ii < nvox; ii++) fdata[ii] = (float)p[ii];
        } break;
        case DT_INT8: {
            signed char *p = (signed char *)nim->data;
            for (ii = 0; ii < nvox; ii++) fdata[ii] = (float)p[ii];
        } break;
        case DT_UINT16: {
            unsigned short *p = (unsigned short *)nim->data;
            for (ii = 0; ii < nvox; ii++) fdata[ii] = (float)p[ii];
        } break;
        case DT_UINT32: {
            unsigned int *p = (unsigned int *)nim->data;
            for (ii = 0; ii < nvox; ii++) fdata[ii] = (float)p[ii];
        } break;
        case DT_INT64: {
            int64_t *p = (int64_t *)nim->data;
            for (ii = 0; ii < nvox; ii++) fdata[ii] = (float)p[ii];
        } break;
        case DT_UINT64: {
            uint64_t *p = (uint64_t *)nim->data;
            for (ii = 0; ii < nvox; ii++) fdata[ii] = (float)p[ii];
        } break;
        default:
            free(fdata);
            return NULL;
    }

    /* Apply scl_slope/scl_inter if set */
    if (nim->scl_slope != 0.0f &&
        !(nim->scl_slope == 1.0f && nim->scl_inter == 0.0f)) {
        for (ii = 0; ii < nvox; ii++)
            fdata[ii] = fdata[ii] * nim->scl_slope + nim->scl_inter;
    }

    return fdata;
}

/* Convert an image to *plain* float32 in place (datatype float32 with no pending
 * scl_slope/scl_inter), updating ALL datatype metadata (datatype/nbyper/swapsize and
 * clearing scale + cal_min/cal_max). Local mirror of miniCoreFLT's nii_ensure_float32
 * kept inside this TU so allineate.c stays a clean niimath drop-in (miniCoreFLT is
 * demo-only). The predicate matches aim_alias/nii_reslice_affine: a float32 image with
 * a PENDING scale is still converted, so callers that read raw samples (e.g. the deface
 * mask fill) see physical intensities and a negative slope can't invert min<->max.
 * Returns 0 on success, -1 (message via `who`) on an unsupported datatype. */
static int al_ensure_float32(nifti_image *nim, const char *who)
{
    if (nim->datatype == DT_FLOAT32 &&
        (nim->scl_slope == 0.0f || (nim->scl_slope == 1.0f && nim->scl_inter == 0.0f)))
        return 0;
    float *fdata = nii_to_float(nim);   /* applies any pending scl_slope/scl_inter */
    if (!fdata) { fprintf(stderr, "%s: unsupported input datatype %d\n", who, nim->datatype); return -1; }
    free(nim->data);
    nim->data      = fdata;
    nim->datatype  = DT_FLOAT32;
    nim->nbyper    = sizeof(float);
    nim->swapsize  = sizeof(float);
    nim->scl_slope = 0.0f;
    nim->scl_inter = 0.0f;
    nim->cal_min   = 0.0f;
    nim->cal_max   = 0.0f;
    return 0;
}

/* Normalized distance between two parameter vectors (in [0..1] space) */
static float param_dist(GA_setup *stup, float *p1, float *p2)
{
    float dmax = 0.0f;
    int jj;
    for (jj = 0; jj < stup->wfunc_numpar; jj++) {
        if (stup->wfunc_param[jj].fixed) continue;
        float s = stup->wfunc_param[jj].siz;
        if (s <= 0.0f) continue;
        float d = fabsf(p1[jj] - p2[jj]) / s;
        if (d > dmax) dmax = d;
    }
    return dmax;
}

/* Internal: run affine registration, return 12 warp parameters.
   source: moving image, base: fixed/reference image.
   match_code: GA_MATCH_* code for cost function.
   do_cmass: if nonzero, compute center-of-mass shift as initial alignment.
   do_src_automask: if nonzero, fill outside source automask with random noise.
   do_dark_automask: if nonzero, drop matched pairs where the base or warped-source
     value is at that image's darkest value (background/pad) from the cost.
   relax_scale: if nonzero, widen the affine scale-parameter range to [0.5,2.0]
     (from the default 0.711..1.406) — the caller (-zoom) is flagging abnormal size.
   wpar_out[12]: filled with optimized affine parameters on success.
   Returns 0 on success, nonzero on error. */
/* free_mask (optional): when non-NULL, a 12-entry 0/1 selector of which affine
   parameters are free (1) vs. permanently fixed (0), allowing a NON-contiguous
   free set (e.g. -sagseed's in-MSP {1,2,4}). When NULL, the common
   contiguous-prefix behavior applies: params[0 .. warp_dof-1] free, the rest
   fixed. */
static int al_register(nifti_image *source, nifti_image *base,
                       int match_code, int do_cmass, int do_src_automask,
                       int do_dark_automask, int relax_scale,
                       int fine_interp_code, int warp_dof, float wpar_out[12],
                       const int *free_mask)
{
    int reg_rc = 0;   /* nonzero -> abort via al_cleanup (setup/optimizer OOM) */
    GA_setup stup;
    GA_param params[12];
    float *bsim, *ajim;
    float *wght;
    unsigned char *smask;
    int bnx, bny, bnz, anx, any, anz;
    float bdx, bdy, bdz, adx, ady, adz;
    mat44 base_cmat, base_imat, targ_cmat, targ_imat;
    int jj, ii, nvox_base, nvox_src;
    int nfunc;

    if (source == NULL || base == NULL) {
        fprintf(stderr, "allineate: NULL input image\n");
        return 1;
    }
    /* Validate the DOF count for the contiguous-prefix path: warp_dof drives
       `for (jj = warp_dof; jj < 12; ...) params[jj]...`, so warp_dof < 0 would write
       out of bounds before params[0] and warp_dof == 0 leaves no free parameters.
       The CLI clamps this, but a direct API caller may not; fail closed. When
       free_mask is supplied it selects the free set and warp_dof is ignored. */
    if (free_mask == NULL && (warp_dof < 1 || warp_dof > 12)) {
        fprintf(stderr, "allineate: invalid warp DOF %d (must be 1..12)\n", warp_dof);
        return 1;
    }

    /* --- 1. Extract float data --- */
    bsim = nii_to_float(base);
    ajim = nii_to_float(source);
    if (!bsim || !ajim) {
        fprintf(stderr, "allineate: failed to extract float data\n");
        free(bsim); free(ajim);
        return 1;
    }

    bnx = base->nx;   bny = base->ny;   bnz = base->nz;
    anx = source->nx;  any = source->ny;  anz = source->nz;
    bdx = (float)fabs(base->dx);   if (bdx <= 0.0f) bdx = 1.0f;
    bdy = (float)fabs(base->dy);   if (bdy <= 0.0f) bdy = 1.0f;
    bdz = (float)fabs(base->dz);   if (bdz <= 0.0f) bdz = 1.0f;
    adx = (float)fabs(source->dx); if (adx <= 0.0f) adx = 1.0f;
    ady = (float)fabs(source->dy); if (ady <= 0.0f) ady = 1.0f;
    adz = (float)fabs(source->dz); if (adz <= 0.0f) adz = 1.0f;

    size_t nvox_base_size, nvox_src_size, tmp_size;
    if (al_mul_size((size_t)bnx, (size_t)bny, &tmp_size) ||
        al_mul_size(tmp_size, (size_t)bnz, &nvox_base_size) ||
        al_mul_size((size_t)anx, (size_t)any, &tmp_size) ||
        al_mul_size(tmp_size, (size_t)anz, &nvox_src_size) ||
        nvox_base_size > INT_MAX || nvox_src_size > INT_MAX) {
        fprintf(stderr, "allineate: image dimensions are too large\n");
        free(bsim); free(ajim);
        return 1;
    }
    nvox_base = (int)nvox_base_size;
    nvox_src  = (int)nvox_src_size;

    /* -dark_automask: the darkest (minimum) value of each image, used to drop
       background/pad matched pairs from the cost (see GA_scalar_fitter). Finite-aware
       (skip NaN/Inf); if an image is entirely non-finite the threshold stays FLT_MAX
       so every pair is masked and the survivor-floor fallback restores the unmasked
       weights each eval — i.e. the feature is effectively off. Using each image's own minimum
       (not a hardcoded 0) keeps it correct for signed data such as CT Hounsfield. */
    float dark_base_val = 0.0f, dark_targ_val = 0.0f;
    if (do_dark_automask) {
        float bmin = FLT_MAX, amin = FLT_MAX;
        for (ii = 0; ii < nvox_base; ii++)
            if (al_finitef(bsim[ii]) && bsim[ii] < bmin) bmin = bsim[ii];
        for (ii = 0; ii < nvox_src; ii++)
            if (al_finitef(ajim[ii]) && ajim[ii] < amin) amin = ajim[ii];
        dark_base_val = bmin; dark_targ_val = amin;
        fprintf(stderr, " + Dark automask: dropping matched pairs at/below base=%.4g, source=%.4g\n",
                bmin, amin);
    }

    /* --- 2. Get index->world matrices (sform preferred when its code >= the
       qform's, else qform; if neither form is coded/usable fall back to a
       pixdim-centered frame — the same no-form policy nii_symmetry and
       nii_center_of_mass use, so a valid NIfTI/ANALYZE with both codes 0 still
       registers rather than erroring out). --- */
    al_image_xform_or_pixdim(base, &base_cmat, NULL);
    base_imat = nifti_mat44_inverse(base_cmat);

    al_image_xform_or_pixdim(source, &targ_cmat, NULL);
    targ_imat = nifti_mat44_inverse(targ_cmat);

    /* --- 3. Compute autoweight --- */
    PROFILE_START(autoweight);
    fprintf(stderr, " + Computing autoweight from base image\n");
    wght = al_autoweight(bsim, bnx, bny, bnz, bdx, bdy, bdz);
    if (wght == NULL) {
        fprintf(stderr, "allineate: failed to allocate autoweight image\n");
        free(bsim); free(ajim);
        return 1;
    }

    /* --- 4. Compute source automask --- */
    smask = al_automask(ajim, nvox_src);
    if (smask == NULL) {
        fprintf(stderr, "allineate: failed to allocate source automask\n");
        free(bsim); free(ajim); free(wght);
        return 1;
    }
    if (smask) {
        int nm = 0;
        for (ii = 0; ii < nvox_src; ii++) if (smask[ii]) nm++;
        fprintf(stderr, " + Source automask: %d of %d voxels (%.0f%%)\n",
                nm, nvox_src, 100.0f * nm / nvox_src);
    }

    /* --- 5. Set up GA_setup structure --- */
    memset(&stup, 0, sizeof(GA_setup));

    stup.bsim = bsim;     stup.bsims = NULL;
    stup.bnx = bnx;       stup.bny = bny;       stup.bnz = bnz;
    stup.bdx = bdx;       stup.bdy = bdy;       stup.bdz = bdz;
    stup.ajim = ajim;     stup.ajims = NULL;
    stup.anx = anx;       stup.any = any;        stup.anz = anz;
    stup.adx = adx;       stup.ady = ady;        stup.adz = adz;
    stup.ajmask = smask;
    stup.do_dark_automask = do_dark_automask;
    stup.dark_base = dark_base_val; stup.dark_targ = dark_targ_val;
    stup.base_cmat = base_cmat; stup.base_imat = base_imat;
    stup.targ_cmat = targ_cmat; stup.targ_imat = targ_imat;
    stup.base_di = mat44_colnorm(base_cmat, 0);
    stup.base_dj = mat44_colnorm(base_cmat, 1);
    stup.base_dk = mat44_colnorm(base_cmat, 2);
    stup.targ_di = mat44_colnorm(targ_cmat, 0);
    stup.targ_dj = mat44_colnorm(targ_cmat, 1);
    stup.targ_dk = mat44_colnorm(targ_cmat, 2);

    /* Weight image */
    if (wght) {
        stup.bwght = wght;
        stup.bmask = (unsigned char *)calloc(nvox_base, sizeof(unsigned char));
        if (stup.bmask == NULL) {
            fprintf(stderr, "allineate: failed to allocate base mask\n");
            free(bsim); free(ajim); free(smask); free(wght);
            return 1;
        }
        stup.nmask = 0;
        /* AFNI's default weight for the histogram/box-mode costs (Hellinger, MI, NMI, CR) is
           `-autobox` (auto_weight=3): BINARIZE the weight so every in-mask voxel is weighted
           equally. Only ls/lpc/lpa keep the graded intensity weight (auto_weight=1). A graded
           weight down-weights the cortical periphery, which lets the fit shrink/expand the
           overlap unpenalized (the T2w->T1 shrink: peripheral voxels pushed out of FOV cost
           little); the binary weight anchors the full brain extent, so the fit contracts to
           fill (det 0.94->0.74, matching AFNI's 0.72). See AGENTS.md. */
        int box_weight = (match_code != GA_MATCH_PEARSON_SCALAR &&
                          match_code != GA_MATCH_PEARSON_LOCALS &&
                          match_code != GA_MATCH_PEARSON_LOCALA);
        float wmx = 0.0f;
        for (ii = 0; ii < nvox_base; ii++) if (wght[ii] > wmx) wmx = wght[ii];
        if (wmx > 0.0f) {
            float inv = 1.0f / wmx;
            for (ii = 0; ii < nvox_base; ii++) {
                wght[ii] = fabsf(wght[ii]) * inv;
                stup.bmask[ii] = (wght[ii] > 0.0f) ? 1 : 0;
                if (stup.bmask[ii]) { stup.nmask++; if (box_weight) wght[ii] = 1.0f; }
            }
        }
    }

    /* --- 5b. Source automask noise-fill setup (AFNI's -source_automask) --- */
    stup.ajmask_ranfill = 0;
    stup.aj_ubot = stup.aj_usiz = 0.0f;
    stup.ajim_orig = NULL;
    if (do_src_automask && smask) {
        /* Compute noise range from 7th and 93rd percentiles of masked voxels */
        int nm = 0;
        for (ii = 0; ii < nvox_src; ii++) if (smask[ii]) nm++;
        if (nm > 10) {
            float *mvals = (float *)malloc(sizeof(float) * nm);
            if (mvals == NULL) {
                fprintf(stderr, "allineate: failed to allocate source automask scratch\n");
                free(bsim); free(ajim); free(smask); free(stup.bmask); free(wght);
                return 1;
            }
            int pp = 0;
            for (ii = 0; ii < nvox_src; ii++)
                if (smask[ii]) mvals[pp++] = ajim[ii];
            /* O(n) selects (was a full qsort+comparator over up to ~millions of
             * masked voxels); two passes for the two percentiles. */
            float q07 = al_select_rank(mvals, nm, (int)(0.07f * (nm - 1)));
            float q93 = al_select_rank(mvals, nm, (int)(0.93f * (nm - 1)));
            free(mvals);
            stup.aj_ubot = q07;
            stup.aj_usiz = (q93 - q07) * 0.5f;
            if (stup.aj_usiz > 0.0f) {
                stup.ajmask_ranfill = 1;
                /* Backup original source data before noise fill */
                stup.ajim_orig = (float *)malloc(sizeof(float) * nvox_src);
                if (stup.ajim_orig == NULL) {
                    fprintf(stderr, "allineate: failed to allocate source automask backup\n");
                    free(bsim); free(ajim); free(smask); free(stup.bmask); free(wght);
                    return 1;
                }
                memcpy(stup.ajim_orig, ajim, sizeof(float) * nvox_src);
                fprintf(stderr, " + Source automask noise fill: range [%.1f, %.1f]\n",
                        q07, q07 + 2.0f * stup.aj_usiz);
            }
        }
    }

    /* Warn if using lpc/lpa without source_automask */
    if (!do_src_automask && METH_USES_BLOKS(match_code)) {
        fprintf(stderr, " + WARNING: '-source_automask' is strongly recommended"
                        " when using -lpc or -lpa\n");
    }

    /* --- 6. Configure cost function --- */
    stup.match_code = match_code;
#ifdef AL_LPC_MICHO
    if (match_code == GA_MATCH_PEARSON_SCALAR  ||
        match_code == GA_MATCH_HELLINGER_SCALAR ||
        match_code == GA_MATCH_PEARSON_LOCALS   ||
        match_code == GA_MATCH_PEARSON_LOCALA) {
        stup.micho_mi = stup.micho_nmi = stup.micho_crA = stup.micho_hel = stup.micho_ov = 0.0;
        stup.micho_zfinal = 0;
    } else if (METH_IS_LPA(match_code)) {
        /* LPA-specific defaults (AFNI 27 May 2021) */
        stup.micho_mi  = 0.0;
        stup.micho_nmi = 0.2;
        stup.micho_crA = 0.4;
        stup.micho_hel = 0.4;
        stup.micho_ov  = 0.4;
        stup.micho_zfinal = 1; /* +ZZ */
    } else {
        /* LPC defaults */
        stup.micho_mi  = 0.2;
        stup.micho_nmi = 0.2;
        stup.micho_crA = 0.4;
        stup.micho_hel = 0.4;
        stup.micho_ov  = 0.4;
        stup.micho_zfinal = 1; /* +ZZ */
    }
#endif

    /* --- 7. Set up affine 12-param warp --- */
    stup.wfunc = al_wfunc_affine;
    stup.wfunc_numpar = 12;
    stup.wfunc_param = params;
    memset(params, 0, sizeof(params));

    /* Set up befafter matrices */
    al_setup_befafter(base_cmat, targ_imat);

    /* AFNI derives shift limits from a zero-padded fixed FOV. Expanding every
       case perturbs this compact port's normalized random search, so preserve the
       established range unless the brightness-centroid displacement approaches
       an unpadded boundary. This changes the allowed search, not the starting pose:
       -cmass/-nocmass retain their explicit semantics. */
    float xxx_m, yyy_m, zzz_m;
    al_shift_range_mm(base_cmat, bnx, bny, bnz, &xxx_m, &yyy_m, &zzz_m);
    int pxm, pxp, pym, pyp, pzm, pzp;
    al_fixed_pad_widths(bsim, bnx, bny, bnz, &pxm, &pxp, &pym, &pyp, &pzm, &pzp);
    size_t pnxz = (size_t)bnx + pxm + pxp;
    size_t pnyz = (size_t)bny + pym + pyp;
    size_t pnzz = (size_t)bnz + pzm + pzp;
    if (pnxz <= INT_MAX && pnyz <= INT_MAX && pnzz <= INT_MAX &&
        (pnxz != (size_t)bnx || pnyz != (size_t)bny || pnzz != (size_t)bnz)) {
        int pnx = (int)pnxz, pny = (int)pnyz, pnz = (int)pnzz;
        double bxc, byc, bzc, axc, ayc, azc;
        float bx, by, bz, ax, ay, az;
        al_center_of_mass(bsim, bnx, bny, bnz, &bxc, &byc, &bzc);
        al_center_of_mass(ajim, anx, any, anz, &axc, &ayc, &azc);
        mat44_vec(base_cmat, (float)bxc, (float)byc, (float)bzc, &bx, &by, &bz);
        mat44_vec(targ_cmat, (float)axc, (float)ayc, (float)azc, &ax, &ay, &az);
        float cmx = ax - bx, cmy = ay - by, cmz = az - bz;
        if (fabsf(cmx) > 0.9f * xxx_m || fabsf(cmy) > 0.9f * yyy_m ||
            fabsf(cmz) > 0.9f * zzz_m) {
            al_shift_range_mm(base_cmat, pnx, pny, pnz, &xxx_m, &yyy_m, &zzz_m);
            if (getenv("AL_VERB"))
                fprintf(stderr, "[AL_VERB range] centroid=(%+.1f,%+.1f,%+.1f) mm; AFNI padded dims %dx%dx%d -> %dx%dx%d\n",
                        cmx, cmy, cmz, bnx, bny, bnz, pnx, pny, pnz);
        }
    }

    /* DEFPAR macro equivalent */
#define SETPAR(p,nm,bb,tt,id) do {          \
    params[p].min = (bb); params[p].max = (tt);    \
    params[p].ident = (id); params[p].val_init = (id); \
    params[p].val_pinit = (id); params[p].val_fixed = (id); \
    params[p].val_out = (id); params[p].fixed = 0;   \
    strncpy(params[p].name, (nm), 31);           \
} while(0)

    /* Optional hard clamp on the shift range (nii_symmetry mirror registration). */
    if (al_shift_max_override > 0.0f) {
        if (xxx_m > al_shift_max_override) xxx_m = al_shift_max_override;
        if (yyy_m > al_shift_max_override) yyy_m = al_shift_max_override;
        if (zzz_m > al_shift_max_override) zzz_m = al_shift_max_override;
    }

    SETPAR(0, "x-shift", -xxx_m, xxx_m, 0.0f);
    SETPAR(1, "y-shift", -yyy_m, yyy_m, 0.0f);
    SETPAR(2, "z-shift", -zzz_m, zzz_m, 0.0f);

    /* --- Center-of-mass initial shift --- */
    if (do_cmass) {
        double bxc, byc, bzc, axc, ayc, azc;
        al_center_of_mass(bsim, bnx, bny, bnz, &bxc, &byc, &bzc);
        al_center_of_mass(ajim, anx, any, anz, &axc, &ayc, &azc);

        float bx_mm, by_mm, bz_mm, ax_mm, ay_mm, az_mm;
        mat44_vec(base_cmat, (float)bxc, (float)byc, (float)bzc, &bx_mm, &by_mm, &bz_mm);
        mat44_vec(targ_cmat, (float)axc, (float)ayc, (float)azc, &ax_mm, &ay_mm, &az_mm);

        float cmx = ax_mm - bx_mm, cmy = ay_mm - by_mm, cmz = az_mm - bz_mm;
        fprintf(stderr, " + Center-of-mass shift: %.2f %.2f %.2f mm\n", cmx, cmy, cmz);

        /* Clamp to parameter range and set as initial value */
        if (cmx < params[0].min) cmx = params[0].min;
        if (cmx > params[0].max) cmx = params[0].max;
        if (cmy < params[1].min) cmy = params[1].min;
        if (cmy > params[1].max) cmy = params[1].max;
        if (cmz < params[2].min) cmz = params[2].min;
        if (cmz > params[2].max) cmz = params[2].max;
        params[0].val_init = params[0].val_pinit = cmx;
        params[1].val_init = params[1].val_pinit = cmy;
        params[2].val_init = params[2].val_pinit = cmz;
    }

    float rval = 30.0f;
    SETPAR(3, "z-angle", -rval, rval, 0.0f);
    SETPAR(4, "x-angle", -rval, rval, 0.0f);
    SETPAR(5, "y-angle", -rval, rval, 0.0f);

    float smin = 0.711f, smax = 1.0f / smin;
    SETPAR(6, "x-scale", smin, smax, 1.0f);
    SETPAR(7, "y-scale", smin, smax, 1.0f);
    SETPAR(8, "z-scale", smin, smax, 1.0f);

    rval = 0.1111f;
    SETPAR(9,  "y/x-shear", -rval, rval, 0.0f);
    SETPAR(10, "z/x-shear", -rval, rval, 0.0f);
    SETPAR(11, "z/y-shear", -rval, rval, 0.0f);
#undef SETPAR

    /* Widen the scale-parameter range when the caller opts into abnormal sizes (-zoom).
       The -sagseed isotropic seed (al_zoom_isotropic) frees only param[6] — 7/8 follow
       it in GA_setup_affine; the main fit (relax_scale) frees all three, so widen 6..8.
       Neither set → the tight default 0.711..1.406 regularization is kept, so ordinary
       (adult) fits are unaffected. */
    if (al_zoom_isotropic || relax_scale) {
        const float zmin = 0.5f, zmax = 2.0f;   /* single source for both copies + the log */
        params[6].min = zmin; params[6].max = zmax;
        if (relax_scale) {   /* main fit frees all three; the isotropic seed only param[6] */
            params[7].min = params[8].min = zmin;
            params[7].max = params[8].max = zmax;
            fprintf(stderr, " + Zoom: relaxed affine scale range to [%.2f, %.2f] (size flagged abnormal)\n",
                    zmin, zmax);
        }
    }

    /* --- 8. BLOK set params (AFNI Jul 2021 defaults: TOHD, ~555 voxels/blok) --- */
    stup.bloktype = GA_BLOK_TOHD;
    { float vvv = stup.base_di * stup.base_dj * stup.base_dk;
      /* TOHD volume factor = 4.0; auto-radius for ~555 voxels per blok */
      stup.blokrad = cbrtf(555.0f * vvv / 4.0f);
    }
    stup.blokmin  = 0; /* auto */
    stup.blokset  = NULL;

    /* ========== COARSE PASS (twopass) ========== */
    PROFILE_END(autoweight, "autoweight + setup");
    PROFILE_START(coarse);
    double _coarse_t0 = al_wtime();   /* always-on coarse timer (see g_last_coarse_ms) */
    fprintf(stderr, " + *** Coarse pass begins ***\n");

    /* LPA gets more twobest candidates (AFNI 27 May 2021) */
    int tbest = 5;
    {
        double vol_src = (double)anx * adx * any * ady * anz * adz;
        double vol_base = (double)bnx * bdx * bny * bdy * bnz * bdz;
        if (vol_src > 1.3 * vol_base) tbest = PARAM_MAXTRIAL;
    }
    if (METH_IS_LPA(match_code) && tbest < DEFAULT_TBEST_LPA)
        tbest = DEFAULT_TBEST_LPA;

    /* Compute dxyz_top for BLOK smoothing floor */
    float dxyz_top = bdx;
    if (bdy > dxyz_top) dxyz_top = bdy;
    if (bdz > dxyz_top) dxyz_top = bdz;
    if (adx > dxyz_top) dxyz_top = adx;
    if (ady > dxyz_top) dxyz_top = ady;
    if (adz > dxyz_top) dxyz_top = adz;

    /* BLOK methods need less coarse smoothing to preserve local structure */
    float sm_rad = 0.0f;
    if (METH_USES_BLOKS(match_code)) {
        sm_rad = 2.222f;
        if (dxyz_top > sm_rad) sm_rad = dxyz_top;
    }

    /* Scale matching points based on source and mask sizes.
       Default (matching AFNI): geometric mean avoids redundant oversampling
       when source is much smaller than base (e.g., 2.4mm fMRI → 1mm T1).
       Compile with -DAL_MATCH_ALL to use all mask points (slower, marginally
       more precise for same-resolution images). */
    int ntask = nvox_base;
#ifdef AL_MATCH_ALL
    /* Use all mask voxels as matching points */
    int npt_match_full = ntask;
#else
    if (stup.nmask > 0) {
        ntask = (nvox_src < stup.nmask)
              ? (int)sqrt((double)nvox_src * stup.nmask) : stup.nmask;
    }
    int npt_match_full = (int)(AL_SPARSE_SAMPLE_FRAC * ntask);
    if (npt_match_full < AL_NPT_MATCH_MIN) npt_match_full = AL_NPT_MATCH_MIN;
    if (npt_match_full > ntask)             npt_match_full = ntask;
#endif

    /* Note: ntask may have been reduced by the geometric-mean scaling above,
       so the coarse-pass count (ntask/15) is also proportionally smaller. */
    stup.interp_code = AL_INTERP_LINEAR;
    stup.npt_match   = ntask / 15;
    if (stup.npt_match < AL_NPT_MATCH_MIN) stup.npt_match = AL_NPT_MATCH_MIN;
    stup.smooth_code = GA_SMOOTH_GAUSSIAN;
    stup.smooth_radius_base = stup.smooth_radius_targ =
        (sm_rad > 0.0f) ? sm_rad : 7.777f;

    al_scalar_setup(&stup);
    if (stup.setup != AL_SMAGIC) {  /* setup OOM (im_ar/bvm/wvm): fail closed */
        fprintf(stderr, "allineate: registration setup failed (out of memory)\n");
        reg_rc = 1; goto al_cleanup;
    }

    /* Permanently fix params outside the free set (fixed=2 survives unfreeze).
       free_mask (when non-NULL) selects an arbitrary free subset for constrained
       fits (e.g. -sagseed's non-contiguous {1,2,4}); NULL keeps the common
       contiguous-prefix warp_dof behavior. nparam_free is the count of free
       params either way, so the coarse-pass freeze below is unchanged. */
    int nparam_free;
    if (free_mask) {
        nparam_free = 0;
        for (jj = 0; jj < 12; jj++) {
            params[jj].fixed = free_mask[jj] ? 0 : 2;
            if (free_mask[jj]) nparam_free++;
        }
    } else {
        nparam_free = warp_dof;
        for (jj = nparam_free; jj < 12; jj++)
            params[jj].fixed = 2;
    }
    /* Fail closed on an empty free set: an all-zero free_mask would hand NEWUOA a
       zero-dimensional problem. Unreachable via the current callers (-sagseed frees
       {1,2,4}; the warp_dof path is validated to [1,12]) but guards the free_mask API. */
    if (nparam_free < 1) {
        fprintf(stderr, "allineate: no free parameters to optimize\n");
        reg_rc = 1; goto al_cleanup;
    }

    /* Temporarily freeze params beyond first 6 for coarse search */
    int nptwo = (nparam_free < 6) ? nparam_free : 6;
    if (nparam_free > nptwo) {
        for (ii = jj = 0; jj < stup.wfunc_numpar; jj++) {
            if (!stup.wfunc_param[jj].fixed) {
                ii++;
                if (ii > nptwo) stup.wfunc_param[jj].fixed = 1;
            }
        }
    }

    powell_set_mfac(1.0f, 3.0f);

    int nrand = 17 + 4 * tbest;
    if (METH_USES_BLOKS(match_code)) nrand += 2 * tbest;
    if (nrand < 31) nrand = 31;
    if (al_scalar_ransetup(&stup, nrand) != 0) {
        reg_rc = 1; goto al_cleanup;
    }
    if (getenv("AL_VERB"))
        fprintf(stderr, "[AL_VERB coarse] best rigid pose: shift=(%.1f,%.1f,%.1f) angle=(%.1f,%.1f,%.1f)\n",
                stup.wfunc_param[0].val_init, stup.wfunc_param[1].val_init, stup.wfunc_param[2].val_init,
                stup.wfunc_param[3].val_init, stup.wfunc_param[4].val_init, stup.wfunc_param[5].val_init);

    /* Unfreeze temporarily frozen params */
    for (jj = 0; jj < stup.wfunc_numpar; jj++)
        if (stup.wfunc_param[jj].fixed == 1) stup.wfunc_param[jj].fixed = 0;

    /* Store trial parameter sets */
    int tb = tbest;
    if (tb > stup.wfunc_ntrial) tb = stup.wfunc_ntrial;
    float tfparm[PARAM_MAXTRIAL + 2][12];
    float ffparm[PARAM_MAXTRIAL + 2][12];  /* scratch for sorting */
    float tfcost[PARAM_MAXTRIAL + 2];
    int tfindx[PARAM_MAXTRIAL + 2];
    int tfdone;

    for (int ib = 0; ib < tb; ib++)
        for (jj = 0; jj < 12; jj++)
            tfparm[ib][jj] = stup.wfunc_param[jj].val_trial[ib];

    /* Add identity transform */
    for (jj = 0; jj < 12; jj++)
        tfparm[tb][jj] = stup.wfunc_param[jj].val_pinit;
    tfdone = tb + 1;

    /* Refinement rounds — parallelize candidates within each round */
    double rad = 0.0555;
    int nrefine = 3;
    for (int rr = 0; rr < nrefine; rr++, rad *= 0.6789) {
        fprintf(stderr, " + Refinement #%d on %d parameter sets\n", rr + 1, tfdone);

        stup.smooth_radius_base *= 0.7777f;
        stup.smooth_radius_targ *= 0.7777f;
        stup.npt_match = (int)(stup.npt_match * 1.5);

        al_scalar_setup(&stup);
        if (stup.setup != AL_SMAGIC) {  /* setup OOM: fail closed */
            fprintf(stderr, "allineate: refinement setup failed (out of memory)\n");
            reg_rc = 1; goto al_cleanup;
        }
        GA_param_setup(&stup);
        gstup = &stup;

        int nfr_ref = stup.wfunc_numfree;
        int maxstep_ref = 99 + 11 * rr;
        if (maxstep_ref <= 4 * nfr_ref + 5) maxstep_ref = 6666;
        double rstart_ref = (rad > 0.2) ? 0.2 : rad;
        double rend_ref = 0.01 * rstart_ref;
        if (rend_ref >= 0.9 * rstart_ref) rend_ref = 0.0666 * rstart_ref;
        float mfac_m = 1.0f, mfac_a = 5.0f + 2.0f * rr;

        /* Convert candidates to normalized 0-1 parameter arrays */
        double *cand_wpar[PARAM_MAXTRIAL + 2] = {0};
        for (int ib = 0; ib < tfdone; ib++) {
            cand_wpar[ib] = (double *)calloc(nfr_ref, sizeof(double));
            if (!cand_wpar[ib]) {
                for (int jb = 0; jb < ib; jb++) free(cand_wpar[jb]);
                fprintf(stderr, "allineate: refinement candidate setup failed (out of memory)\n");
                reg_rc = 1; goto al_cleanup;
            }
            int qi = 0;
            for (jj = 0; jj < stup.wfunc_numpar; jj++) {
                if (!stup.wfunc_param[jj].fixed) {
                    double v = (tfparm[ib][jj] - stup.wfunc_param[jj].min)
                              / stup.wfunc_param[jj].siz;
                    cand_wpar[ib][qi++] = PRED01(v);
                }
            }
        }

        /* Consume need_hist_setup (armed by al_scalar_setup above) with ONE sequential eval at
           candidate 0 — refreshing the sample-based edge-bin clips before the parallel region,
           so the shared write is race-free (AFNI refreshes on candidate 0, whose refinement it
           runs first). */
        if (tfdone > 0 && stup.need_hist_setup) (void)GA_scalar_fitter(nfr_ref, cand_wpar[0]);

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic)
#endif
        for (int ib = 0; ib < tfdone; ib++) {
            powell_set_mfac(mfac_m, mfac_a);
            int prc = powell_newuoa(nfr_ref, cand_wpar[ib], rstart_ref, rend_ref,
                         maxstep_ref, GA_scalar_fitter);
            tfcost[ib] = (prc < 0) ? AL_BIGVAL : (float)GA_scalar_fitter(nfr_ref, cand_wpar[ib]);
        }

        /* Unpack results back to tfparm */
        for (int ib = 0; ib < tfdone; ib++) {
            int qi = 0;
            for (jj = 0; jj < stup.wfunc_numpar; jj++) {
                if (stup.wfunc_param[jj].fixed) {
                    tfparm[ib][jj] = stup.wfunc_param[jj].val_fixed;
                } else {
                    tfparm[ib][jj] = stup.wfunc_param[jj].min
                                    + stup.wfunc_param[jj].siz * (float)PRED01(cand_wpar[ib][qi]);
                    qi++;
                }
            }
            free(cand_wpar[ib]);
        }

        /* Sort by cost and deduplicate close parameter sets */
        if (tfdone > 2) {
            /* Copy into ffparm for sorting */
            for (int ib = 0; ib < tfdone; ib++) {
                memcpy(ffparm[ib], tfparm[ib], sizeof(float) * 12);
                tfindx[ib] = ib;
            }
            qsort_floatint(tfdone, tfcost, tfindx);
            for (int ib = 0; ib < tfdone; ib++)
                memcpy(tfparm[ib], ffparm[tfindx[ib]], sizeof(float) * 12);
            if (getenv("AL_VERB"))
                fprintf(stderr, "[AL_VERB refine#%d] best cost=%.5f  scale=(%.3f,%.3f,%.3f) "
                        "shear=(%.3f,%.3f,%.3f) shift=(%.1f,%.1f,%.1f) angle=(%.1f,%.1f,%.1f)\n", rr + 1, tfcost[0],
                        tfparm[0][6], tfparm[0][7], tfparm[0][8], tfparm[0][9], tfparm[0][10], tfparm[0][11],
                        tfparm[0][0], tfparm[0][1], tfparm[0][2], tfparm[0][3], tfparm[0][4], tfparm[0][5]);

            /* Cast out parameter sets too close to the best */
#define CTHRESH 0.02f
            for (int ib = 1; ib < tfdone; ib++) {
                float pdist = param_dist(&stup, tfparm[0], tfparm[ib]);
                if (tfdone > 2 && pdist < CTHRESH) {
                    for (int jb = ib + 1; jb < tfdone; jb++) {
                        memcpy(tfparm[jb - 1], tfparm[jb], sizeof(float) * 12);
                        tfcost[jb - 1] = tfcost[jb];
                    }
                    tfdone--;
                    ib--;  /* re-check the entry shifted into this position */
                }
            }
#undef CTHRESH
        }
    }

    /* ========== FINE PASS ========== */
    PROFILE_END(coarse, "coarse pass");
    g_last_coarse_ms = (al_wtime() - _coarse_t0) * 1000.0;
    PROFILE_START(fine);
    double _fine_t0 = al_wtime();   /* always-on fine timer (see g_last_fine_ms) */
    fprintf(stderr, " + *** Fine pass begins ***\n");

    /* Fine pass: full resolution, no smoothing.
       Default is linear (matching AFNI); cubic is slower but marginally more precise. */
    stup.interp_code = fine_interp_code;
    stup.smooth_code = 0;
    stup.smooth_radius_base = 0.0f;
    stup.smooth_radius_targ = 0.0f;
    stup.npt_match = npt_match_full;

    al_scalar_setup(&stup);
    if (stup.setup != AL_SMAGIC) {  /* setup OOM: fail closed */
        fprintf(stderr, "allineate: final setup failed (out of memory)\n");
        reg_rc = 1; goto al_cleanup;
    }
    powell_set_mfac(0.0f, 0.0f);

    /* Refine all candidates at full resolution then pick the best */
    PROFILE_START(fine_cand);
    {
        int kb = 0;
        float cbest = 1.e+33f;
        int num_rtb = 99;
        float cand_cost[PARAM_MAXTRIAL + 2];

        rad = (tfdone > 2) ? 0.0333 : 0.0444;

        /* Prepare per-candidate parameter arrays (normalized 0-1 for powell) */
        int nfr = stup.wfunc_numfree;
        double *cand_wpar[PARAM_MAXTRIAL + 2] = {0};
        int cand_rtb[PARAM_MAXTRIAL + 2];
        for (int ib = 0; ib < tfdone; ib++) {
            cand_wpar[ib] = (double *)calloc(nfr, sizeof(double));
            if (!cand_wpar[ib]) {
                for (int jb = 0; jb < ib; jb++) free(cand_wpar[jb]);
                fprintf(stderr, "allineate: fine candidate setup failed (out of memory)\n");
                reg_rc = 1; goto al_cleanup;
            }
            cand_rtb[ib] = (ib == tfdone - 1) ? 2 * num_rtb : num_rtb;
            /* Convert from val_init space to normalized 0-1 for powell */
            int qi = 0;
            for (jj = 0; jj < stup.wfunc_numpar; jj++) {
                if (!stup.wfunc_param[jj].fixed) {
                    double v = (tfparm[ib][jj] - stup.wfunc_param[jj].min)
                              / stup.wfunc_param[jj].siz;
                    cand_wpar[ib][qi++] = PRED01(v);
                }
            }
        }

        gstup = &stup;
        /* Consume need_hist_setup (armed by the fine-pass al_scalar_setup) with ONE sequential
           eval at candidate 0 — refresh the sample-based clips before the parallel region so the
           shared write stays race-free and p1==pN bit-identical. */
        if (tfdone > 0 && stup.need_hist_setup) (void)GA_scalar_fitter(nfr, cand_wpar[0]);
        /* re-apply the main thread's thread-local sampling factors per worker */
        float fc_mfac, fc_afac; powell_get_mfac(&fc_mfac, &fc_afac);
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic)
#endif
        for (int ib = 0; ib < tfdone; ib++) {
            int maxstep = cand_rtb[ib];
            if (maxstep <= 4 * nfr + 5) maxstep = 6666;
            powell_set_mfac(fc_mfac, fc_afac);
            int prc = powell_newuoa(nfr, cand_wpar[ib], rad, 0.01 * rad, maxstep, GA_scalar_fitter);
            cand_cost[ib] = (prc < 0) ? AL_BIGVAL : (float)GA_scalar_fitter(nfr, cand_wpar[ib]);
        }

        /* Unpack results back to ffparm */
        for (int ib = 0; ib < tfdone; ib++) {
            int qi = 0;
            for (jj = 0; jj < stup.wfunc_numpar; jj++) {
                if (stup.wfunc_param[jj].fixed) {
                    ffparm[ib][jj] = stup.wfunc_param[jj].val_fixed;
                } else {
                    ffparm[ib][jj] = stup.wfunc_param[jj].min
                                    + stup.wfunc_param[jj].siz * PRED01(cand_wpar[ib][qi]);
                    qi++;
                }
            }
            free(cand_wpar[ib]);
            if (cand_cost[ib] < cbest) { cbest = cand_cost[ib]; kb = ib; }
        }
        for (jj = 0; jj < stup.wfunc_numpar; jj++)
            stup.wfunc_param[jj].val_init = ffparm[kb][jj];
        fprintf(stderr, " + Best fine cost = %f (set #%d of %d)\n", cbest, kb + 1, tfdone);
    }
    FITTER_PROFILE_DUMP("fine candidates fitter");
    PROFILE_END(fine_cand, "fine candidate refinement");

    /* Final optimization: full resolution, user-selected interpolation */
    PROFILE_START(fine_final);
    rad = 0.0333;
    nfunc = al_scalar_optim(&stup, rad, 0.001, 6666);
    if (nfunc < 0) {   /* OOM in the optimizer: fail closed, do not emit an unrefined result */
        fprintf(stderr, "allineate: final optimization failed (out of memory)\n");
        reg_rc = 1; goto al_cleanup;
    }
    fprintf(stderr, " + Fine cost = %f (%d evaluations)\n",
            stup.vbest, nfunc);
    FITTER_PROFILE_DUMP("fine final fitter");
    PROFILE_END(fine_final, "fine final optimization");

#ifdef AL_LPC_MICHO
    /* --- +ZZ refinal: zero micho weights, re-optimize with pure local correlation --- */
    if (METH_USES_BLOKS(match_code) && stup.micho_zfinal) {
        for (jj = 0; jj < stup.wfunc_numpar; jj++)
            stup.wfunc_param[jj].val_init = stup.wfunc_param[jj].val_out;

        double save_mi  = stup.micho_mi;
        double save_nmi = stup.micho_nmi;
        double save_crA = stup.micho_crA;
        double save_hel = stup.micho_hel;
        double save_ov  = stup.micho_ov;
        stup.micho_mi = stup.micho_nmi = stup.micho_crA = stup.micho_hel = stup.micho_ov = 0.0;
        fprintf(stderr, " + +ZZ refinal (pure %s)\n",
                METH_IS_LPA(match_code) ? "lpa" : "lpc");
        rad = 0.0666;
        powell_set_mfac(4.0f, 4.0f);
        nfunc = al_scalar_optim(&stup, rad, 0.001, 6666);
        if (nfunc < 0) {
            fprintf(stderr, "allineate: +ZZ refinal failed (out of memory)\n");
            reg_rc = 1; goto al_cleanup;
        }
        fprintf(stderr, " + Refinal cost = %f (%d evaluations)\n", stup.vbest, nfunc);

        stup.micho_mi  = save_mi;
        stup.micho_nmi = save_nmi;
        stup.micho_crA = save_crA;
        stup.micho_hel = save_hel;
        stup.micho_ov  = save_ov;
    }
#endif

    /* --- Extract final parameters --- */
    for (jj = 0; jj < 12; jj++)
        wpar_out[jj] = stup.wfunc_param[jj].val_out;

    fprintf(stderr, " + Final parameters:");
    for (jj = 0; jj < 12; jj++) fprintf(stderr, " %.4f", wpar_out[jj]);
    fprintf(stderr, "\n");

    /* Fine-pass timing lives here (success path only): the coarse-pass OOM gotos
       jump over PROFILE_START(fine) straight to al_cleanup, so ending the timer in
       the shared cleanup would read an uninitialized _prof_fine. */
    PROFILE_END(fine, "fine pass total");
    g_last_fine_ms = (al_wtime() - _fine_t0) * 1000.0;

    /* --- Cleanup --- */
al_cleanup:
    /* Restore original source data if noise-filled, before freeing */
    if (stup.ajim_orig) {
        memcpy(ajim, stup.ajim_orig, sizeof(float) * nvox_src);
        free(stup.ajim_orig);
    }
    free(bsim);
    free(ajim);
    if (stup.bmask) free(stup.bmask);
    if (stup.bsims) free(stup.bsims);
    if (stup.ajims) free(stup.ajims);
    if (stup.bvm) free(stup.bvm);
    if (stup.wvm) free(stup.wvm);
    if (stup.im_ar) free(stup.im_ar);
    if (stup.jm_ar) free(stup.jm_ar);
    if (stup.km_ar) free(stup.km_ar);
    if (stup.blokset) free_GA_BLOK_set(stup.blokset);
    if (smask) free(smask);
    if (stup.bwght) free(stup.bwght);
    /* Free thread-local buffers for all threads (not just the main thread) */
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        clear_2Dhist();  /* each thread has its own AL_TLOCAL histogram */
        free(tl_avm);  tl_avm = NULL;  tl_avm_len = 0;
        free(tl_weff); tl_weff = NULL; tl_weff_len = 0;
        free(tl_wpar); tl_wpar = NULL; tl_wpar_len = 0;
        free(tl_wbuf); tl_wbuf = NULL; tl_wbuf_len = 0;
        powell_newuoa_free_threadlocal();  /* free this thread's NEWUOA workspace */
    }

    if (reg_rc == 0)
        fprintf(stderr, " + Registration complete\n");
    return reg_rc;
}

/* Adopt base grid geometry into source (dims + sform/qform); source->data must
 * already hold the float volume sampled on the base grid. */
static void al_adopt_geometry(nifti_image *s, const nifti_image *b)
{
    s->datatype = DT_FLOAT32;
    s->nbyper = sizeof(float);
    s->swapsize = sizeof(float);   /* keep element swap width consistent for write */
    s->scl_slope = 0.0f;
    s->scl_inter = 0.0f;
    s->cal_min = 0.0f;
    s->cal_max = 0.0f;
    s->ndim = 3;
    s->nx = b->nx; s->ny = b->ny; s->nz = b->nz;
    s->nt = s->nu = s->nv = s->nw = 1;
    s->dx = b->dx; s->dy = b->dy; s->dz = b->dz;
    s->dt = b->dt; s->du = b->du; s->dv = b->dv; s->dw = b->dw;
    s->dim[0] = 3;
    s->dim[1] = b->nx; s->dim[2] = b->ny; s->dim[3] = b->nz;
    s->dim[4] = s->dim[5] = s->dim[6] = s->dim[7] = 1;
    s->pixdim[0] = b->pixdim[0];
    s->pixdim[1] = b->pixdim[1];
    s->pixdim[2] = b->pixdim[2];
    s->pixdim[3] = b->pixdim[3];
    s->pixdim[4] = b->pixdim[4];
    s->pixdim[5] = b->pixdim[5];
    s->pixdim[6] = b->pixdim[6];
    s->pixdim[7] = b->pixdim[7];
    s->nvox = (size_t)b->nx * b->ny * b->nz;
    s->sform_code = b->sform_code; s->sto_xyz = b->sto_xyz; s->sto_ijk = b->sto_ijk;
    s->qform_code = b->qform_code; s->qto_xyz = b->qto_xyz; s->qto_ijk = b->qto_ijk;
    s->quatern_b = b->quatern_b; s->quatern_c = b->quatern_c; s->quatern_d = b->quatern_d;
    s->qoffset_x = b->qoffset_x; s->qoffset_y = b->qoffset_y; s->qoffset_z = b->qoffset_z;
    s->qfac = b->qfac;
    s->xyz_units = b->xyz_units;
    s->time_units = b->time_units;
    s->freq_dim = b->freq_dim;
    s->phase_dim = b->phase_dim;
    s->slice_dim = b->slice_dim;
    s->slice_code = b->slice_code;
    s->slice_start = b->slice_start;
    s->slice_end = b->slice_end;
    s->slice_duration = b->slice_duration;
}

/* Validate a 3D image for resampling: an overflow-safe voxel count that fits int
 * (via al_safe_nvox) and that the 3D product equals nim->nvox (a genuine single-volume
 * 3D image, not 4D). Returns 0 if usable. Guards malformed/extreme NIfTI-2 headers. */
static int al_dims_ok(const nifti_image *n, const char *who)
{
    size_t nv;
    if (al_safe_nvox(n, who, &nv)) return 1;
    if ((long long)n->nvox != (long long)nv) {
        fprintf(stderr, "%s: 4D/multi-volume input not supported (3D only)\n", who); return 1; }
    return 0;
}

/* Reslice `source` onto `base`'s grid using an explicit base-index -> source-index
 * affine `gam` (0-based NIfTI voxel indices). Replaces source->data with the
 * resliced float volume and adopts base dims + sform/qform. interp: AL_INTERP_*
 * (NN/LINEAR/CUBIC); out-of-FOV voxels set to `fillv` (0 or NaN). Returns 0 on
 * success. The interpolation is the same BSD code -allineate uses; -spm_coreg
 * reuses it so the GPL module only computes the transform, never resamples. */
int nii_reslice_affine(nifti_image *source, const nifti_image *base,
                       mat44 gam, int interp, float fillv)
{
    if (source == NULL || base == NULL) { fprintf(stderr, "reslice: NULL image\n"); return 1; }
    if (al_dims_ok(source, "reslice source") || al_dims_ok(base, "reslice grid")) return 1;
    /* Fast path: alias float32 source data instead of copying via nii_to_float.
     * Only safe when no intensity scaling is pending — nii_to_float bakes in
     * scl_slope/scl_inter, so a scaled float32 source (e.g. a raw-read -spm_deface
     * mask) must take the converting path or the reslice would use raw values. */
    int aim_alias = (source->datatype == DT_FLOAT32 && source->data != NULL &&
                     (source->scl_slope == 0.0f ||
                      (source->scl_slope == 1.0f && source->scl_inter == 0.0f)));
    float *aim = aim_alias ? (float *)source->data : nii_to_float(source);
    if (!aim) { fprintf(stderr, "reslice: cannot extract source data\n"); return 1; }
    int anx = source->nx, any = source->ny, anz = source->nz;
    int bnx = base->nx, bny = base->ny, bnz = base->nz;
    size_t bnvox = (size_t)bnx * bny * bnz;
    float *out = (float *)malloc(sizeof(float) * bnvox);
    if (!out) { if (!aim_alias) free(aim); fprintf(stderr, "reslice: out of memory\n"); return 1; }
    int nxy_a = anx * any;
    float nxh = anx - 0.501f, nyh = any - 0.501f, nzh = anz - 0.501f;
    int anx1 = anx - 1, any1 = any - 1, anz1 = anz - 1;
    float mx0=gam.m[0][0], mx1=gam.m[0][1], mx2=gam.m[0][2], mx3=gam.m[0][3];
    float my0=gam.m[1][0], my1=gam.m[1][1], my2=gam.m[1][2], my3=gam.m[1][3];
    float mz0=gam.m[2][0], mz1=gam.m[2][1], mz2=gam.m[2][2], mz3=gam.m[2][3];
    /* Cubic interpolation can overshoot the source range (ringing); clamp each
     * in-FOV interpolated value to [clo,chi]. This must NOT touch out-of-FOV
     * voxels: those carry the caller's requested `fillv`, and clamping them to a
     * positive source minimum would corrupt the fill (turning "remove"/0 into the
     * source min — a defacing mask-safety hazard). So the clamp is applied inline
     * to the interpolated value only, never to `fillv`. */
    int clip_cubic = (interp == AL_INTERP_CUBIC);
    float clo = 0.0f, chi = 0.0f;
    if (clip_cubic) {
        size_t nsrc = (size_t)anx * any * anz;
        clo = chi = aim[0];
        for (size_t i = 1; i < nsrc; i++) {
            if (aim[i] < clo) clo = aim[i];
            else if (aim[i] > chi) chi = aim[i];
        }
    }
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) if(bnvox > 100000)
#endif
    for (int kk = 0; kk < bnz; kk++) {
        for (int jj = 0; jj < bny; jj++) {
            size_t ob = ((size_t)kk * bny + jj) * bnx;
            float bx = mx1*jj + mx2*kk + mx3;
            float by = my1*jj + my2*kk + my3;
            float bz = mz1*jj + mz2*kk + mz3;
            for (int ii = 0; ii < bnx; ii++) {
                float xx = mx0*ii + bx, yy = my0*ii + by, zz = mz0*ii + bz;
                float v = fillv; int oob = 0;
                if (interp == AL_INTERP_NN) {
                    if (xx<-0.499f||xx>nxh||yy<-0.499f||yy>nyh||zz<-0.499f||zz>nzh) oob = 1;
                    else { int ix=(int)(xx+0.5f), jy=(int)(yy+0.5f), kz=(int)(zz+0.5f);
                           if(ix<0)ix=0; else if(ix>anx1)ix=anx1;
                           if(jy<0)jy=0; else if(jy>any1)jy=any1;
                           if(kz<0)kz=0; else if(kz>anz1)kz=anz1;
                           v = aim[ix + jy*anx + (long)kz*nxy_a]; }
                } else if (interp == AL_INTERP_CUBIC) {
                    v = al_interp_cubic_checked(xx,yy,zz, nxh,nyh,nzh, anx1,any1,anz1, aim,anx,nxy_a, &oob);
                    if (!oob) { if (v < clo) v = clo; else if (v > chi) v = chi; }
                } else {
                    v = al_interp_linear_checked(xx,yy,zz, nxh,nyh,nzh, anx1,any1,anz1, aim,anx,nxy_a, &oob);
                }
                out[ob+ii] = oob ? fillv : v;
            }
        }
    }
    if (!aim_alias) free(aim);
    free(source->data);
    source->data = out;
    al_adopt_geometry(source, base);
    return 0;
}

/* Apply a saved world-space FIXED->MOVING affine (as written by -savemat's
   `fixed_to_moving`) to reslice `input` (an image in the MOVING/source space of the
   prior registration) onto `target`'s grid. Because `fixed_to_moving` is world-mm ->
   world-mm, `target` may have any resolution / FOV / origin that shares the fixed
   world frame. Builds the index->index map gam = S_input^{-1} * fixed_to_moving *
   S_target and calls nii_reslice_affine (which replaces input->data with the resliced
   volume and adopts target geometry). interp: AL_INTERP_*; out-of-FOV -> fillv.
   Returns 0 on success, nonzero on error. */
int nii_apply_affine(nifti_image *input, const nifti_image *target,
                     mat44 fixed_to_moving, int interp, float fillv)
{
    if (input == NULL || target == NULL) { fprintf(stderr, "applymat: NULL image\n"); return 1; }
    if (!al_mat44_usable(fixed_to_moving)) {
        fprintf(stderr, "applymat: matrix is not finite/invertible\n");
        return 1;
    }
    /* A valid affine's homogeneous bottom row is [0,0,0,1]. A corrupt or hand-edited
       -savemat matrix (e.g. last row [1,2,3,4]) can be finite with an invertible 3x3 yet
       is not an affine transform — it would silently produce garbage (an all-zero image).
       Reject it here, once, at the apply boundary (not duplicated in the JSON parser).
       Internally-computed transforms — GA_setup_affine / al_world_xform_from_params /
       coreg_fast / seeded compositions of vox->world matrices — all carry [0,0,0,1]. */
    /* al_mat44_usable only checks the upper 3 rows for finiteness, and `fabsf(nan) > x` is
       false — so validate the bottom row's finiteness explicitly (a NaN there would otherwise
       pass and write an all-NaN image) as well as the [0,0,0,1] affine pattern. */
    if (!al_finitef(fixed_to_moving.m[3][0]) || !al_finitef(fixed_to_moving.m[3][1]) ||
        !al_finitef(fixed_to_moving.m[3][2]) || !al_finitef(fixed_to_moving.m[3][3]) ||
        fabsf(fixed_to_moving.m[3][0]) > 1e-4f || fabsf(fixed_to_moving.m[3][1]) > 1e-4f ||
        fabsf(fixed_to_moving.m[3][2]) > 1e-4f || fabsf(fixed_to_moving.m[3][3] - 1.0f) > 1e-4f) {
        fprintf(stderr, "applymat: matrix bottom row is not a finite [0,0,0,1] (not an affine transform)\n");
        return 1;
    }
    mat44 S_in, S_tg;
    al_image_xform_or_pixdim(input, &S_in, NULL);   /* no-form -> pixdim frame (unified policy) */
    al_image_xform_or_pixdim(target, &S_tg, NULL);
    mat44 gam = nifti_mat44_mul(nifti_mat44_inverse(S_in),
                                nifti_mat44_mul(fixed_to_moving, S_tg));
    return nii_reslice_affine(input, target, gam, interp, fillv);
}

static mat44 al_world_xform_from_params(int npar, float *wpar);  /* defined below */

/* Estimate-only registration: fit the affine and return BOTH the raw 12-DOF warp
   params (wpar_out, for the fused final reslice) AND the world-mm FIXED(base)->MOVING
   (source) "pull" transform (world_out, == what -savemat/nii_apply_affine consume),
   WITHOUT reslicing or mutating either image. Prints the cost line; al_register prints
   its own progress. This is the single estimation seam shared by nii_allineate (which
   then reslices onto base), the -master wrapper, and future -savemat — so estimation and
   application are separable and the moving image is never resliced-then-discarded.
   Returns 0 on success, nonzero on failure (leaves *world_out untouched). */
static int al_estimate(nifti_image *source, nifti_image *base, al_opts opts,
                       float wpar_out[12], mat44 *world_out)
{
    if (source == NULL || base == NULL) {
        fprintf(stderr, "allineate: NULL input image\n");
        return 1;
    }
    /* 3D only. al_register uses volume 0 of a 4D image and the warped output is a
     * single 3D buffer, so a 4D source would write a header that advertises 4D over
     * 3D bytes (corrupt file); a 4D base is equally unsupported. Fail closed. */
    if (al_dims_ok(source, "allineate source") || al_dims_ok(base, "allineate base")) return 1;

    int match_code;
    const char *cost_name;
    al_resolve_cost(opts.cost, &match_code, &cost_name);
    fprintf(stderr, " + Cost function: %s, cmass: %s, warp: %s (%d DOF)\n", cost_name,
            opts.cmass ? "yes" : "no", al_warp_name(opts.warp), opts.warp);

    int ok = al_register(source, base, match_code, opts.cmass,
                         opts.source_automask, opts.dark_automask, opts.zoom /*relax_scale*/,
                         opts.interp, opts.warp, wpar_out, NULL);
    if (ok) return ok;
    *world_out = al_world_xform_from_params(12, wpar_out);
    return 0;
}

/* Public estimate-only entry: fit source->base and return the world-mm FIXED->MOVING
   affine without touching image data. Used by -master (estimate once, apply once onto
   the requested grid) and available for -savemat. Serial-only. Returns 0 on success. */
int nii_allineate_estimate(nifti_image *source, nifti_image *base, al_opts opts,
                           mat44 *fixed_to_moving)
{
    if (fixed_to_moving == NULL) { fprintf(stderr, "allineate: NULL output matrix\n"); return 1; }
    float wpar[12];
    return al_estimate(source, base, opts, wpar, fixed_to_moving);
}

int nii_allineate(nifti_image *source, nifti_image *base, al_opts opts)
{
    double t_start = al_wtime();
    g_last_affine_valid = 0;   /* invalidate up front: nii_last_affine() must reflect only
                                  a COMPLETED fit, never a stale one or a partial failure */
    float wpar[12];
    mat44 world;
    int ok = al_estimate(source, base, opts, wpar, &world);
    if (ok) return ok;

    /* Extract source float data for warping */
    float *ajim = nii_to_float(source);
    if (!ajim) {
        fprintf(stderr, "allineate: failed to extract source data for warp\n");
        return 1;
    }

    int anx = source->nx, any = source->ny, anz = source->nz;
    int bnx = base->nx, bny = base->ny, bnz = base->nz;

    /* Apply final warp (default: cubic for registration) */
    int final_ic = (opts.final_interp == AL_INTERP_DEFAULT) ? AL_INTERP_CUBIC : opts.final_interp;
    const char *interp_name = (final_ic == AL_INTERP_NN) ? "nearest" :
                              (final_ic == AL_INTERP_LINEAR) ? "linear" : "cubic";
    long reg_ms = (long)((al_wtime() - t_start) * 1000.0 + 0.5);
    fprintf(stderr, " + Registration completed in %ldms (coarse %ldms; fine %ldms)\n",
            reg_ms, (long)(g_last_coarse_ms + 0.5), (long)(g_last_fine_ms + 0.5));
#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
    if (nthreads > 1)
        fprintf(stderr, " + Applying final warp with %s interpolation (%d threads)\n", interp_name, nthreads);
    else
#endif
        fprintf(stderr, " + Applying final warp with %s interpolation\n", interp_name);
    PROFILE_START(warp);
    float *warped = al_scalar_warpone(12, wpar, al_wfunc_affine,
                                      ajim, anx, any, anz,
                                      bnx, bny, bnz, final_ic);
    PROFILE_END(warp, "final warp");
    free(ajim);
    if (warped == NULL) {
        fprintf(stderr, "allineate: failed to warp image\n");
        return 1;
    }

    /* Replace source->data with warped result, adopt base dims/transforms */
    free(source->data);
    source->data = warped;
    al_adopt_geometry(source, base);

    /* Record the fitted world-space FIXED(base)->MOVING(source) affine only now that
       the whole fit+warp succeeded, so nii_last_affine()/-savemat never expose a
       partial result. `world` was computed by al_estimate (== al_world_xform_from_params
       (12, wpar)). g_last_affine is retained for the shared standalone allineate project's
       wired -savemat (which reads it via nii_last_affine); niimath's -master no longer
       uses it — it takes the matrix directly from nii_allineate_estimate. */
    g_last_affine = world;
    g_last_affine_valid = 1;
    return 0;
}

/* Copy the most recent nii_allineate() fit's world-space FIXED(base)->MOVING(source)
   affine into *out. Returns 0 if a fit has run (valid), nonzero otherwise. Reflects
   the moving image as passed to registration (after any -com/-sym header fold).
   Serial-only (reads a process-global; see nii_allineate). */
int nii_last_affine(mat44 *out)
{
    if (out == NULL || !g_last_affine_valid) return 1;
    *out = g_last_affine;
    return 0;
}

/* Zero (to the image minimum) every input voxel whose warped mask < 0.5. The
   warped mask is sampled on input's grid (e.g. from al_scalar_warpone or
   nii_reslice_affine). Ensures input is float32 first. BSD; shared by -deface
   and the GPL -spm_deface. Returns #voxels masked, or -1 on error. */
long nii_apply_deface_mask(nifti_image *input, const float *warped_mask)
{
    /* Exported shared helper (BSD -deface + GPL -spm_deface): validate the
     * preconditions callers rely on rather than trusting them. */
    if (input == NULL || warped_mask == NULL) { fprintf(stderr, "deface: NULL input/mask\n"); return -1; }
    if (al_dims_ok(input, "deface apply")) return -1;
    if (al_ensure_float32(input, "deface")) return -1;
    size_t nvox = (size_t)input->nx * input->ny * input->nz;
    float *idata = (float *)input->data;
    if (input->data == NULL) { fprintf(stderr, "deface: input has no data\n"); return -1; }
    /* Finite minimum: a NaN-first voxel must not poison every masked voxel with NaN.
     * Magnitude guard (not isfinite()) because this TU is built with -ffast-math;
     * `-FLT_MAX <= v <= FLT_MAX` excludes NaN and ±inf. Fallback 0 if all nonfinite. */
    float min_val = 0.0f; int have_min = 0;
    for (size_t i = 0; i < nvox; i++) {
        float v = idata[i];
        if (al_finitef(v) && (!have_min || v < min_val)) { min_val = v; have_min = 1; }
    }
    long nmasked = 0;
    long nv = (long)nvox;
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) reduction(+:nmasked) if(nv > 100000)
#endif
    for (long i = 0; i < nv; i++) {
        float m = warped_mask[i];
        if (!(m >= 0.5f && m <= FLT_MAX)) { idata[i] = min_val; nmasked++; }
    }
    return nmasked;
}

/* Deface: register input to template, INVERT the transform, warp the
   template-space mask onto the input's native grid, zero masked voxels.
   input: the image to deface (modified in-place, stays in its own native space)
   tmpl: template image — the registration FIXED/base (input is the moving image,
         the well-posed direction; the transform is then inverted for the mask warp)
   mask: mask in template space (>=0.5 = keep, <0.5 = remove). The mask, not any
         command name, determines what is removed (brain mask -> keep brain;
         face mask -> remove face).
   opts: registration options (cost function, cmass)
   Returns 0 on success, nonzero on error. */
int nii_deface(nifti_image *input, nifti_image *tmpl, nifti_image *mask, al_opts opts)
{
    double t_start = al_wtime();
    const char *label = "Deface";
    if (input == NULL || tmpl == NULL || mask == NULL) {
        fprintf(stderr, "%s: NULL input image\n", label);
        return 1;
    }
    /* 3D only — for the subject AND the template/mask. The face mask covers a
     * single volume (nx*ny*nz); a 4D input would zero only volume 0 and silently
     * leave identifiable faces in volumes 1..N while reporting success (a privacy
     * failure). A 4D template/mask would silently register/apply only their volume
     * 0 (al_register / nii_to_float use volume 0). Reject all three; fail closed. */
    if (al_dims_ok(input, label) || al_dims_ok(tmpl, "deface template")
        || al_dims_ok(mask, "deface mask")) return 1;
    if (al_ensure_float32(input, label)) return 1;

    int match_code;
    const char *cost_name;
    al_resolve_cost(opts.cost, &match_code, &cost_name);

    /* Default: linear for mask warping (cubic can ring). Out-of-FOV -> 0 so
     * anything the mask does not cover is treated as "remove". */
    int mask_interp = (opts.final_interp == AL_INTERP_DEFAULT) ? AL_INTERP_LINEAR : opts.final_interp;
    const char *minterp_name = mask_interp == AL_INTERP_NN ? "NN" :
                               mask_interp == AL_INTERP_CUBIC ? "cubic" : "linear";

    /* Register the SUBJECT (moving) to the TEMPLATE (fixed) — the well-posed
     * direction, identical to -allineate. It is robust because the template is
     * the base: a clean target image with a good autoweight. Registering the
     * brain-only template ONTO the full-head subject (base = subject) converges
     * poorly and mislocates the mask, because the subject's neck/shoulders/FOV
     * dominate the cost. We then INVERT the transform to pull the template-space
     * mask back onto the subject's native grid (so brain voxels are never
     * interpolated — the subject stays in its own space). Both engines reslice
     * the mask onto input's grid IN PLACE (mask->data becomes float on the input
     * lattice), then the common nii_apply_deface_mask() below zeros the masked
     * voxels. */
    if (opts.fast) {
        /* Fast SPM/FLIRT-inspired engine (default for -deface). coreg_fast_estimate
         * returns the world-mm FIXED(template)->MOVING(subject) transform WITHOUT
         * mutating inputs; INVERT it to pull the template-space mask onto the subject
         * grid (nii_apply_affine consumes a world FIXED->MOVING affine and reslices its
         * first arg — here the mask, treated as the FIXED-space image — onto the second
         * arg's grid, so the inverse is the moving->fixed direction it needs). */
        coreg_fast_opts cfo = coreg_fast_opts_default();
        cfo.cost = (opts.fast == AL_ENGINE_FAST_HEL) ? CF_COST_HEL : CF_COST_CR;
        /* -nocmass forces the supplied affine; otherwise auto-select (the fast engine
         * scores the supplied-affine and COM starts and keeps the better one). -deface
         * exposes no -com seed, so there is no recentered-header case here. */
        cfo.use_cmass = !((opts.cli_set & AL_CLI_CMASS) && opts.cmass == AL_CMASS_NONE);
        const char *fname = (opts.fast == AL_ENGINE_FAST_HEL) ? "fast (Hellinger)"
                                                              : "fast (correlation-ratio)";
        fprintf(stderr, " + %s: registering input to template using %s (mask warped back to native space)\n", label, fname);
        coreg_fast_result res;
        if (coreg_fast_estimate(input, tmpl, &cfo, &res)) {
            fprintf(stderr, "%s: fast registration failed\n", label);
            return 1;
        }
        fprintf(stderr, " + Warping mask to input space (%s interpolation)\n", minterp_name);
        if (nii_apply_affine(mask, input, nifti_mat44_inverse(res.fixed_to_moving), mask_interp, 0.0f)) {
            fprintf(stderr, "%s: failed to warp mask\n", label);
            return 1;
        }
    } else {
        fprintf(stderr, " + %s: registering input to template using %s (mask warped back to native space)\n", label, cost_name);
        float wpar[12];
        int ok = al_register(input, tmpl, match_code, opts.cmass,
                             opts.source_automask, opts.dark_automask, 0 /*relax_scale*/,
                             opts.interp, opts.warp, wpar, NULL);
        if (ok) {
            fprintf(stderr, "%s: registration failed\n", label);
            return 1;
        }
        /* gam maps template(base) voxels -> subject(source) voxels; invert to get
         * subject -> template, the pull-warp that resamples the mask onto the input
         * grid. nifti_mat44_inverse handles the translation/origin offset. */
        mat44 gam = GA_setup_affine(12, wpar);
        mat44 gam_inv = nifti_mat44_inverse(gam);
        fprintf(stderr, " + Warping mask to input space (%s interpolation)\n", minterp_name);
        if (nii_reslice_affine(mask, input, gam_inv, mask_interp, 0.0f)) {
            fprintf(stderr, "%s: failed to warp mask\n", label);
            return 1;
        }
    }

    int inx = input->nx, iny = input->ny, inz = input->nz;
    int nvox = inx * iny * inz;

    long nmasked = nii_apply_deface_mask(input, (const float *)mask->data);
    if (nmasked < 0) return 1;
    long elapsed_ms = (long)((al_wtime() - t_start) * 1000.0 + 0.5);
#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
    if (nthreads > 1)
        fprintf(stderr, " + %s complete: %ld of %d voxels masked (%.1f%%; %ldms, %d threads)\n",
                label, nmasked, nvox, 100.0f * nmasked / nvox, elapsed_ms, nthreads);
    else
#endif
        fprintf(stderr, " + %s complete: %ld of %d voxels masked (%.1f%%; %ldms)\n",
                label, nmasked, nvox, 100.0f * nmasked / nvox, elapsed_ms);
    return 0;
}

/*==========================================================================*/
/*=================== -sym: midsagittal alignment (MSP) ====================*/
/*==========================================================================*/

/* Compute H = T^(1/2) of a rigid mat44 T (rotation + translation), so H*H == T.
   Rotation: quaternion half-angle, q_half = normalize(q + identity). Translation:
   solve (R_half + I) * t_half = t. The two guards below are defensive: for a
   proper-rotation input (which GA_setup_affine(6,·) always yields, scalar part
   qw >= 0) neither can actually fire — hw = qw+1 >= 1 so nrm >= 1, and R_half+I is
   non-singular for any real rotation. Returns 0 on success, 1 if degenerate (the
   caller falls back to the identity correction). */
static int al_half_transform(mat44 T, mat44 *H_out)
{
    /* Pure-rotation matrix from T's upper-left 3x3. */
    mat44 R; memset(&R, 0, sizeof(R)); R.m[3][3] = 1.0f;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) R.m[i][j] = T.m[i][j];

    float qb, qc, qd, qx, qy, qz, dx, dy, dz, qfac;
    nifti_mat44_to_quatern(R, &qb, &qc, &qd, &qx, &qy, &qz, &dx, &dy, &dz, &qfac);
    float qw = 1.0f - qb*qb - qc*qc - qd*qd;
    qw = (qw > 0.0f) ? sqrtf(qw) : 0.0f;

    /* Half-angle quaternion = normalize(q + identity). */
    float hw = qw + 1.0f, hb = qb, hc = qc, hd = qd;
    float nrm = sqrtf(hw*hw + hb*hb + hc*hc + hd*hd);
    if (nrm <= 1.0e-7f) return 1;   /* defensive: unreachable for qw >= 0 */
    float inv = 1.0f / nrm;
    hb *= inv; hc *= inv; hd *= inv;
    mat44 Rhalf = nifti_quatern_to_mat44(hb, hc, hd, 0.0f, 0.0f, 0.0f,
                                         1.0f, 1.0f, 1.0f, qfac);

    /* M = R_half + I (3x3, zero translation). */
    mat44 M; memset(&M, 0, sizeof(M)); M.m[3][3] = 1.0f;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) M.m[i][j] = Rhalf.m[i][j];
        M.m[i][i] += 1.0f;
    }
    float det = M.m[0][0]*(M.m[1][1]*M.m[2][2] - M.m[1][2]*M.m[2][1])
              - M.m[0][1]*(M.m[1][0]*M.m[2][2] - M.m[1][2]*M.m[2][0])
              + M.m[0][2]*(M.m[1][0]*M.m[2][1] - M.m[1][1]*M.m[2][0]);
    if (fabsf(det) < 1.0e-6f) return 1;

    mat44 Minv = nifti_mat44_inverse(M);   /* zero translation → pure linear */
    float hx, hy, hz;
    mat44_vec(Minv, T.m[0][3], T.m[1][3], T.m[2][3], &hx, &hy, &hz);

    mat44 H = Rhalf;
    H.m[0][3] = hx; H.m[1][3] = hy; H.m[2][3] = hz;
    *H_out = H;
    return 0;
}

/* Fold a world-space correction C into nim's header as an initial-estimate seed
   for a subsequent registration. `S` is nim's currently-selected index->world
   transform (the one C multiplies on the left: newS = C*S). sform is always
   updated; qform is updated only when the input carried a valid qform
   (qform_code > 0) — an sform-only input stays sform-only (no synthesized qform).
   Shared by the -sym and -sagseed pre-steps. */
/* Extract the pure world-space transform from fitted warp params: GA_setup_affine
   normally folds the base/target befafter matrices into its result (index-space
   warp); -sym and -sagseed want the bare world-space parameter matrix, so the
   process-global aff_use_before/after must be transiently disabled around the call.
   Centralized here so the subtle save/restore has one home. */
static mat44 al_world_xform_from_params(int npar, float *wpar)
{
    int sb = aff_use_before, sa = aff_use_after;
    aff_use_before = aff_use_after = 0;
    mat44 M = GA_setup_affine(npar, wpar);
    aff_use_before = sb; aff_use_after = sa;
    return M;
}

/* Write world-space transform M into nim's sform (sto_xyz/ijk); the caller sets
   sform_code per its own policy. Also syncs voxel sizes (dx/dy/dz + pixdim) to M's
   column norms so pixdim stays consistent when M carries a scale (e.g. a -zoom
   fold) — a no-op for a pure rigid M (column norms unchanged). Shared by
   al_fold_correction_into_header and al_deoblique_frame. */
static void al_write_sform(nifti_image *nim, mat44 M)
{
    nim->sto_xyz = mat44_to_dmat44(M);
    nim->sto_ijk = nifti_dmat44_inverse(nim->sto_xyz);
    nim->dx = nim->pixdim[1] = mat44_colnorm(M, 0);
    nim->dy = nim->pixdim[2] = mat44_colnorm(M, 1);
    nim->dz = nim->pixdim[3] = mat44_colnorm(M, 2);
}

/* Write world-space transform M into nim's qform (quaternion + qto_xyz/ijk); the
   caller sets qform_code per its own policy. Pair to al_write_sform (which syncs
   pixdim, the field the serialized qform scale is derived from). */
static void al_write_qform(nifti_image *nim, mat44 M)
{
    float qb, qc, qd, qx, qy, qz, dx, dy, dz, qfac;
    nifti_mat44_to_quatern(M, &qb, &qc, &qd, &qx, &qy, &qz, &dx, &dy, &dz, &qfac);
    nim->quatern_b = qb; nim->quatern_c = qc; nim->quatern_d = qd;
    nim->qoffset_x = qx; nim->qoffset_y = qy; nim->qoffset_z = qz;
    nim->qfac = qfac;
    nim->qto_xyz = mat44_to_dmat44(M);
    nim->qto_ijk = nifti_dmat44_inverse(nim->qto_xyz);
}

static void al_fold_correction_into_header(nifti_image *nim, mat44 C, mat44 S)
{
    al_write_sform(nim, nifti_mat44_mul(C, S));
    if (nim->sform_code <= 0) nim->sform_code = 1;
    /* Update the qform only when the input carried a *usable* qform. A coded but
       degenerate/NaN qform (the reader does not validate srow/quaternion) would
       otherwise fold C into garbage and write back an invalid coded qform; drop
       its code instead so downstream selectors fall through to the folded sform. */
    if (nim->qform_code > 0) {
        if (al_mat44_usable(dmat44_to_mat44(nim->qto_xyz)))
            al_write_qform(nim, nifti_mat44_mul(C, dmat44_to_mat44(nim->qto_xyz)));
        else
            nim->qform_code = 0;   /* unusable coded qform: drop it, keep folded sform */
    }
}

/* De-oblique: snap the image's index->world transform to the nearest axis-aligned
   frame — a signed permutation of the voxel axes scaled by the voxel sizes — while
   keeping the grid-center world position fixed. Header-only (voxel data untouched);
   both sform and qform are rewritten to the snapped frame. This treats the
   ACQUISITION GRID as the anatomical frame, appropriate for scans deliberately
   acquired oblique-to-world (AC-PC angling, tilted infant/kyphotic positioning)
   where the voxel axes are the anatomical axes and the oblique sform merely records
   the scanner angle. Used by -symd so the mirror fit is not fooled into rotating an
   already-grid-symmetric head onto the oblique world frame. */
static void al_deoblique_frame(nifti_image *nim)
{
    mat44 S;
    al_image_xform_or_pixdim(nim, &S, NULL);

    float vs[3];
    for (int j = 0; j < 3; j++)
        vs[j] = sqrtf(S.m[0][j]*S.m[0][j] + S.m[1][j]*S.m[1][j] + S.m[2][j]*S.m[2][j]);

    /* Snap S's 3x3 to the nearest ORIENTATION-PRESERVING signed permutation of the
       voxel axes. Score all six permutations and fold the handedness constraint INTO
       the score: a permutation whose natural per-column signs give the wrong
       determinant sign must flip one column to stay proper, costing 2x that column's
       alignment confidence, so its best achievable signed score is
       (unsigned_score - 2*weakest_conf). Choosing the max over these *penalized* scores
       — rather than the max unsigned score with a post-selection flip — yields the
       globally nearest proper frame: a strong unsigned match with the wrong handedness
       can lose to a slightly weaker permutation that is already correctly oriented.
       (det S == 0: degenerate input, no handedness to preserve, so never penalize.) */
    static const int   P[6][3]  = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    static const float Ppar[6]  = { 1.0f, -1.0f, -1.0f, 1.0f, 1.0f, -1.0f }; /* perm parity */
    float detS = S.m[0][0] * (S.m[1][1]*S.m[2][2] - S.m[1][2]*S.m[2][1])
               - S.m[0][1] * (S.m[1][0]*S.m[2][2] - S.m[1][2]*S.m[2][0])
               + S.m[0][2] * (S.m[1][0]*S.m[2][1] - S.m[1][1]*S.m[2][0]);
    float sdetS = (detS > 0.0f) ? 1.0f : (detS < 0.0f ? -1.0f : 0.0f);

    int bp = 0, bflip = -1; float bscore = -FLT_MAX;
    for (int p = 0; p < 6; p++) {
        float uscore = 0.0f, weakc = FLT_MAX, sprod = 1.0f; int weakj = 0;
        for (int j = 0; j < 3; j++) {
            int bk = P[p][j];
            float conf = (vs[j] > 0.0f) ? fabsf(S.m[bk][j]) / vs[j] : 0.0f;
            uscore += conf;
            sprod  *= (S.m[bk][j] >= 0.0f) ? 1.0f : -1.0f;
            if (conf < weakc) { weakc = conf; weakj = j; }
        }
        int flip = (sdetS != 0.0f && Ppar[p] * sprod * sdetS < 0.0f) ? weakj : -1;
        float sc = (flip >= 0) ? uscore - 2.0f * weakc : uscore;
        if (sc > bscore) { bscore = sc; bp = p; bflip = flip; }
    }

    float sgn[3];
    for (int j = 0; j < 3; j++)
        sgn[j] = (S.m[P[bp][j]][j] >= 0.0f) ? 1.0f : -1.0f;
    if (bflip >= 0) sgn[bflip] = -sgn[bflip];   /* restore input handedness */

    mat44 Sd; memset(&Sd, 0, sizeof(Sd)); Sd.m[3][3] = 1.0f;
    for (int j = 0; j < 3; j++)
        Sd.m[P[bp][j]][j] = sgn[j] * (vs[j] > 0.0f ? vs[j] : 1.0f);

    /* Keep the grid-center voxel at the same world position it had under S. */
    float cx = (nim->nx - 1) * 0.5f, cy = (nim->ny - 1) * 0.5f, cz = (nim->nz - 1) * 0.5f;
    float wx, wy, wz, dx, dy, dz;
    mat44_vec(S,  cx, cy, cz, &wx, &wy, &wz);
    mat44_vec(Sd, cx, cy, cz, &dx, &dy, &dz);   /* Sd translation still 0 here */
    Sd.m[0][3] = wx - dx; Sd.m[1][3] = wy - dy; Sd.m[2][3] = wz - dz;

    al_write_sform(nim, Sd);   /* also syncs pixdim to the snapped column norms */
    nim->sform_code = 1;
    al_write_qform(nim, Sd);
    nim->qform_code = 1;
    fprintf(stderr, " + [sym] de-obliqued starting frame (voxel grid treated as anatomical)\n");
}

/* Run the mirror-registration fit in nim's CURRENT frame and return the world-space
   MSP correction C (= T^{-1/2}) and the selected index->world S (both used by the
   caller to apply). Does not modify nim's data or header. Returns 0 on success. */
static int al_sym_correction(nifti_image *nim, int dark_automask, mat44 *C_out, mat44 *S_out)
{
    mat44 S;
    al_image_xform_or_pixdim(nim, &S, "sym");

    /* Diagnostic: which voxel axis is most left-right (largest |world-X| column). */
    {
        int lr = 0; float best = fabsf(S.m[0][0]);
        for (int j = 1; j < 3; j++) if (fabsf(S.m[0][j]) > best) { best = fabsf(S.m[0][j]); lr = j; }
        fprintf(stderr, " + [sym] left-right voxel axis: %d (world-X polarity %+.0f)\n",
                lr, S.m[0][lr] >= 0.0f ? 1.0f : -1.0f);
    }

    /* Source view: pin the selected frame S onto a shallow copy (al_register errors
       on an image with no usable sform/qform; the pixdim fallback above needs this). */
    nifti_image src = *nim;
    src.sform_code = 1;
    src.sto_xyz = mat44_to_dmat44(S);
    src.sto_ijk = nifti_dmat44_inverse(src.sto_xyz);
    src.qform_code = 0;

    /* Mirror view: same voxel data, world-X reflected sform (F = diag(-1,1,1)). */
    mat44 F = mat44_diag(-1.0f, 1.0f, 1.0f);
    mat44 FS = nifti_mat44_mul(F, S);
    nifti_image mir = *nim;
    mir.sform_code = 1;
    mir.sto_xyz = mat44_to_dmat44(FS);
    mir.sto_ijk = nifti_dmat44_inverse(mir.sto_xyz);
    mir.qform_code = 0;

    /* Clamp the mirror-registration shift to a quarter of the longest FOV. */
    float fov_x = nim->nx * (float)fabs(nim->dx);
    float fov_y = nim->ny * (float)fabs(nim->dy);
    float fov_z = nim->nz * (float)fabs(nim->dz);
    float fov = fov_x; if (fov_y > fov) fov = fov_y; if (fov_z > fov) fov = fov_z;
    al_shift_max_override = fov * 0.25f;

    fprintf(stderr, " + [sym] mirror registration (ls, 6 DOF, cmass, shift<=%.1fmm)\n",
            al_shift_max_override);
    float wpar[12];
    int ok = al_register(&src, &mir, GA_MATCH_PEARSON_SCALAR, 1 /*cmass*/,
                         0 /*source_automask*/, dark_automask, 0 /*relax_scale*/,
                         AL_INTERP_LINEAR, 6 /*warp*/, wpar, NULL);
    al_shift_max_override = 0.0f;
    if (ok) { fprintf(stderr, "sym: mirror registration failed\n"); return 1; }

    /* Pure world transform T (base=mirror -> source=original); C = H^(-1) = T^(-1/2). */
    mat44 T = al_world_xform_from_params(6, wpar);
    mat44 H, C;
    if (al_half_transform(T, &H)) {
        fprintf(stderr, " + [sym] WARNING: half-transform degenerate; using identity correction\n");
        C = mat44_diag(1.0f, 1.0f, 1.0f);
    } else {
        C = nifti_mat44_inverse(H);
    }
    if (fabsf(wpar[3]) > 29.0f || fabsf(wpar[4]) > 29.0f || fabsf(wpar[5]) > 29.0f)
        fprintf(stderr, " + [sym] WARNING: fitted rotation near the ±30° guardrail;"
                        " MSP correction may be incomplete\n");
    *C_out = C; *S_out = S;
    return 0;
}

/* -symb tie tolerance (deg): de-oblique must beat the world frame's correction
   rotation by more than this to be chosen; otherwise the original frame is kept. */
#define AL_SYMB_ROT_TOL_DEG 0.5f

/* Rotation magnitude (deg) of a rigid correction's 3x3, for -symb frame selection. */
static float al_correction_rot_deg(mat44 M)
{
    float tr = M.m[0][0] + M.m[1][1] + M.m[2][2];
    float c = (tr - 1.0f) * 0.5f;
    if (c > 1.0f) c = 1.0f; else if (c < -1.0f) c = -1.0f;
    return (float)(acos((double)c) * 180.0 / 3.14159265358979323846);
}

/* Template-free midsagittal-plane (MSP) alignment.

   Registers the image (source) to its world-X mirror (base) with a 6-DOF rigid
   fit, takes the half transform H = T^(1/2) of the recovered rigid T, and forms
   the correction C = H^(-1) (= T^(-1/2)). The direction is proven by the
   translation-only phantom regression: a symmetric image shifted +10 mm in world
   X yields C = a -10 mm shift (re-centering the MSP to X = 0).

   nim:     image to reorient (must be 3D; float32 recommended).
   C_out:   if non-NULL, receives the 4x4 world-space correction C.
   reslice: nonzero -> resample the data so it is symmetric about world X = 0
            (standalone use). zero -> fold C into the header (sform always; qform
            only when the input carried a *usable* coded qform) as an
            initial-estimate seed for a subsequent registration (pre-step use).
   deoblique: 0 = -sym (fit in the image's own world frame); 1 = -symd (first snap the
            frame to axis-aligned, treating the voxel grid as anatomical, so an
            obliquely-acquired but grid-symmetric head is not rotated onto the oblique
            world frame — see al_deoblique_frame); 2 = -symb (auto-compete: fit BOTH
            frames and keep the one whose correction rotation is SMALLER — that frame's
            X=0 is already closer to the true MSP, i.e. it is the native anatomical
            frame, so -symb picks de-oblique for an oblique grid-symmetric head and the
            world frame when the sform genuinely encodes anatomical axes).
   The mirror fit uses the `ls`/Pearson cost (a Hellinger variant was tried and
   removed — it gave the same result; the roll it was meant to fix was frame
   obliquity, which -symd addresses).
   dark_automask: nonzero (-dark_automask) -> drop background/pad matched pairs
            (at the image minimum) from the mirror-fit cost.
   Returns 0 on success, nonzero on error. */
int nii_symmetry(nifti_image *nim, mat44 *C_out, int reslice, int deoblique, int dark_automask)
{
    if (nim == NULL) { fprintf(stderr, "sym: NULL image\n"); return 1; }
    if (al_dims_ok(nim, "sym")) return 1;

    /* De-oblique (modes 1/2) mutates nim's header up front; snapshot it so any
       failure — or the -symb decision to keep the world frame — rolls the header
       back. The data pointer is preserved across every restore, so nothing is
       double-freed. */
    nifti_image nim_hdr0 = *nim;
    mat44 C, S;

    if (deoblique == 2) {   /* -symb: auto-compete world vs de-obliqued frame */
        mat44 C_ob, S_ob;
        if (al_sym_correction(nim, dark_automask, &C_ob, &S_ob)) return 1;
        float rot_ob = al_correction_rot_deg(C_ob);

        al_deoblique_frame(nim);
        mat44 C_de, S_de;
        if (al_sym_correction(nim, dark_automask, &C_de, &S_de)) {
            void *_d = nim->data; *nim = nim_hdr0; nim->data = _d;   /* undo de-oblique */
            return 1;
        }
        float rot_de = al_correction_rot_deg(C_de);

        /* Prefer the ORIGINAL world frame unless de-oblique is meaningfully better by a
           tolerance. The frames' fit costs are near-tied (same data, mirror geometry),
           so correction rotation is the discriminator — but on an already axis-aligned
           image both rotations are identical, and switching to de-oblique there would
           needlessly rewrite the header (e.g. MNI sform_code 4 -> scanner-anatomical 1)
           for zero geometric gain. The tolerance makes the world frame win on ties. */
        if (rot_de < rot_ob - AL_SYMB_ROT_TOL_DEG) {
            C = C_de; S = S_de;   /* keep the de-obliqued frame (nim stays de-obliqued) */
            fprintf(stderr, " + [symb] chose de-obliqued frame (correction %.1f deg vs world %.1f deg)\n",
                    rot_de, rot_ob);
        } else {
            void *_d = nim->data; *nim = nim_hdr0; nim->data = _d;   /* restore world frame */
            C = C_ob; S = S_ob;
            fprintf(stderr, " + [symb] chose world frame (correction %.1f deg vs de-obliqued %.1f deg)\n",
                    rot_ob, rot_de);
        }
    } else {
        if (deoblique) al_deoblique_frame(nim);
        if (al_sym_correction(nim, dark_automask, &C, &S)) {
            if (deoblique) { void *_d = nim->data; *nim = nim_hdr0; nim->data = _d; }
            return 1;
        }
    }

    fprintf(stderr, " + [sym] correction matrix C:\n");
    for (int i = 0; i < 3; i++)
        fprintf(stderr, "     % 9.4f % 9.4f % 9.4f % 9.4f\n",
                C.m[i][0], C.m[i][1], C.m[i][2], C.m[i][3]);

    if (C_out) *C_out = C;

    if (reslice) {
        /* Standalone: resample so data is symmetric about world X = 0.
           gam maps output(base)-index -> input(source)-index = S^(-1) * C^(-1) * S. */
        mat44 Sinv = nifti_mat44_inverse(S);
        mat44 Cinv = nifti_mat44_inverse(C);   /* = H */
        mat44 gam = nifti_mat44_mul(Sinv, nifti_mat44_mul(Cinv, S));
        if (nii_reslice_affine(nim, nim, gam, AL_INTERP_CUBIC, 0.0f)) {
            fprintf(stderr, "sym: reslice failed\n");
            void *_d = nim->data; *nim = nim_hdr0; nim->data = _d;   /* undo any de-oblique */
            return 1;
        }
    } else {
        /* Pre-step seed: fold C into the header (sform always; qform iff the input
         * carried a valid qform). Shared with -sagseed. */
        al_fold_correction_into_header(nim, C, S);
    }
    return 0;
}

/* -sagseed: in-MSP rigid seed (the complement of -sym). Precondition: -sym has
   already folded its MSP correction into nim's header so the midsagittal plane
   sits at world X = 0, and `tmpl` is an MSP-aligned template. Runs a 3-DOF-
   constrained fit of nim to tmpl freeing ONLY the MSP-preserving isometries
   {y-shift(1), z-shift(2), pitch = x-rotation(4)} — exactly the rigid DOF -sym is
   blind to — and folds the resulting correction P^(-1) into nim's header, yielding
   a full-rigid seed for the subsequent unconstrained nii_allineate. No reslice, no
   standalone output; uses the user's -cost / -source_automask / -interp / -cmass.
   Uses the process-global registration workspaces (serial-only, see nii_allineate).
   Returns 0 on success, nonzero on error. */
int nii_sagseed(nifti_image *nim, nifti_image *tmpl, al_opts opts)
{
    if (nim == NULL || tmpl == NULL) { fprintf(stderr, "sagseed: NULL image\n"); return 1; }
    if (al_dims_ok(nim, "sagseed") || al_dims_ok(tmpl, "sagseed template")) return 1;

    /* nim's current (post-sym) index->world transform; C is folded as newS = C*S. */
    mat44 S;
    al_image_xform_or_pixdim(nim, &S, "sagseed");

    int match_code;
    const char *cost_name;
    al_resolve_cost(opts.cost, &match_code, &cost_name);

    /* Free the in-MSP isometries: y-shift, z-shift, x-rotation (pitch); with -zoom,
       also free param 6 (x-scale), which GA_setup_affine ties to y/z under
       al_zoom_isotropic -> a single global isotropic zoom (helps extreme size
       mismatches, e.g. an infant brain vs an adult template). */
    static const int free_mask_rigid[12] = { 0,1,1, 0,1,0, 0,0,0, 0,0,0 };
    static const int free_mask_zoom[12]  = { 0,1,1, 0,1,0, 1,0,0, 0,0,0 };
    const int *free_mask = opts.zoom ? free_mask_zoom : free_mask_rigid;
    fprintf(stderr, " + [sagseed] in-MSP seed (%s; free: y-shift, z-shift, pitch%s)\n",
            cost_name, opts.zoom ? ", global zoom" : "");

    /* Enable the isotropic-zoom tie for the seed fit AND the transform extraction
       below; reset immediately after so the main nii_allineate() fit is unaffected. */
    int prev_zoom = al_zoom_isotropic;
    if (opts.zoom) al_zoom_isotropic = 1;

    float wpar[12];
    int ok = al_register(nim, tmpl, match_code, opts.cmass,
                         opts.source_automask, opts.dark_automask, 0 /*relax_scale: seed uses the isotropic tie*/,
                         opts.interp, 6 /*warp: ignored w/ mask*/,
                         wpar, free_mask);
    if (ok) { al_zoom_isotropic = prev_zoom; fprintf(stderr, "sagseed: constrained registration failed\n"); return 1; }

    /* Pure world transform P (base=tmpl -> source=nim); the seed correction that
       aligns nim to tmpl is P^(-1) (full inverse — direct fit, no -sym half). */
    mat44 P = al_world_xform_from_params(12, wpar);
    al_zoom_isotropic = prev_zoom;
    mat44 C = nifti_mat44_inverse(P);

    if (opts.zoom)
        fprintf(stderr, " + [sagseed] recovered seed: y-shift %+.2f mm, z-shift %+.2f mm, pitch %+.2f deg, zoom %.3f\n",
                wpar[1], wpar[2], wpar[4], wpar[6]);
    else
        fprintf(stderr, " + [sagseed] recovered seed: y-shift %+.2f mm, z-shift %+.2f mm, pitch %+.2f deg\n",
                wpar[1], wpar[2], wpar[4]);

    al_fold_correction_into_header(nim, C, S);
    return 0;
}

/*==========================================================================*/
/*==================== -com: center-of-mass origin =========================*/
/*==========================================================================*/

/* Set the image origin to its brightness center of mass: compute the
   intensity-weighted centroid of the positive voxels (same estimator as -cmass),
   map it to world coordinates through the selected index->world transform, and
   fold the pure translation C = translate(-centroid_world) into the header so the
   centroid lands at world (0,0,0). Header-only (no reslice), like the -sym
   pre-step; a cheap origin reset that gives a good starting point for symmetric
   images, intended to run early (right after -robustfov). Template-free, so a
   pixdim-centered frame is the documented fallback when the input carries no usable
   sform/qform. Returns 0 on success, nonzero on error. */
int nii_center_of_mass(nifti_image *nim)
{
    if (nim == NULL) { fprintf(stderr, "com: NULL image\n"); return 1; }
    if (al_dims_ok(nim, "com")) return 1;

    /* S: index -> world (same selector/fallback as -sym; -com is template-free). */
    mat44 S;
    al_image_xform_or_pixdim(nim, &S, "com");

    float *fdata = nii_to_float(nim);
    if (!fdata) { fprintf(stderr, "com: failed to extract float data\n"); return 1; }
    double cx, cy, cz;
    al_center_of_mass(fdata, nim->nx, nim->ny, nim->nz, &cx, &cy, &cz);
    free(fdata);

    float wx, wy, wz;
    mat44_vec(S, (float)cx, (float)cy, (float)cz, &wx, &wy, &wz);
    fprintf(stderr, " + [com] brightness centroid: voxel (%.1f %.1f %.1f) -> world (%+.2f %+.2f %+.2f) mm\n",
            cx, cy, cz, wx, wy, wz);

    mat44 C = mat44_diag(1.0f, 1.0f, 1.0f);
    C.m[0][3] = -wx; C.m[1][3] = -wy; C.m[2][3] = -wz;

    al_fold_correction_into_header(nim, C, S);
    return 0;
}
