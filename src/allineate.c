/* Affine image registration engine - adapted from AFNI 3dAllineate (public domain)
   Original by RW Cox, National Institutes of Health.
   Ported to standalone NIfTI by the niimath project.

   This file implements affine (12 DOF) registration with lpc+ZZ cost function,
   twopass coarse-to-fine optimization, autoweight, source automask, and cubic
   output interpolation. The Powell/NEWUOA optimizer is in a separate file. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "nifti_io.h"
#include "allineate.h"

/* Convert nifti_dmat44 (double) to mat44 (float) */
static mat44 dmat44_to_mat44(nifti_dmat44 d) {
    mat44 f;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            f.m[i][j] = (float)d.m[i][j];
    return f;
}

/*==========================================================================*/
/*============================== CONSTANTS ================================*/
/*==========================================================================*/

#define AL_BIGVAL   1.e+38f
#define AL_SMAGIC   208921148

/* Cost function method codes */
#define GA_MATCH_PEARSON_SCALAR        1
#define GA_MATCH_SPEARMAN_SCALAR       2
#define GA_MATCH_KULLBACK_SCALAR       3
#define GA_MATCH_CORRATIO_SCALAR       4
#define GA_MATCH_NORMUTIN_SCALAR       5
#define GA_MATCH_JOINTENT_SCALAR       6
#define GA_MATCH_HELLINGER_SCALAR      7
#define GA_MATCH_CORRATIO_ADD          8
#define GA_MATCH_CORRATIO_UNS          9
#define GA_MATCH_PEARSON_SIGNED       10
#define GA_MATCH_PEARSON_LOCALS       11
#define GA_MATCH_PEARSON_LOCALA       12
#define GA_MATCH_LPC_MICHO_SCALAR     13
#define GA_MATCH_LPA_MICHO_SCALAR     14
#define GA_MATCH_METHNUM_SCALAR       14

/* Interpolation codes */
#define AL_INTERP_NN      0
#define AL_INTERP_LINEAR  1
#define AL_INTERP_CUBIC   3

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
#define D2R (3.14159265358979f/180.0f)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define CMAX 0.9999f /* max correlation magnitude to allow */

#define NPER 262144  /* max points to warp at once */

#define PARAM_MAXTRIAL 15 /* max number of trial parameter sets to save */

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
    unsigned char *bmask;  /* base mask */
    float *bwght;          /* base weights */
    int nmask;

    /* Source/target image */
    float *ajim;           /* source image data */
    float *ajims;          /* smoothed source */
    float *ajimor;         /* original source (unused here) */
    int anx, any, anz;     /* source dimensions */
    float adx, ady, adz;   /* source voxel sizes */
    float ajbot, ajtop, ajclip;
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

    /* 2D histogram for MI/NMI */
    int need_hist_setup;

    /* micho (lpc+ZZ) parameters */
    double micho_mi, micho_nmi, micho_crA, micho_hel, micho_ov;
    int micho_zfinal;
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
extern void powell_set_verbose(int v);

/*==========================================================================*/
/*========================== GLOBAL STATE =================================*/
/*==========================================================================*/

static GA_setup *gstup = NULL;  /* current setup for optimizer callback */

static int aff_use_before = 0, aff_use_after = 0;
static mat44 aff_before, aff_after;

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

/* Rank-order a float array, ties get average rank. Output overwrites input. */
static void rank_order_float(int n, float *a)
{
    int ii, ns, n1, ib;
    int   *b;
    float *c;
    float cs;

    if (a == NULL || n < 1) return;
    if (n == 1) { a[0] = 0.0f; return; }

    b = (int *)  malloc(sizeof(int)   * n);
    c = (float *)malloc(sizeof(float) * n);
    if (b == NULL || c == NULL) { free(b); free(c); return; }

    for (ii = 0; ii < n; ii++) c[ii] = b[ii] = ii;
    qsort_floatint(n, a, b);

    n1 = n - 1;
    for (ii = 0; ii < n1; ii++) {
        if (a[ii] == a[ii + 1]) {
            cs = 2 * ii + 1;
            ns = 2;
            ib = ii; ii++;
            while (ii < n1 && a[ii] == a[ii + 1]) { ii++; ns++; cs += ii; }
            for (cs /= ns; ib <= ii; ib++) c[ib] = cs;
        }
    }
    for (ii = 0; ii < n; ii++) a[b[ii]] = c[ii];
    free(c); free(b);
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
static void mat44_vec(mat44 m, float x, float y, float z,
                      float *ox, float *oy, float *oz)
{
    *ox = m.m[0][0]*x + m.m[0][1]*y + m.m[0][2]*z + m.m[0][3];
    *oy = m.m[1][0]*x + m.m[1][1]*y + m.m[1][2]*z + m.m[1][3];
    *oz = m.m[2][0]*x + m.m[2][1]*y + m.m[2][2]*z + m.m[2][3];
}

/* Load a diagonal mat44 */
static mat44 mat44_diag(float a, float b, float c)
{
    mat44 m;
    memset(&m, 0, sizeof(mat44));
    m.m[0][0] = a; m.m[1][1] = b; m.m[2][2] = c; m.m[3][3] = 1.0f;
    return m;
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
static void gaussian_blur_3d(float *data, int nx, int ny, int nz,
                             float dx, float dy, float dz, float sigma)
{
    int krad, nxy, ii, jj, kk, rr;
    float *kernel, *buf;
    float sum, sigma_v, ksum;

    if (data == NULL || sigma <= 0.0f) return;
    nxy = nx * ny;

    /* Blur along X */
    sigma_v = sigma / dx;
    if (sigma_v > 0.25f) {
        krad = (int)(4.0f * sigma_v + 0.5f);
        if (krad < 1) krad = 1;
        kernel = (float *)calloc(2 * krad + 1, sizeof(float));
        buf    = (float *)calloc(nx, sizeof(float));
        if (kernel && buf) {
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
        if (kernel && buf) {
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
        if (kernel && buf) {
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
        }
        free(kernel); free(buf);
    }
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
    gaussian_blur_3d(om, nx, ny, nz, dx, dy, dz, sigma);
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
                            nalm[dd] = (int)(1.5f * nalm[dd]) + 16;
                            elm[dd] = (int *)realloc(elm[dd], sizeof(int) * nalm[dd]);
                            if (elm[dd] == NULL) { nalm[dd] = nelm[dd] = 0; continue; }
                        }
                        elm[dd][nelm[dd]++] = ss;
                        ntot++; nsaved++;
                    }
                }
            }
        }
        if (nsaved > 1) ndup++;
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
                elm[dd] = (int *)realloc(elm[dd], sizeof(int) * nelm[dd]);
            }
            nsav++;
        }
    }
    free(nalm);

    if (nsav == 0) {
        fprintf(stderr, "allineate: BLOK set has 0 surviving bloks\n");
        free(nelm); free(elm); return NULL;
    }

    /* Build output struct */
    gbs = (GA_BLOK_set *)malloc(sizeof(GA_BLOK_set));
    if (gbs == NULL) { free(nelm); free(elm); return NULL; }
    gbs->num  = nsav;
    gbs->nelm = (int *)  calloc(nsav, sizeof(int));
    gbs->elm  = (int **) calloc(nsav, sizeof(int *));
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
/*========== SECTION 7: 2D HISTOGRAM AND CORRELATION METRICS ==============*/
/*==========================================================================*/

/* Shannon entropy function */
#define SHANENT(z) (((z) <= 0.0f) ? 0.0f : (z)*logf(z))

/* Static global 2D histogram state */
static float *al_xc = NULL, *al_yc = NULL, *al_xyc = NULL;
static float al_nww = 0.0f;
static int al_nbin = 0, al_nbp = 0, al_nbm = 0;
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

static void clear_2Dhist(void)
{
    if (al_xc)  { free(al_xc);  al_xc = NULL; }
    if (al_yc)  { free(al_yc);  al_yc = NULL; }
    if (al_xyc) { free(al_xyc); al_xyc = NULL; }
    al_nbin = al_nbp = al_nbm = 0; al_nww = 0.0f;
}

/* Build 2D histogram with linear spread (LINHIST mode from AFNI) */
static void build_2Dhist(int n, float xbot, float xtop, float *x,
                                float ybot, float ytop, float *y, float *w)
{
    int ii, jj, kk, ngood;
    float xb, xi, yb, yi, xx, yy, x1, y1, ww;
    unsigned char *good;

    if (n <= 9 || x == NULL || y == NULL) return;

    clear_2Dhist();

    good = (unsigned char *)malloc(n);
    if (!good) return;
    for (ii = 0; ii < n; ii++)
        good[ii] = GOODVAL(x[ii]) && GOODVAL(y[ii]);

    /* Find range if not specified */
    if (xbot >= xtop) {
        xbot = WAY_BIG; xtop = -WAY_BIG;
        for (ii = 0; ii < n; ii++)
            if (good[ii]) {
                if (x[ii] > xtop) xtop = x[ii];
                else if (x[ii] < xbot) xbot = x[ii];
            }
    }
    if (xbot >= xtop) { free(good); return; }

    if (ybot >= ytop) {
        ybot = WAY_BIG; ytop = -WAY_BIG;
        for (ii = 0; ii < n; ii++)
            if (good[ii]) {
                if (y[ii] > ytop) ytop = y[ii];
                else if (y[ii] < ybot) ybot = y[ii];
            }
    }
    if (ybot >= ytop) { free(good); return; }

    /* Count good values in range */
    memset(good, 0, n);
    for (ngood = ii = 0; ii < n; ii++) {
        if (RANGVAL(x[ii], xbot, xtop) && RANGVAL(y[ii], ybot, ytop) && WW(ii) > 0.0f) {
            good[ii] = 1; ngood++;
        }
    }
    if (ngood == 0) { free(good); return; }

    /* Compute number of bins */
    al_nbin = (int)pow((double)n, al_hpow);
    if (al_nbin > 255) al_nbin = 255;
    else if (al_nbin < 3) al_nbin = 3;
    al_nbp = al_nbin + 1;
    al_nbm = al_nbin - 1;

    al_xc  = (float *)calloc(al_nbp, sizeof(float));
    al_yc  = (float *)calloc(al_nbp, sizeof(float));
    al_xyc = (float *)calloc(al_nbp * al_nbp, sizeof(float));
    if (!al_xc || !al_yc || !al_xyc) { clear_2Dhist(); free(good); return; }
    al_nww = 0.0f;

    /* Equal-size bins with linear spread */
    xb = xbot; xi = al_nbm / (xtop - xbot);
    yb = ybot; yi = al_nbm / (ytop - ybot);
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
    free(good);
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

/* Pearson correlation (nondestructive) */
static float al_pearson_corr(int n, float *x, float *y)
{
    float xv = 0.0f, yv = 0.0f, xy = 0.0f, vv, ww;
    float xm = 0.0f, ym = 0.0f;
    int ii;
    if (n < 2) return 0.0f;
    for (ii = 0; ii < n; ii++) { xm += x[ii]; ym += y[ii]; }
    xm /= n; ym /= n;
    for (ii = 0; ii < n; ii++) {
        vv = x[ii] - xm; ww = y[ii] - ym;
        xv += vv * vv; yv += ww * ww; xy += vv * ww;
    }
    if (xv <= 0.0f || yv <= 0.0f) return 0.0f;
    return xy / sqrtf(xv * yv);
}

/* Weighted Pearson */
static float al_pearson_corr_wt(int n, float *x, float *y, float *wt)
{
    float xv = 0.0f, yv = 0.0f, xy = 0.0f, vv, ww;
    float xm = 0.0f, ym = 0.0f, ws = 0.0f;
    int ii;
    if (wt == NULL) return al_pearson_corr(n, x, y);
    for (ii = 0; ii < n; ii++) {
        xm += wt[ii] * x[ii]; ym += wt[ii] * y[ii]; ws += wt[ii];
    }
    if (ws <= 0.0f) return 0.0f;
    xm /= ws; ym /= ws;
    for (ii = 0; ii < n; ii++) {
        vv = x[ii] - xm; ww = y[ii] - ym;
        xv += wt[ii] * vv * vv; yv += wt[ii] * ww * ww; xy += wt[ii] * vv * ww;
    }
    if (xv <= 0.0f || yv <= 0.0f) return 0.0f;
    return xy / sqrtf(xv * yv);
}

/* Combined Hellinger, MI, NMI, CRA from built histogram
   Returns: .a=hel, .b=mi, .c=nmi, .d=crA */
typedef struct { float a, b, c, d; } float_quad;

static float_quad al_helmicra(int n, float xbot, float xtop, float *x,
                              float ybot, float ytop, float *y, float *w)
{
    int ii, jj;
    float hel, pq, vv, uu;
    float cyvar, uyvar, yrat, xrat;
    float_quad hmc = {0.0f, 0.0f, 0.0f, 0.0f};

    build_2Dhist(n, xbot, xtop, x, ybot, ytop, y, w);
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
    if (npar >= 8) { b = parvec[7]; if (b <= 0.10f || b >= 10.0f) b = 1.0f; }
    if (npar >= 9) { c = parvec[8]; if (c <= 0.10f || c >= 10.0f) c = 1.0f; }
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
    static mat44 gam;
    int ii;

    if (npar > 0 && wpar != NULL)
        gam = GA_setup_affine(npar, wpar);

    if (npt <= 0 || xi == NULL || xo == NULL) return;

    for (ii = 0; ii < npt; ii++)
        mat44_vec(gam, xi[ii], yi[ii], zi[ii], &xo[ii], &yo[ii], &zo[ii]);
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

/*--- Get warped target values at control points (or all points) ---*/
static void GA_get_warped_values(int nmpar, double *mpar, float *avm)
{
    int npar, ii, jj, kk, qq, pp, npp, mm, nx, ny, nxy, npt, nall, nper, clip = 0;
    float *wpar, v;
    float *imf = NULL, *jmf = NULL, *kmf = NULL;
    float *imw, *jmw, *kmw;
    float *aim;

    npar = gstup->wfunc_numpar;
    wpar = (float *)calloc(npar, sizeof(float));
    if (!wpar) return;
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

    /* Space for control points */
    if (mpar == NULL || gstup->im_ar == NULL) {
        npt = gstup->bnx * gstup->bny * gstup->bnz;
        nall = (nper < npt) ? nper : npt;
        imf = (float *)calloc(nall, sizeof(float));
        jmf = (float *)calloc(nall, sizeof(float));
        kmf = (float *)calloc(nall, sizeof(float));
    } else {
        npt = gstup->npt_match;
        nall = (nper < npt) ? nper : npt;
    }

    imw = (float *)calloc(nall, sizeof(float));
    jmw = (float *)calloc(nall, sizeof(float));
    kmw = (float *)calloc(nall, sizeof(float));
    if (!imw || !jmw || !kmw || (imf && (!imf || !jmf || !kmf))) {
        free(imf); free(jmf); free(kmf); free(imw); free(jmw); free(kmw);
        free(wpar); return;
    }

    nx = gstup->bnx; ny = gstup->bny; nxy = nx * ny;

    /* Send parameters to warp function for setup */
    gstup->wfunc(npar, wpar, 0, NULL, NULL, NULL, NULL, NULL, NULL);

    /* Choose source image (smoothed if available) */
    aim = (gstup->ajims != NULL && mpar != NULL) ? gstup->ajims : gstup->ajim;

    /* Process in chunks */
    for (pp = 0; pp < npt; pp += nall) {
        npp = nall;
        if (npp > npt - pp) npp = npt - pp;

        if (mpar == NULL || gstup->im_ar == NULL) {
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

    free(kmw); free(jmw); free(imw);
    if (mpar == NULL || gstup->im_ar == NULL) {
        free(kmf); free(jmf); free(imf);
    }
    free(wpar);

    /* Clip interpolated values to source range */
    if (gstup->interp_code == AL_INTERP_CUBIC) {
        float bb = gstup->ajbot, tt = gstup->ajtop;
        for (pp = 0; pp < npt; pp++)
            if (avm[pp] < bb) avm[pp] = bb;
            else if (avm[pp] > tt) avm[pp] = tt;
    }
}

/*--- Local Pearson correlation using BLOKs ---*/
static float GA_pearson_local(int npt, float *avm, float *bvm, float *wvm)
{
    GA_BLOK_set *gbs;
    int nblok, nelm, *elm, dd, ii, jj;
    float xv, yv, xy, xm, ym, vv, ww, ws, wss, pcor, wt, psum = 0.0f, pabs;

    if (gstup->blokset == NULL) {
        float rad = gstup->blokrad, mrad;
        if (gstup->smooth_code > 0 && gstup->smooth_radius_base > 0.0f)
            rad = sqrtf(rad * rad + gstup->smooth_radius_base * gstup->smooth_radius_base);
        mrad = 1.2345f * (gstup->base_di + gstup->base_dj + gstup->base_dk);
        if (rad < mrad) rad = mrad;
        gstup->blokset = create_GA_BLOK_set(
            gstup->bnx, gstup->bny, gstup->bnz,
            gstup->base_di, gstup->base_dj, gstup->base_dk,
            gstup->npt_match, gstup->im_ar, gstup->jm_ar, gstup->km_ar,
            gstup->bloktype, rad, gstup->blokmin, 1.0f, 0);
        if (gstup->blokset == NULL) return AL_BIGVAL;
    }

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
        psum += ws * pcor * pabs;
        wss  += ws;
    }

    if (wss <= 0.0f) return AL_BIGVAL;
    return 0.25f * psum / wss;
}

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

    for (pp = qq = 0; pp < npt; pp++) {
        if (bsar[pp]) {
            ii = pp % nx; kk = pp / nxy; jj = (pp - kk * nxy) / nx;
            imf[qq] = (float)ii; jmf[qq] = (float)jj; kmf[qq] = (float)kk; qq++;
        }
    }

    imw = (float *)malloc(sizeof(float) * nqq);
    jmw = (float *)malloc(sizeof(float) * nqq);
    kmw = (float *)malloc(sizeof(float) * nqq);

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

/*--- Compute cost function ---*/
static double GA_scalar_costfun(int meth, int npt,
                                float *avm, float *bvm, float *wvm)
{
    double val = 0.0;

    switch (meth) {
        default:
        case GA_MATCH_PEARSON_SCALAR:
            val = (double)al_pearson_corr_wt(npt, avm, bvm, wvm);
            val = 1.0 - fabs(val);
            break;

        case GA_MATCH_PEARSON_LOCALS: /* LPC (signed) */
            val = (double)GA_pearson_local(npt, avm, bvm, wvm);
            break;

        case GA_MATCH_PEARSON_LOCALA: /* LPA */
            val = (double)GA_pearson_local(npt, avm, bvm, wvm);
            val = 1.0 - fabs(val);
            break;

        case GA_MATCH_LPC_MICHO_SCALAR: { /* lpc+ZZ */
            val = (double)GA_pearson_local(npt, avm, bvm, wvm);
            if (gstup->micho_hel != 0.0 || gstup->micho_mi  != 0.0 ||
                gstup->micho_nmi != 0.0 || gstup->micho_crA != 0.0) {
                float_quad hmc;
                float ovv;
                hmc = al_helmicra(npt, gstup->ajbot, gstup->ajclip, avm,
                                       gstup->bsbot, gstup->bsclip, bvm, wvm);
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
    }

    if (!isfinite(val)) val = (double)AL_BIGVAL;
    return val;
}

/*--- NEWUOA callback: evaluate cost from scaled parameters ---*/
static double GA_scalar_fitter(int npar, double *mpar)
{
    double val;
    float *avm, *bvm, *wvm;

    avm = (float *)calloc(gstup->npt_match, sizeof(float));
    if (!avm) return (double)AL_BIGVAL;

    GA_get_warped_values(npar, mpar, avm);

    bvm = gstup->bvm;
    wvm = gstup->wvm;

    if (gstup->need_hist_setup) {
        /* Just set flag to 0; histogram will be built inside costfun */
        gstup->need_hist_setup = 0;
    }

    val = GA_scalar_costfun(gstup->match_code, gstup->npt_match, avm, bvm, wvm);
    free(avm);
    return val;
}

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
static void al_scalar_setup(GA_setup *stup, int new_base, int new_targ)
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
    } else {
        if (stup->bsims) { free(stup->bsims); stup->bsims = NULL; }
    }

    /* Smooth source if needed */
    if (stup->smooth_code > 0 && stup->smooth_radius_targ > 0.0f) {
        if (stup->ajims) free(stup->ajims);
        stup->ajims = al_smooth(stup->ajim, stup->anx, stup->any, stup->anz,
                                stup->adx, stup->ady, stup->adz,
                                stup->smooth_radius_targ);
    } else {
        if (stup->ajims) { free(stup->ajims); stup->ajims = NULL; }
    }

    /* Get min/max of source image */
    {
        float *src = (stup->ajims != NULL) ? stup->ajims : stup->ajim;
        int nvox = stup->anx * stup->any * stup->anz;
        stup->ajbot = src[0]; stup->ajtop = src[0];
        for (ii = 1; ii < nvox; ii++) {
            if (src[ii] < stup->ajbot) stup->ajbot = src[ii];
            if (src[ii] > stup->ajtop) stup->ajtop = src[ii];
        }
        stup->ajclip = stup->ajtop;
    }

    /* Get min/max of base image */
    {
        float *bas = (stup->bsims != NULL) ? stup->bsims : stup->bsim;
        int nvox = nx * ny * nz;
        stup->bsbot = bas[0]; stup->bstop = bas[0];
        for (ii = 1; ii < nvox; ii++) {
            if (bas[ii] < stup->bsbot) stup->bsbot = bas[ii];
            if (bas[ii] > stup->bstop) stup->bstop = bas[ii];
        }
        stup->bsclip = stup->bstop;
    }

    /* Determine number of matching points */
    nmatch = stup->npt_match;
    if (nmatch <= 9 || nmatch > nx * ny * nz)
        nmatch = nx * ny * nz;
    if (stup->nmask > 0 && nmatch > stup->nmask)
        nmatch = stup->nmask;
    stup->npt_match = nmatch;

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
        if (!stup->im_ar || !stup->jm_ar || !stup->km_ar) return;
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
        }
    } else {
        /* Subset of points */
        int nvox = nx * ny * nz;
        int dm = ga_find_relprime_fixed(nvox);
        stup->im_ar = (float *)malloc(sizeof(float) * nmatch);
        stup->jm_ar = (float *)malloc(sizeof(float) * nmatch);
        stup->km_ar = (float *)malloc(sizeof(float) * nmatch);
        if (!stup->im_ar || !stup->jm_ar || !stup->km_ar) return;
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

    stup->need_hist_setup = 1;
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
    stup->vbest = (float)GA_scalar_fitter(stup->wfunc_numfree, wpar);

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

/*--- Random startup search (coarse pass) ---*/
static void al_scalar_ransetup(GA_setup *stup, int nrand)
{
#define NKEEP (3*PARAM_MAXTRIAL+1)
    double *wpar, *spar, val, vbest, *bpar;
    double *kpar[NKEEP], kval[NKEEP];
    int ii, qq, ss, nfr, kk, jj, ngrid, ngtot, ngood, twof, maxstep;
    int ival[NKEEP], rval[NKEEP];
    float fval[NKEEP];

    if (stup == NULL || stup->setup != AL_SMAGIC) return;
    if (nrand < NKEEP) nrand = NKEEP + 13;

    GA_param_setup(stup); gstup = stup;
    if (stup->wfunc_numfree <= 0) return;

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
    int icod = stup->interp_code;
    stup->interp_code = AL_INTERP_LINEAR;

    wpar = (double *)calloc(nfr, sizeof(double));
    spar = (double *)calloc(nfr, sizeof(double));
    for (kk = 0; kk < NKEEP; kk++)
        kpar[kk] = (double *)calloc(nfr, sizeof(double));

    /* Evaluate center of parameter space */
    for (qq = 0; qq < nfr; qq++) wpar[qq] = 0.5;
    val = GA_scalar_fitter(nfr, wpar);
    memcpy(kpar[0], wpar, sizeof(double) * nfr);
    kval[0] = val; rval[0] = 0;
    for (kk = 1; kk < NKEEP; kk++) { rval[kk] = 0; kval[kk] = AL_BIGVAL; }

    twof = 1 << nfr;
    myunif_reset(3456789);

    fprintf(stderr, " + Coarse search: %d grid + %d random trials\n", ngtot, nrand);

    for (ii = 0; ii < nrand + ngtot; ii++) {
        if (ii < ngtot) {
            double kp;
            val = 0.5 / (ngrid + 1.0); ss = ii;
            for (qq = 0; qq < nfr; qq++) {
                kk = ss % ngrid; ss = ss / ngrid;
                kp = (kk == 0) ? 0.5 : (kk + 1.0);
                wpar[qq] = 0.5 + kp * val;
            }
        } else {
            for (qq = 0; qq < nfr; qq++) wpar[qq] = 0.5 * (1.05 + 0.90 * myunif());
        }

        for (ss = 0; ss < twof; ss++) {
            for (qq = 0; qq < nfr; qq++)
                spar[qq] = (ss & (1 << qq)) ? 1.0 - wpar[qq] : wpar[qq];
            val = GA_scalar_fitter(nfr, spar);
            for (kk = 0; kk < NKEEP; kk++) {
                if (val < kval[kk]) {
                    for (jj = NKEEP - 2; jj >= kk; jj--) {
                        memcpy(kpar[jj + 1], kpar[jj], sizeof(double) * nfr);
                        kval[jj + 1] = kval[jj]; rval[jj + 1] = rval[jj];
                    }
                    memcpy(kpar[kk], spar, sizeof(double) * nfr);
                    kval[kk] = val; rval[kk] = (ii >= ngtot);
                    break;
                }
            }
        }
    }

    for (ngood = kk = 0; kk < NKEEP && kval[kk] < AL_BIGVAL; kk++, ngood++) ;
    if (ngood < 1) {
        fprintf(stderr, "allineate: no good starting locations found\n");
        free(wpar); free(spar);
        for (kk = 0; kk < NKEEP; kk++) free(kpar[kk]);
        stup->interp_code = icod;
        return;
    }

    /* Make sure all in 0..1 range */
    for (kk = 0; kk < ngood; kk++)
        for (ii = 0; ii < nfr; ii++) kpar[kk][ii] = PRED01(kpar[kk][ii]);

    /* Little optimization on each set */
    vbest = AL_BIGVAL; jj = 0;
    stup->interp_code = AL_INTERP_LINEAR;
    maxstep = 11 * nfr + 17;
    fprintf(stderr, " + Refining %d candidate parameter sets\n", ngood);
    for (kk = 0; kk < ngood; kk++) {
        if (kval[kk] >= AL_BIGVAL) continue;
        powell_newuoa(nfr, kpar[kk], 0.05, 0.001, maxstep, GA_scalar_fitter);
        kval[kk] = GA_scalar_fitter(nfr, kpar[kk]);
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

    free(wpar); free(spar);
    for (kk = 0; kk < NKEEP; kk++) free(kpar[kk]);
    stup->interp_code = icod;
#undef NKEEP
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

    /* Send parameters to warp function */
    wfunc(npar, wpar, 0, NULL, NULL, NULL, NULL, NULL, NULL);

    nx = nnx; ny = nny; nz = nnz; nxy = nx * ny; npt = nxy * nz;
    war = (float *)calloc(npt, sizeof(float));
    if (!war) return NULL;

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
    if (icode == AL_INTERP_CUBIC) {
        int nvox_t = tnx * tny * tnz;
        float bb = imtarg[0], tt = imtarg[0];
        for (pp = 1; pp < nvox_t; pp++) {
            if (imtarg[pp] < bb) bb = imtarg[pp];
            if (imtarg[pp] > tt) tt = imtarg[pp];
        }
        for (pp = 0; pp < npt; pp++) {
            if (war[pp] < bb) war[pp] = bb;
            else if (war[pp] > tt) war[pp] = tt;
        }
    }

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

    /* Squash large values */
    mx = 0.0f;
    for (ii = 0; ii < nxyz; ii++) if (wf[ii] > mx) mx = wf[ii];
    clip = 0.5f * mx;
    for (ii = 0; ii < nxyz; ii++) if (wf[ii] > clip) wf[ii] = clip;

    /* Gaussian blur */
    sigma = 2.25f * (dx + dy + dz) / 3.0f;
    gaussian_blur_3d(wf, nx, ny, nz, dx, dy, dz, sigma);

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

/* Compute simple source automask (binary) */
static unsigned char *al_automask(float *srcim, int nvox)
{
    unsigned char *mask;
    float *sorted, pct98;
    int ii, np;

    mask = (unsigned char *)calloc(nvox, sizeof(unsigned char));
    if (!mask) return NULL;

    sorted = (float *)malloc(sizeof(float) * nvox);
    if (!sorted) { free(mask); return NULL; }
    memcpy(sorted, srcim, sizeof(float) * nvox);
    qsort_floatint(nvox, sorted, NULL); /* sort ascending */

    np = (int)(0.98f * nvox);
    if (np >= nvox) np = nvox - 1;
    pct98 = sorted[np];
    free(sorted);

    float thresh = 0.10f * pct98;
    for (ii = 0; ii < nvox; ii++)
        mask[ii] = (srcim[ii] > thresh) ? 1 : 0;

    return mask;
}

/*==========================================================================*/
/*============= SECTION 10: BEFAFTER COORDINATE SETUP =====================*/
/*==========================================================================*/

static void al_setup_befafter(mat44 base_cmat, mat44 base_imat,
                              mat44 targ_cmat, mat44 targ_imat)
{
    /* before = targ_imat * (transform from base_xyz to targ_ijk)
       after  = base_imat converted to base_ijk from base_xyz

       The affine parameters operate in xyz (DICOM) space.
       before: base_ijk -> base_xyz (via base_cmat)
       after:  targ_xyz -> targ_ijk (via targ_imat)

       So the full chain is:
       base_ijk -> base_xyz [before] -> targ_xyz [affine] -> targ_ijk [after]

       In AFNI's convention:
       aff_before = base_cmat (base index -> base xyz)
       aff_after  = targ_imat (targ xyz -> targ index)
    */

    aff_before = base_cmat;
    aff_after  = targ_imat;
    aff_use_before = 1;
    aff_use_after  = 1;
}

/*==========================================================================*/
/*================ SECTION 11: PUBLIC FUNCTION nii_allineate() ============*/
/*==========================================================================*/

/* Extract float data from NIfTI image. Returns malloc'd float array. */
static float *nii_to_float(nifti_image *nim)
{
    int nvox, ii;
    float *fdata;
    if (nim == NULL || nim->data == NULL) return NULL;
    nvox = (int)(nim->nx * nim->ny * nim->nz);
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
        default:
            /* Best effort: just zero-fill */
            break;
    }

    /* Apply scl_slope/scl_inter if not already handled and if set */
    if (nim->datatype != DT_INT16 && nim->scl_slope != 0.0f &&
        !(nim->scl_slope == 1.0f && nim->scl_inter == 0.0f)) {
        for (ii = 0; ii < nvox; ii++)
            fdata[ii] = fdata[ii] * nim->scl_slope + nim->scl_inter;
    }

    return fdata;
}

int nii_allineate(nifti_image *source, nifti_image *base)
{
    GA_setup stup;
    GA_param params[12];
    float *bsim, *ajim;   /* base and source float data */
    float *wght;           /* autoweight */
    unsigned char *smask;  /* source automask */
    int bnx, bny, bnz, anx, any, anz;
    float bdx, bdy, bdz, adx, ady, adz;
    mat44 base_cmat, base_imat, targ_cmat, targ_imat;
    int jj, ii, nvox_base, nvox_src;
    int nfunc;

    if (source == NULL || base == NULL) {
        fprintf(stderr, "allineate: NULL input image\n");
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

    nvox_base = bnx * bny * bnz;
    nvox_src  = anx * any * anz;

    /* --- 2. Get sform matrices --- */
    if (base->sform_code > 0)
        base_cmat = dmat44_to_mat44(base->sto_xyz);
    else {
        base_cmat = mat44_diag(bdx, bdy, bdz);
        base_cmat.m[0][3] = -(bnx - 1) * 0.5f * bdx;
        base_cmat.m[1][3] = -(bny - 1) * 0.5f * bdy;
        base_cmat.m[2][3] = -(bnz - 1) * 0.5f * bdz;
    }
    base_imat = nifti_mat44_inverse(base_cmat);

    if (source->sform_code > 0)
        targ_cmat = dmat44_to_mat44(source->sto_xyz);
    else {
        targ_cmat = mat44_diag(adx, ady, adz);
        targ_cmat.m[0][3] = -(anx - 1) * 0.5f * adx;
        targ_cmat.m[1][3] = -(any - 1) * 0.5f * ady;
        targ_cmat.m[2][3] = -(anz - 1) * 0.5f * adz;
    }
    targ_imat = nifti_mat44_inverse(targ_cmat);

    /* --- 3. Compute autoweight --- */
    fprintf(stderr, " + Computing autoweight from base image\n");
    wght = al_autoweight(bsim, bnx, bny, bnz, bdx, bdy, bdz);

    /* --- 4. Compute source automask --- */
    fprintf(stderr, " + Computing source automask\n");
    smask = al_automask(ajim, nvox_src);

    /* --- 5. Set up GA_setup structure --- */
    memset(&stup, 0, sizeof(GA_setup));

    stup.bsim = bsim;     stup.bsims = NULL;
    stup.bnx = bnx;       stup.bny = bny;       stup.bnz = bnz;
    stup.bdx = bdx;       stup.bdy = bdy;       stup.bdz = bdz;
    stup.ajim = ajim;     stup.ajims = NULL;     stup.ajimor = NULL;
    stup.anx = anx;       stup.any = any;        stup.anz = anz;
    stup.adx = adx;       stup.ady = ady;        stup.adz = adz;
    stup.ajmask = smask;
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
        stup.nmask = 0;
        float wmx = 0.0f;
        for (ii = 0; ii < nvox_base; ii++) if (wght[ii] > wmx) wmx = wght[ii];
        if (wmx > 0.0f) {
            float inv = 1.0f / wmx;
            for (ii = 0; ii < nvox_base; ii++) {
                wght[ii] = fabsf(wght[ii]) * inv;
                stup.bmask[ii] = (wght[ii] > 0.0f) ? 1 : 0;
                if (stup.bmask[ii]) stup.nmask++;
            }
        }
    }

    /* --- 6. Configure lpc+ZZ --- */
    stup.match_code = GA_MATCH_LPC_MICHO_SCALAR;
    stup.micho_mi  = 0.2;
    stup.micho_nmi = 0.2;
    stup.micho_crA = 0.4;
    stup.micho_hel = 0.4;
    stup.micho_ov  = 0.4;
    stup.micho_zfinal = 1; /* +ZZ */

    /* --- 7. Set up affine 12-param warp --- */
    stup.wfunc = al_wfunc_affine;
    stup.wfunc_numpar = 12;
    stup.wfunc_param = params;
    memset(params, 0, sizeof(params));

    /* Set up befafter matrices */
    al_setup_befafter(base_cmat, base_imat, targ_cmat, targ_imat);

    /* Compute shift ranges (about 1/3 of base FOV in each direction) */
    float xxx = 0.321f * (bnx - 1);
    float yyy = 0.321f * (bny - 1);
    float zzz = 0.321f * (bnz - 1);
    float sm, sp;

    /* Transform to mm space for shift range */
    float xxx_m = 0.01f, yyy_m = 0.01f, zzz_m = 0.01f;
    for (ii = -1; ii <= 1; ii += 2)
        for (jj = -1; jj <= 1; jj += 2)
            for (int kk2 = -1; kk2 <= 1; kk2 += 2) {
                float xp, yp, zp;
                xp = base_cmat.m[0][0]*(ii*xxx) + base_cmat.m[0][1]*(jj*yyy) + base_cmat.m[0][2]*(kk2*zzz);
                yp = base_cmat.m[1][0]*(ii*xxx) + base_cmat.m[1][1]*(jj*yyy) + base_cmat.m[1][2]*(kk2*zzz);
                zp = base_cmat.m[2][0]*(ii*xxx) + base_cmat.m[2][1]*(jj*yyy) + base_cmat.m[2][2]*(kk2*zzz);
                xp = fabsf(xp); yp = fabsf(yp); zp = fabsf(zp);
                if (xp > xxx_m) xxx_m = xp;
                if (yp > yyy_m) yyy_m = yp;
                if (zp > zzz_m) zzz_m = zp;
            }

    /* DEFPAR macro equivalent */
#define SETPAR(p,nm,bb,tt,id) do {          \
    params[p].min = (bb); params[p].max = (tt);    \
    params[p].ident = (id); params[p].val_init = (id); \
    params[p].val_pinit = (id); params[p].val_fixed = (id); \
    params[p].val_out = (id); params[p].fixed = 0;   \
    strncpy(params[p].name, (nm), 31);           \
} while(0)

    SETPAR(0, "x-shift", -xxx_m, xxx_m, 0.0f);
    SETPAR(1, "y-shift", -yyy_m, yyy_m, 0.0f);
    SETPAR(2, "z-shift", -zzz_m, zzz_m, 0.0f);

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

    /* --- 8. BLOK set params --- */
    stup.bloktype = GA_BLOK_RHDD;
    stup.blokrad  = 6.54321f * (stup.base_di + stup.base_dj + stup.base_dk) / 3.0f;
    stup.blokmin  = 0; /* auto */
    stup.blokset  = NULL;

    /* ========== COARSE PASS (twopass) ========== */
    fprintf(stderr, " + *** Coarse pass begins ***\n");

    int ntask = nvox_base;
    stup.interp_code = AL_INTERP_LINEAR;
    stup.npt_match   = ntask / 15;
    if (stup.npt_match < 9999) stup.npt_match = 9999;
    stup.smooth_code = GA_SMOOTH_GAUSSIAN;
    stup.smooth_radius_base = stup.smooth_radius_targ = 7.777f;

    al_scalar_setup(&stup, 1, 1);

    /* Temporarily freeze params beyond first 6 for coarse search */
    int nparam_free = 12;
    int nptwo = 6;
    if (nparam_free > nptwo) {
        for (ii = jj = 0; jj < stup.wfunc_numpar; jj++) {
            if (!stup.wfunc_param[jj].fixed) {
                ii++;
                if (ii > nptwo) stup.wfunc_param[jj].fixed = 1;
            }
        }
    }

    powell_set_mfac(1.0f, 3.0f);

    int tbest = 5;
    int nrand = 17 + 4 * tbest + 2 * tbest; /* extra for BLOK methods */
    if (nrand < 31) nrand = 31;
    al_scalar_ransetup(&stup, nrand);

    /* Unfreeze temporarily frozen params */
    for (jj = 0; jj < stup.wfunc_numpar; jj++)
        if (stup.wfunc_param[jj].fixed == 1) stup.wfunc_param[jj].fixed = 0;

    /* Store trial parameter sets */
    int tb = tbest;
    if (tb > stup.wfunc_ntrial) tb = stup.wfunc_ntrial;
    float tfparm[PARAM_MAXTRIAL + 2][12];
    float tfcost[PARAM_MAXTRIAL + 2];
    int tfdone;

    for (int ib = 0; ib < tb; ib++)
        for (jj = 0; jj < 12; jj++)
            tfparm[ib][jj] = stup.wfunc_param[jj].val_trial[ib];

    /* Add identity transform */
    for (jj = 0; jj < 12; jj++)
        tfparm[tb][jj] = stup.wfunc_param[jj].val_pinit;
    tfdone = tb + 1;

    /* Refinement rounds */
    double rad = 0.0555;
    int nrefine = 3;
    for (int rr = 0; rr < nrefine; rr++, rad *= 0.6789) {
        fprintf(stderr, " + Refinement #%d on %d parameter sets\n", rr + 1, tfdone);
        powell_set_mfac(1.0f, 5.0f + 2.0f * rr);

        stup.smooth_radius_base *= 0.7777f;
        stup.smooth_radius_targ *= 0.7777f;
        stup.npt_match = (int)(stup.npt_match * 1.5);

        al_scalar_setup(&stup, 0, 0);

        for (int ib = 0; ib < tfdone; ib++) {
            for (jj = 0; jj < stup.wfunc_numpar; jj++)
                stup.wfunc_param[jj].val_init = tfparm[ib][jj];
            nfunc = al_scalar_optim(&stup, rad, 0.01 * rad, 99 + 11 * rr);
            for (jj = 0; jj < stup.wfunc_numpar; jj++)
                tfparm[ib][jj] = stup.wfunc_param[jj].val_out;
            tfcost[ib] = stup.vbest;
        }
    }

    /* ========== FINE PASS ========== */
    fprintf(stderr, " + *** Fine pass begins ***\n");

    stup.interp_code = AL_INTERP_CUBIC;
    stup.smooth_code = 0;
    stup.smooth_radius_base = 0.0f;
    stup.smooth_radius_targ = 0.0f;
    stup.npt_match = nvox_base;

    al_scalar_setup(&stup, 0, 0);
    powell_set_mfac(0.0f, 0.0f);

    /* Pick best from coarse */
    {
        int kb = 0;
        float cbest = 1.e+33f;
        GA_param_setup(&stup);
        gstup = &stup;

        for (int ib = 0; ib < tfdone; ib++) {
            for (jj = 0; jj < stup.wfunc_numpar; jj++)
                stup.wfunc_param[jj].val_init = tfparm[ib][jj];
            /* Evaluate cost */
            double *wpar = (double *)calloc(stup.wfunc_numfree, sizeof(double));
            for (ii = jj = 0; jj < stup.wfunc_numpar; jj++) {
                if (!stup.wfunc_param[jj].fixed) {
                    wpar[ii] = (stup.wfunc_param[jj].val_init - stup.wfunc_param[jj].min)
                               / stup.wfunc_param[jj].siz;
                    if (wpar[ii] < 0.0 || wpar[ii] > 1.0) wpar[ii] = PRED01(wpar[ii]);
                    ii++;
                }
            }
            double cost = GA_scalar_fitter(stup.wfunc_numfree, wpar);
            free(wpar);
            if (cost < cbest) { cbest = (float)cost; kb = ib; }
        }
        for (jj = 0; jj < stup.wfunc_numpar; jj++)
            stup.wfunc_param[jj].val_init = tfparm[kb][jj];
        fprintf(stderr, " + Best coarse cost = %f (set #%d)\n", cbest, kb + 1);
    }

    /* Final optimization */
    rad = 0.0444;
    nfunc = al_scalar_optim(&stup, rad, 0.001, 6666);
    fprintf(stderr, " + Fine cost = %f (%d evaluations)\n", stup.vbest, nfunc);

    /* --- 11. +ZZ refinal: zero micho weights, re-optimize with pure lpc --- */
    for (jj = 0; jj < stup.wfunc_numpar; jj++)
        stup.wfunc_param[jj].val_init = stup.wfunc_param[jj].val_out;
    stup.need_hist_setup = 1;

    /* Zero micho weights for final pure-lpc pass */
    double save_mi  = stup.micho_mi;
    double save_nmi = stup.micho_nmi;
    double save_crA = stup.micho_crA;
    double save_hel = stup.micho_hel;
    double save_ov  = stup.micho_ov;
    stup.micho_mi = stup.micho_nmi = stup.micho_crA = stup.micho_hel = stup.micho_ov = 0.0;

    fprintf(stderr, " + +ZZ refinal (pure lpc)\n");
    rad = 0.0666;
    powell_set_mfac(4.0f, 4.0f);
    nfunc = al_scalar_optim(&stup, rad, 0.001, 6666);
    fprintf(stderr, " + Refinal cost = %f (%d evaluations)\n", stup.vbest, nfunc);

    /* Restore micho weights */
    stup.micho_mi  = save_mi;
    stup.micho_nmi = save_nmi;
    stup.micho_crA = save_crA;
    stup.micho_hel = save_hel;
    stup.micho_ov  = save_ov;

    /* --- 12. Apply final warp --- */
    fprintf(stderr, " + Applying final warp with cubic interpolation\n");
    {
        float wpar[12];
        for (jj = 0; jj < 12; jj++)
            wpar[jj] = stup.wfunc_param[jj].val_out;

        fprintf(stderr, " + Final parameters:");
        for (jj = 0; jj < 12; jj++) fprintf(stderr, " %.4f", wpar[jj]);
        fprintf(stderr, "\n");

        float *warped = al_scalar_warpone(12, wpar, al_wfunc_affine,
                                          ajim, anx, any, anz,
                                          bnx, bny, bnz, AL_INTERP_CUBIC);
        if (warped == NULL) {
            fprintf(stderr, "allineate: failed to warp image\n");
            free(bsim); free(ajim);
            if (wght) free(wght);
            if (smask) free(smask);
            if (stup.bmask) free(stup.bmask);
            if (stup.bsims) free(stup.bsims);
            if (stup.ajims) free(stup.ajims);
            if (stup.bvm) free(stup.bvm);
            if (stup.wvm) free(stup.wvm);
            if (stup.im_ar) free(stup.im_ar);
            if (stup.jm_ar) free(stup.jm_ar);
            if (stup.km_ar) free(stup.km_ar);
            if (stup.blokset) free_GA_BLOK_set(stup.blokset);
            clear_2Dhist();
            return 1;
        }

        /* --- 13. Replace source->data with warped result, update dims/transforms --- */
        free(source->data);
        source->data = warped;
        source->datatype = DT_FLOAT32;
        source->nbyper = sizeof(float);
        source->scl_slope = 0.0f;
        source->scl_inter = 0.0f;
        source->nx = bnx; source->ny = bny; source->nz = bnz;
        source->dx = base->dx; source->dy = base->dy; source->dz = base->dz;
        source->dim[1] = bnx; source->dim[2] = bny; source->dim[3] = bnz;
        source->pixdim[1] = base->pixdim[1];
        source->pixdim[2] = base->pixdim[2];
        source->pixdim[3] = base->pixdim[3];
        source->nvox = (size_t)bnx * bny * bnz;
        source->nbyper = sizeof(float);
        source->datatype = DT_FLOAT32;
        /* Copy sform from base */
        source->sform_code = base->sform_code;
        source->sto_xyz = base->sto_xyz;
        source->sto_ijk = base->sto_ijk;
        source->qform_code = base->qform_code;
        source->qto_xyz = base->qto_xyz;
        source->qto_ijk = base->qto_ijk;
        source->quatern_b = base->quatern_b;
        source->quatern_c = base->quatern_c;
        source->quatern_d = base->quatern_d;
        source->qoffset_x = base->qoffset_x;
        source->qoffset_y = base->qoffset_y;
        source->qoffset_z = base->qoffset_z;
    }

    /* --- 14. Cleanup --- */
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
    clear_2Dhist();

    fprintf(stderr, " + Registration complete\n");
    return 0;
}
