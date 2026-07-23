/*----------------------------------------------------------------------------
 * qwarp.c — minimal nonlinear (3dQwarp) registration stage.
 *
 * ATTRIBUTED PUBLIC-DOMAIN PORT of AFNI's nonlinear-warp code (public domain,
 * NIMH/NIH; Robert W. Cox / "Zhark"): 3dQwarp.c, mri_nwarp.c, thd_incorrelate.c,
 * thd_cliplevel.c, thd_automask.c, edt_volpad.c, mri_genalign_util.c, pinned at AFNI
 * rev 506e48403 (AFNI_26.1.00-6). Optimizer: Powell NEWUOA (powell_newuoa.c).
 *
 * LICENSING BOUNDARY: AFNI's LICENSE.txt is a NIH Public Domain Notice with EXCEPTIONS
 * for bundled third-party code. Three utilities on this path are GPL-v2 / Medical College
 * of Wisconsin exceptions, NOT public domain — edt_blur.c (FIR Gaussian blur),
 * edt_buildmask.c (sphere mask), cs_qmed.c (median). These are NOT ported; instead their
 * function is provided by niimath's own BSD/public-domain coreFLT.c (Gaussian blur via
 * nifti_smooth_gauss_f32) or by small clean-room routines here. Only genuinely public-
 * domain AFNI files are ported. Never copied from any GPL AFNI file or from
 * niimath/src/GPL/. See qwarp.h and AGENTS.md.
 *
 * Ports the single operation:
 *     3dQwarp -blur 0 3 -source <moving> -base <stationary> -prefix <output>
 *
 * PIPELINE (qwarp_run): float extract -> data-dependent zeropad -> FWHM-3 source blur ->
 * mri_weightize of the padded base -> warpomatic (level-0 basis escalation SINCC/CUBIC/
 * QUINTIC + shrinking-patch schedule + hexahedron-energy penalty at level >= 3, clipped-
 * Pearson similarity, NEWUOA patch optimizer) -> single WSINC5 reslice of the unblurred
 * padded source -> crop back to the stationary grid. Equivalent (not bit-identical) to AFNI:
 * Pearson 0.98-0.996 vs 3dQwarp on 1mm MNI-space brains. Serial + OpenMP; the internal
 * IndexWarp3D and all state are private to this file. See AGENTS.md.
 *--------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <limits.h>
#include "qwarp.h"
#include "core32.h"   /* nifti_smooth_gauss_f32 — niimath's BSD/public-domain Gaussian */
#include "al_size_guard.h"

/* FWHM (voxels) -> Gaussian sigma (voxels): AFNI FWHM_TO_SIGMA (editvol.h). */
#define QW_FWHM_TO_SIGMA(f) (0.42466090 * (f))
#define QW_SRC_BLUR_FWHM    3.0    /* `-blur 0 3`: source blurred FWHM 3 voxels, base unblurred */

/*--- finite / usable-geometry helpers -------------------------------------*/

static int qw_finite(double v) { return (v == v) && (v > -1e300) && (v < 1e300); }

/* Reject NaN/Inf without isfinite() (unusable under -ffast-math): v==v kills NaN, the magnitude
 * bound kills +/-Inf (Inf compares outside +/-FLT_MAX ~ 3.4e38). Returns 1 if all finite. */
static int qw_all_finite(const float *f, int64_t n) {
    for (int64_t i = 0; i < n; i++) { float v = f[i]; if (!((v == v) && v >= -3.0e38f && v <= 3.0e38f)) return 0; }
    return 1;
}

/* Copy a nifti_dmat44 into M if it is finite and non-singular; return 0 if usable. */
static int qw_try_xform(const nifti_dmat44 *s, double M[4][4]) {
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
            if (!qw_finite(s->m[i][j])) return 1;
            M[i][j] = s->m[i][j];
        }
    double det =
        M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
      - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
      + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);
    if (!qw_finite(det) || fabs(det) < 1e-12) return 1;   /* singular/degenerate 3x3 */
    return 0;
}

/* Pick the usable voxel->world transform using this project's al_image_xform() policy: prefer
 * the coded form with the HIGHER code (sform on a tie), then fall back to the other coded form
 * if the preferred one is unusable (non-finite/singular). Unlike al_image_xform_or_pixdim(),
 * qwarp does NOT synthesize a pixdim frame — decision 3 rejects unusable geometry outright.
 *
 * This intentionally DUPLICATES allineate.c's al_image_xform() rather than calling it: qwarp.c is
 * designed to be independently linkable and promotable to niimath as a standalone module, so it does not depend on
 * allineate.c internals (see AGENTS.md). Keep the two in sync
 * if the selection policy ever changes; the qwarp C-API test locks the two-form fallback. */
static int qw_get_xform(const nifti_image *n, double M[4][4]) {
    const nifti_dmat44 *pref = NULL, *alt = NULL;
    if (n->sform_code > 0 && n->sform_code >= n->qform_code) {          /* sform preferred */
        pref = &n->sto_xyz; if (n->qform_code > 0) alt = &n->qto_xyz;
    } else if (n->qform_code > 0) {                                     /* qform preferred */
        pref = &n->qto_xyz; if (n->sform_code > 0) alt = &n->sto_xyz;
    } else return 1;                                                    /* neither coded */
    if (qw_try_xform(pref, M) == 0) return 0;
    if (alt && qw_try_xform(alt, M) == 0) return 0;
    return 1;
}

static void qw_apply(const double M[4][4], double x, double y, double z,
                     double *ox, double *oy, double *oz) {
    *ox = M[0][0]*x + M[0][1]*y + M[0][2]*z + M[0][3];
    *oy = M[1][0]*x + M[1][1]*y + M[1][2]*z + M[1][3];
    *oz = M[2][0]*x + M[2][1]*y + M[2][2]*z + M[2][3];
}

/* Single-volume 3D check: nvox must equal nx*ny*nz and every higher dim <= 1. */
static int qw_is_single_vol_3d(const nifti_image *n) {
    if (n->nx < 1 || n->ny < 1 || n->nz < 1) return 0;
    if (n->nt > 1 || n->nu > 1 || n->nv > 1 || n->nw > 1) return 0;
    int64_t n3 = n->nx * n->ny * n->nz;
    return n->nvox == n3;
}

/* Validate that moving and stationary share a grid (decision 3). Returns 0 if OK,
 * else prints a concise reason and returns nonzero. */
static int qw_grid_compat(const nifti_image *mov, const nifti_image *sta) {
    if (!mov || !sta || !mov->data || !sta->data) {
        fprintf(stderr, "qwarp: null image or image data\n"); return 1; }
    /* Bound each axis of BOTH inputs FIRST — before any nx*ny*nz product is formed. nifti dims are
     * int64; an oversize dimension would overflow the signed product (UB) in qw_is_single_vol_3d
     * and later truncate under the engine's `int` casts. 20000 >> any real brain grid; with dims
     * <= 20000 the int64 product cannot overflow. (Lower bound 5: the warpomatic level-0 bounds
     * assume a >=1-voxel interior and the optimizer's smallest patch is 5, so 1/2-voxel dims — incl.
     * a 2D nz==1 slab, which is unsupported/untested — would make reversed/degenerate patches.) */
    const nifti_image *im2[2] = { mov, sta };
    for (int t = 0; t < 2; t++) {
        const nifti_image *n = im2[t];
        if (n->nx > 20000 || n->ny > 20000 || n->nz > 20000) {
            fprintf(stderr, "qwarp: dimension too large for the int-based engine\n"); return 1; }
        if (n->nx < 5 || n->ny < 5 || n->nz < 5) {
            fprintf(stderr, "qwarp: grid too small / not 3D (each dimension must be >= 5 voxels)\n"); return 1; }
    }
    if (!qw_is_single_vol_3d(mov) || !qw_is_single_vol_3d(sta)) {   /* product now overflow-safe */
        fprintf(stderr, "qwarp: inputs must be single-volume 3D (no 4D/multi-volume)\n"); return 1; }
    /* The unpadded product must fit both the int-based engine and an addressable
       float buffer (padding only grows it; the padded product is re-checked). */
    int64_t unv = (int64_t)mov->nx * mov->ny * mov->nz;
    if (unv > (int64_t)INT_MAX || !al_float_nvox_fits((uint64_t)unv)) {
        fprintf(stderr, "qwarp: voxel count exceeds the supported maximum\n"); return 1; }
    if (mov->nx != sta->nx || mov->ny != sta->ny || mov->nz != sta->nz) {
        fprintf(stderr, "qwarp: moving %ldx%ldx%ld and stationary %ldx%ldx%ld dims differ "
                        "(inputs must already share the grid)\n",
                (long)mov->nx,(long)mov->ny,(long)mov->nz,
                (long)sta->nx,(long)sta->ny,(long)sta->nz); return 1; }
    if (!(qw_finite(mov->dx) && mov->dx > 0 && qw_finite(mov->dy) && mov->dy > 0 &&
          qw_finite(mov->dz) && mov->dz > 0 &&
          qw_finite(sta->dx) && sta->dx > 0 && qw_finite(sta->dy) && sta->dy > 0 &&
          qw_finite(sta->dz) && sta->dz > 0)) {
        fprintf(stderr, "qwarp: non-finite or non-positive voxel size\n"); return 1; }
    double Mm[4][4], Ms[4][4];
    if (qw_get_xform(mov, Mm) || qw_get_xform(sta, Ms)) {
        fprintf(stderr, "qwarp: unusable NIfTI geometry (no coded/usable sform or qform)\n"); return 1; }
    /* The two frames must agree: map the 8 index corners through each and require
     * every world position within a small fraction of the smallest voxel edge. */
    double vmin = sta->dx; if (sta->dy < vmin) vmin = sta->dy; if (sta->dz < vmin) vmin = sta->dz;
    double tol = 0.05 * vmin; if (!(tol > 0.0)) tol = 0.05;
    for (int c = 0; c < 8; c++) {
        double ix = (c&1)?mov->nx-1:0, iy = (c&2)?mov->ny-1:0, iz = (c&4)?mov->nz-1:0;
        double mx,my,mz, sx,sy,sz;
        qw_apply(Mm, ix,iy,iz, &mx,&my,&mz);
        qw_apply(Ms, ix,iy,iz, &sx,&sy,&sz);
        double dx=mx-sx, dy=my-sy, dz=mz-sz;
        if (sqrt(dx*dx+dy*dy+dz*dz) > tol) {
            fprintf(stderr, "qwarp: moving and stationary voxel->world transforms disagree "
                            "(same dims but different orientation/origin)\n"); return 1; }
    }
    return 0;
}

/*--- non-mutating datatype -> float32 extraction (qwarp-private) -----------
 * Applies the NIfTI scale (slope/intercept) like nii_to_float: value =
 * raw*slope + inter, with slope==0 meaning "no scaling" (slope 1, inter 0).
 * Returns a newly malloc'd float[nvox], or NULL on unsupported type / OOM. */
static float *qw_extract_float(const nifti_image *n) {
    int64_t nv = n->nvox;
    float *out = (float *)malloc((size_t)(nv ? nv : 1) * sizeof(float));
    if (!out) { fprintf(stderr, "qwarp: out of memory extracting float image\n"); return NULL; }
    double sl = (n->scl_slope != 0.0) ? n->scl_slope : 1.0;
    double in = n->scl_inter;
    const void *d = n->data;
#define QW_CONV(CTYPE) do { const CTYPE *p = (const CTYPE *)d; \
        for (int64_t i = 0; i < nv; i++) out[i] = (float)((double)p[i]*sl + in); } while (0)
    switch (n->datatype) {
        case DT_UINT8:   QW_CONV(uint8_t);  break;
        case DT_INT8:    QW_CONV(int8_t);   break;
        case DT_INT16:   QW_CONV(int16_t);  break;
        case DT_UINT16:  QW_CONV(uint16_t); break;
        case DT_INT32:   QW_CONV(int32_t);  break;
        case DT_UINT32:  QW_CONV(uint32_t); break;
        case DT_INT64:   QW_CONV(int64_t);  break;
        case DT_UINT64:  QW_CONV(uint64_t); break;
        case DT_FLOAT32: QW_CONV(float);    break;
        case DT_FLOAT64: QW_CONV(double);   break;
        default:
            fprintf(stderr, "qwarp: unsupported input datatype %d\n", n->datatype);
            free(out); return NULL;
    }
#undef QW_CONV
    return out;
}

/*==========================================================================
 * PORTED AFNI IMAGE-PREP ROUTINES (public domain; faithful to the pinned rev).
 * Operate on plain float buffers + int dims (our images are single-channel
 * float on a 3D grid), preserving AFNI's operation order, constants, float
 * precision, and boundary handling.
 *========================================================================*/

/* THD_cliplevel (thd_cliplevel.c), the MRI_float path (our images are float).
 * Iterative histogram cut = mfrac * median of values above the cut. */
static float qw_cliplevel(const float *far, int64_t nvox, float mfrac) {
    if (mfrac <= 0.0f || mfrac >= 0.99f) mfrac = 0.50f;
    const int nhist = 10000;
    double fac = 0.0;
    for (int64_t i = 0; i < nvox; i++) if (far[i] > (float)fac) fac = far[i];  /* mri_max */
    if (fac < 1.e-100) return 0.0f;
    double sfac = nhist / fac;
    int *hist = (int *)calloc((size_t)nhist + 1, sizeof(int));
    if (!hist) return 0.0f;
    double dsum = 0.0; int npos = 0, kk, ii;
    for (int64_t i = 0; i < nvox; i++) {
        if (far[i] > 0.0f) {
            kk = (int)(sfac * far[i] + 0.499);
            if (kk <= nhist) { hist[kk]++; dsum += (double)kk * (double)kk; npos++; }
        }
    }
    if (npos <= 222) { free(hist); return 0.0f; }
    int qq = (int)(0.65f * npos);
    int ib = (int)rint(0.5 * sqrt(dsum / npos));
    for (kk = 0, ii = nhist - 1; ii >= ib && kk < qq; ii--) kk += hist[ii];
    int ncut = ii, nold, nhalf; qq = 0;
    do {
        for (npos = 0, ii = ncut; ii < nhist; ii++) npos += hist[ii];
        nhalf = npos / 2;
        for (kk = 0, ii = ncut; ii < nhist && kk < nhalf; ii++) kk += hist[ii];
        nold = ncut;
        ncut = (int)(mfrac * ii);
        qq++;
    } while (qq < 66 && ncut != nold);
    free(hist);
    fac = ncut / sfac;
    if (fac > 1.e38) fac = 1.e38;
    return (float)fac;
}

/* Gaussian blur — the SOURCE `-blur 0 3` step (and the sigma=4.5 vox step in weightize).
 * AFNI's blur (edt_blur.c, FIR_blur_volume_3d) is GPL-v2 / Medical College of Wisconsin —
 * a bundled EXCEPTION to AFNI's public-domain core — so we do NOT port it. Instead we call
 * niimath's own BSD/public-domain separable Gaussian (coreFLT.c via core32.h), the same
 * backend coreg_fast uses. Voxel units (dx=dy=dz=1) to match AFNI's dx=1 blur convention;
 * the kernel differs from AFNI's, so the result is EQUIVALENT, not bit-identical (allowed). */
static int qw_gauss_blur_vox(float *f, int nx, int ny, int nz, double sigma_vox) {
    if (sigma_vox <= 0.0) return 0;
    float s = (float)sigma_vox;
    return nifti_smooth_gauss_f32(f, nx, ny, nz, 1, 1.0f, 1.0f, 1.0f, s, s, s, -6.0f);
}

/* ---- byte-mask morphology (thd_automask.c, public domain) ---------------- */

/* THD_mask_clust: keep only the largest 6-connected (NN1) cluster of set voxels.
 * Faithful to AFNI's largest-component semantics via an index flood fill. */
static void qw_mask_clust(int nx, int ny, int nz, unsigned char *mmm) {
    int nxy = nx * ny; int64_t nv = (int64_t)nxy * nz;
    int *best = NULL, nbest = 0, cap = 4096, *cur = (int *)malloc(sizeof(int) * cap);
    if (!cur) return;
    for (int64_t s = 0; s < nv; s++) {
        if (!mmm[s]) continue;
        int ncur = 0, head = 0;
        cur[ncur++] = (int)s; mmm[s] = 0;
        while (head < ncur) {
            int ijk = cur[head++];
            int k = ijk / nxy, r = ijk - k * nxy, j = r / nx, i = r - j * nx;
            int nb[6], nn = 0;
            if (i > 0)      nb[nn++] = ijk - 1;
            if (i < nx - 1) nb[nn++] = ijk + 1;
            if (j > 0)      nb[nn++] = ijk - nx;
            if (j < ny - 1) nb[nn++] = ijk + nx;
            if (k > 0)      nb[nn++] = ijk - nxy;
            if (k < nz - 1) nb[nn++] = ijk + nxy;
            for (int t = 0; t < nn; t++) {
                int q = nb[t];
                if (mmm[q]) {
                    mmm[q] = 0;
                    if (ncur >= cap) { cap *= 2; int *nc = (int *)realloc(cur, sizeof(int) * cap);
                        if (!nc) { free(cur); free(best); return; } cur = nc; }
                    cur[ncur++] = q;
                }
            }
        }
        if (ncur > nbest) {
            int *nb2 = (int *)malloc(sizeof(int) * ncur);
            if (!nb2) { free(cur); free(best); return; }
            memcpy(nb2, cur, sizeof(int) * ncur); free(best); best = nb2; nbest = ncur;
        }
    }
    for (int t = 0; t < nbest; t++) mmm[best[t]] = 1;
    free(best); free(cur);
}

/* THD_mask_erode(...,redilate=1,NN=2): erode any in-mask voxel lacking all 18 of its
 * NN1+NN2 neighbors (edges always erode), then redilate eroded voxels adjacent to a
 * survivor. Faithful port of the thd_automask.c body (public domain). */
static void qw_mask_erode(int nx, int ny, int nz, unsigned char *mmm) {
    int nxy = nx * ny; int64_t nxyz = (int64_t)nxy * nz;
    unsigned char *nnn = (unsigned char *)calloc((size_t)nxyz, 1); if (!nnn) return;
    for (int kk = 0; kk < nz; kk++) {
        int kz = kk * nxy, km = (kk == 0) ? kz : kz - nxy, kp = (kk == nz - 1) ? kz : kz + nxy;
        for (int jj = 0; jj < ny; jj++) {
            int jy = jj * nx, jm = (jj == 0) ? jy : jy - nx, jp = (jj == ny - 1) ? jy : jy + nx;
            int jykz = jy + kz;
            for (int ii = 0; ii < nx; ii++) {
                if (!mmm[ii + jykz]) continue;
                if (ii == 0 || jj == 0 || kk == 0 || ii == nx - 1 || jj == ny - 1 || kk == nz - 1) {
                    nnn[ii + jykz] = 1; continue;
                }
                int im = ii - 1, ip = ii + 1, victim = 0;
                int num = mmm[ii + jy + km] + mmm[im + jy + kz] + mmm[ii + jm + kz]
                        + mmm[ii + jp + kz] + mmm[ip + jykz] + mmm[ii + jy + kp];
                if (num < 6) victim = 1;
                if (!victim) {
                    num += mmm[im + jy + km] + mmm[ii + jm + km] + mmm[ii + jp + km] + mmm[ip + jy + km]
                         + mmm[im + jm + kz] + mmm[im + jp + kz] + mmm[ip + jm + kz] + mmm[ip + jp + kz]
                         + mmm[im + jy + kp] + mmm[ii + jm + kp] + mmm[ii + jp + kp] + mmm[ip + jy + kp];
                    if (num < 18) victim = 1;
                }
                nnn[ii + jykz] = (unsigned char)victim;
            }
        }
    }
    for (int64_t i = 0; i < nxyz; i++) if (nnn[i]) mmm[i] = 0;
    for (int kk = 0; kk < nz; kk++) {
        int kz = kk * nxy, km = (kk == 0) ? kz : kz - nxy, kp = (kk == nz - 1) ? kz : kz + nxy;
        for (int jj = 0; jj < ny; jj++) {
            int jy = jj * nx, jm = (jj == 0) ? jy : jy - nx, jp = (jj == ny - 1) ? jy : jy + nx;
            for (int ii = 0; ii < nx; ii++) {
                if (!nnn[ii + jy + kz]) continue;
                int im = (ii == 0) ? 0 : ii - 1, ip = (ii == nx - 1) ? ii : ii + 1;
                int victim = mmm[ii + jy + km] || mmm[im + jy + kz] || mmm[ii + jm + kz]
                          || mmm[ii + jp + kz] || mmm[ip + jy + kz] || mmm[ii + jy + kp];
                if (!victim)
                    victim = mmm[im + jy + km] || mmm[ii + jm + km] || mmm[ii + jp + km] || mmm[ip + jy + km]
                          || mmm[im + jm + kz] || mmm[im + jp + kz] || mmm[ip + jm + kz] || mmm[ip + jp + kz]
                          || mmm[im + jy + kp] || mmm[ii + jm + kp] || mmm[ii + jp + kp] || mmm[ip + jy + kp];
                nnn[ii + jy + kz] = (unsigned char)victim;
            }
        }
    }
    for (int64_t i = 0; i < nxyz; i++) if (nnn[i]) mmm[i] = 1;
    free(nnn);
}

/* THD_mask_erodemany(...,npeel=1): peel layers where an in-mask voxel has < peelthr(17) of
 * 18 NN2 neighbors, then unpeel next to survivors. Faithful port (public domain). */
static void qw_mask_erodemany(int nx, int ny, int nz, unsigned char *mmm, int npeel) {
    int nxy = nx * ny; int64_t nxyz = (int64_t)nxy * nz;
    if (npeel < 1 || nxyz < 27) return;
    const int realpeelthr = 17;
    unsigned char *nnn = (unsigned char *)calloc((size_t)nxyz, 1); if (!nnn) return;
    unsigned char *qqq = (unsigned char *)malloc((size_t)nxyz); if (!qqq) { free(nnn); return; }
    for (int pp = 1; pp <= npeel; pp++) {
        for (int kk = 0; kk < nz; kk++) {
            int kz = kk * nxy, km = (kk == 0) ? kz : kz - nxy, kp = (kk == nz - 1) ? kz : kz + nxy;
            for (int jj = 0; jj < ny; jj++) {
                int jy = jj * nx, jm = (jj == 0) ? jy : jy - nx, jp = (jj == ny - 1) ? jy : jy + nx;
                for (int ii = 0; ii < nx; ii++) {
                    if (!mmm[ii + jy + kz]) continue;
                    int im = (ii == 0) ? 0 : ii - 1, ip = (ii == nx - 1) ? ii : ii + 1;
                    int num = mmm[im + jy + km]
                            + mmm[ii + jm + km] + mmm[ii + jy + km] + mmm[ii + jp + km] + mmm[ip + jy + km]
                            + mmm[im + jm + kz] + mmm[im + jy + kz] + mmm[im + jp + kz]
                            + mmm[ii + jm + kz]                     + mmm[ii + jp + kz]
                            + mmm[ip + jm + kz] + mmm[ip + jy + kz] + mmm[ip + jp + kz]
                            + mmm[im + jy + kp]
                            + mmm[ii + jm + kp] + mmm[ii + jy + kp] + mmm[ii + jp + kp] + mmm[ip + jy + kp];
                    if (num < realpeelthr) nnn[ii + jy + kz] = (unsigned char)pp;
                }
            }
        }
        for (int64_t i = 0; i < nxyz; i++) if (nnn[i]) mmm[i] = 0;
    }
    for (int pp = npeel; pp >= 1; pp--) {
        memset(qqq, 0, (size_t)nxyz);
        int bth = (pp == npeel) ? 0 : 1;
        for (int kk = 0; kk < nz; kk++) {
            int kz = kk * nxy, km = (kk == 0) ? kz : kz - nxy, kp = (kk == nz - 1) ? kz : kz + nxy;
            for (int jj = 0; jj < ny; jj++) {
                int jy = jj * nx, jm = (jj == 0) ? jy : jy - nx, jp = (jj == ny - 1) ? jy : jy + nx;
                for (int ii = 0; ii < nx; ii++) {
                    if (!(nnn[ii + jy + kz] >= pp && !mmm[ii + jy + kz])) continue;
                    int im = (ii == 0) ? 0 : ii - 1, ip = (ii == nx - 1) ? ii : ii + 1;
                    qqq[ii + jy + kz] = (unsigned char)(mmm[im + jy + km]
                        + mmm[ii + jm + km] + mmm[ii + jy + km] + mmm[ii + jp + km] + mmm[ip + jy + km]
                        + mmm[im + jm + kz] + mmm[im + jy + kz] + mmm[im + jp + kz]
                        + mmm[ii + jm + kz]                     + mmm[ii + jp + kz]
                        + mmm[ip + jm + kz] + mmm[ip + jy + kz] + mmm[ip + jp + kz]
                        + mmm[im + jy + kp]
                        + mmm[ii + jm + kp] + mmm[ii + jy + kp] + mmm[ii + jp + kp] + mmm[ip + jy + kp]);
                }
            }
        }
        for (int64_t i = 0; i < nxyz; i++) if (qqq[i] > bth) mmm[i] = 1;
    }
    free(qqq); free(nnn);
}

/* MRI_autobbox (thd_automask.c, PD): inclusive first/last set-voxel index per axis, after
 * the default clust -> erodemany(peelcount=1) -> clust despeckle of the nonzero mask. */
static void qw_autobbox(const float *far, int nx, int ny, int nz,
                        int *xm, int *xp, int *ym, int *yp, int *zm, int *zp) {
    int nxy = nx * ny; int64_t nv = (int64_t)nxy * nz;
    *xm = *xp = *ym = *yp = *zm = *zp = 0;
    unsigned char *mmm = (unsigned char *)calloc((size_t)nv, 1); if (!mmm) return;
    int64_t nmm = 0;
    for (int64_t i = 0; i < nv; i++) if (far[i] != 0.0f) { mmm[i] = 1; nmm++; }
    if (nmm == 0) { free(mmm); return; }
    qw_mask_clust(nx, ny, nz, mmm);
    qw_mask_erodemany(nx, ny, nz, mmm, 1);
    qw_mask_clust(nx, ny, nz, mmm);
    int ii, jj, kk;
    for (ii = 0; ii < nx; ii++) { for (kk = 0; kk < nz; kk++) for (jj = 0; jj < ny; jj++) if (mmm[ii + jj*nx + kk*nxy]) goto X0; } X0: *xm = (ii < nx) ? ii : 0;
    for (ii = nx-1; ii >= 0; ii--) { for (kk = 0; kk < nz; kk++) for (jj = 0; jj < ny; jj++) if (mmm[ii + jj*nx + kk*nxy]) goto X1; } X1: *xp = (ii >= 0) ? ii : 0;
    for (jj = 0; jj < ny; jj++) { for (kk = 0; kk < nz; kk++) for (ii = 0; ii < nx; ii++) if (mmm[ii + jj*nx + kk*nxy]) goto Y0; } Y0: *ym = (jj < ny) ? jj : 0;
    for (jj = ny-1; jj >= 0; jj--) { for (kk = 0; kk < nz; kk++) for (ii = 0; ii < nx; ii++) if (mmm[ii + jj*nx + kk*nxy]) goto Y1; } Y1: *yp = (jj >= 0) ? jj : 0;
    for (kk = 0; kk < nz; kk++) { for (jj = 0; jj < ny; jj++) for (ii = 0; ii < nx; ii++) if (mmm[ii + jj*nx + kk*nxy]) goto Z0; } Z0: *zm = (kk < nz) ? kk : 0;
    for (kk = nz-1; kk >= 0; kk--) { for (jj = 0; jj < ny; jj++) for (ii = 0; ii < nx; ii++) if (mmm[ii + jj*nx + kk*nxy]) goto Z1; } Z1: *zp = (kk >= 0) ? kk : 0;
    free(mmm);
}

/* EDIT_volpad (edt_volpad.c, PD), float only: grow (bot/top >=0) or crop (negative) each
 * face; new voxels are zero. Returns a new float[new nvox] (caller frees) + new dims, or
 * NULL on error. bot/top are per-face voxel counts (negative crops). */
static float *qw_zeropad(const float *f, int nx, int ny, int nz,
                         int xb, int xt, int yb, int yt, int zb, int zt,
                         int *onx, int *ony, int *onz) {
    int nxn = nx + xb + xt, nyn = ny + yb + yt, nzn = nz + zb + zt;
    int ib = (0 > -xb) ? 0 : -xb, it = (nx < nx + xt) ? nx : nx + xt;
    int jb = (0 > -yb) ? 0 : -yb, jt = (ny < ny + yt) ? ny : ny + yt;
    int kb = (0 > -zb) ? 0 : -zb, kt = (nz < nz + zt) ? nz : nz + zt;
    if (nxn < 2 || ib >= it || nyn < 2 || jb >= jt || nzn < 2 || kb >= kt) return NULL;
    int64_t nvn = (int64_t)nxn * nyn * nzn;
    float *o = (float *)calloc((size_t)nvn, sizeof(float)); if (!o) return NULL;
    int64_t nxyo = (int64_t)nx * ny, nxyn = (int64_t)nxn * nyn;
    for (int kk = kb; kk < kt; kk++)
        for (int jj = jb; jj < jt; jj++)
            for (int ii = ib; ii < it; ii++)
                o[(ii + xb) + (int64_t)(jj + yb) * nxn + (int64_t)(kk + zb) * nxyn]
                    = f[ii + (int64_t)jj * nx + (int64_t)kk * nxyo];
    *onx = nxn; *ony = nyn; *onz = nzn;
    return o;
}

/* Data-dependent pad sizes (3dQwarp.c, PD): base-only bounding box at threshold
 * 0.33*cliplevel(base,0.22); mpad_a = rintf(0.1234*n)+1 floored at 9; pad_?m = mpad - bbox_lo,
 * pad_?p = mpad - (n-1 - bbox_hi), each floored at MINPAD=3; z-pad 0 for 2D. */
static void qw_compute_pad(const float *base, int nx, int ny, int nz,
                           int *xm, int *xp, int *ym, int *yp, int *zm, int *zp) {
    int64_t nv = (int64_t)nx * ny * nz;
    float cv = 0.33f * qw_cliplevel(base, nv, 0.22f);
    float *q = (float *)malloc(sizeof(float) * (size_t)nv);
    if (!q) { *xm=*xp=*ym=*yp=*zm=*zp=3; return; }
    for (int64_t i = 0; i < nv; i++) q[i] = (base[i] < cv) ? 0.0f : base[i];
    int bxm, bxp, bym, byp, bzm, bzp;
    qw_autobbox(q, nx, ny, nz, &bxm, &bxp, &bym, &byp, &bzm, &bzp);
    free(q);
    int mpx = (int)rintf(0.1234f * nx) + 1; if (mpx < 9) mpx = 9;
    int mpy = (int)rintf(0.1234f * ny) + 1; if (mpy < 9) mpy = 9;
    int mpz = (int)rintf(0.1234f * nz) + 1; if (mpz < 9) mpz = 9;
    const int MINPAD = 3;
    *xm = mpx - bxm;                 if (*xm <= MINPAD - 1) *xm = MINPAD;
    *ym = mpy - bym;                 if (*ym <= MINPAD - 1) *ym = MINPAD;
    *zm = mpz - bzm;                 if (*zm <= MINPAD - 1) *zm = MINPAD;
    *xp = mpx - (nx - 1 - bxp);      if (*xp <= MINPAD - 1) *xp = MINPAD;
    *yp = mpy - (ny - 1 - byp);      if (*yp <= MINPAD - 1) *yp = MINPAD;
    *zp = mpz - (nz - 1 - bzp);      if (*zp <= MINPAD - 1) *zp = MINPAD;
    if (nz == 1) { *zm = *zp = 0; }
}

/* ---- clean-room spherical median (replaces the GPL mri_medianfilter chain) ----
 * mri_medianfilter.c is public domain but calls the GPL MCW_build_mask (sphere offsets)
 * and qmed_float (median). Both are standard algorithms reimplemented here independently:
 * the sphere is the geometric set {d : |d|^2 <= r^2, d != 0} plus the center; the median is
 * a plain sort of the in-grid, in-mask neighborhood (even count -> mean of the two middle,
 * matching AFNI's convention). Applied only where mask!=0; elsewhere the output is 0. */
static int qw_cmp_flt(const void *a, const void *b) {
    float x = *(const float *)a, y = *(const float *)b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}
static float *qw_median_sphere(const float *data, int nx, int ny, int nz,
                               const unsigned char *mask, double radius) {
    int64_t nv = (int64_t)nx * ny * nz;
    int rr = (int)radius; double r2 = radius * radius;
    int cap = (2 * rr + 1) * (2 * rr + 1) * (2 * rr + 1);
    int *di = malloc(sizeof(int) * cap), *dj = malloc(sizeof(int) * cap), *dk = malloc(sizeof(int) * cap);
    if (!di || !dj || !dk) { free(di); free(dj); free(dk); return NULL; }
    int nd = 0;
    for (int kk = -rr; kk <= rr; kk++) for (int jj = -rr; jj <= rr; jj++) for (int ii = -rr; ii <= rr; ii++) {
        double q = (double)ii*ii + (double)jj*jj + (double)kk*kk;
        if (q <= r2 && q > 0.0) { di[nd] = ii; dj[nd] = jj; dk[nd] = kk; nd++; }
    }
    di[nd] = 0; dj[nd] = 0; dk[nd] = 0; nd++;   /* add the center (AFNI ADDTO_CLUSTER 0,0,0) */
    float *out = (float *)calloc((size_t)nv, sizeof(float));
    float *tmp = (float *)malloc(sizeof(float) * nd);
    if (!out || !tmp) { free(out); free(tmp); free(di); free(dj); free(dk); return NULL; }
    int nxy = nx * ny;
    for (int kk = 0; kk < nz; kk++) for (int jj = 0; jj < ny; jj++) for (int ii = 0; ii < nx; ii++) {
        int64_t ijk = ii + (int64_t)jj * nx + (int64_t)kk * nxy;
        if (!mask[ijk]) continue;
        int nt = 0;
        for (int d = 0; d < nd; d++) {
            int ip = ii + di[d]; if (ip < 0 || ip >= nx) continue;
            int jp = jj + dj[d]; if (jp < 0 || jp >= ny) continue;
            int kp = kk + dk[d]; if (kp < 0 || kp >= nz) continue;
            int64_t pjk = ip + (int64_t)jp * nx + (int64_t)kp * nxy;
            if (!mask[pjk]) continue;
            tmp[nt++] = data[pjk];
        }
        if (nt == 0) continue;
        qsort(tmp, nt, sizeof(float), qw_cmp_flt);
        out[ijk] = (nt & 1) ? tmp[nt/2] : 0.5f * (tmp[nt/2] + tmp[nt/2 - 1]);
    }
    free(tmp); free(di); free(dj); free(dk);
    return out;
}

/* mri_weightize(base, acod=1, ndil=5, aclip=0, apow=1) — 3dQwarp.c, public domain.
 * Graded base weight in [0,1]: |base| -> edge fade -> squash>3*cliplevel(0.5) -> median r=2.25
 * -> Gaussian sigma 4.5 vox -> lower-clip + largest-cluster (nz>2) -> normalize. The GPL
 * median/blur are substituted (clean-room median, niimath Gaussian); the rest is the PD
 * structure. Returns float[nvox] in [0,1] (caller frees), or NULL on failure. */
static float *qw_weightize(const float *base, int nx, int ny, int nz) {
    int64_t nv = (int64_t)nx * ny * nz;
    float *wf = (float *)malloc(sizeof(float) * (size_t)nv);
    if (!wf) return NULL;
    for (int64_t i = 0; i < nv; i++) wf[i] = fabsf(base[i]);
    /* edge fade-to-zero */
    int nxy = nx * ny;
    int xf = (int)(0.04 * nx + 2.0); if (6 * xf >= nx) xf = (nx - 1) / 6;
    int yf = (int)(0.04 * ny + 2.0); if (6 * yf >= ny) yf = (ny - 1) / 6;
    int zf = (int)(0.04 * nz + 2.0); if (6 * zf >= nz) zf = (nz - 1) / 6;
    for (int jj = 0; jj < ny; jj++) for (int ii = 0; ii < nx; ii++) for (int ff = 0; ff < zf; ff++) {
        wf[ii + (int64_t)jj*nx + (int64_t)ff*nxy] = 0.0f;
        wf[ii + (int64_t)jj*nx + (int64_t)(nz-1-ff)*nxy] = 0.0f; }
    for (int kk = 0; kk < nz; kk++) for (int jj = 0; jj < ny; jj++) for (int ff = 0; ff < xf; ff++) {
        wf[ff + (int64_t)jj*nx + (int64_t)kk*nxy] = 0.0f;
        wf[(nx-1-ff) + (int64_t)jj*nx + (int64_t)kk*nxy] = 0.0f; }
    for (int kk = 0; kk < nz; kk++) for (int ii = 0; ii < nx; ii++) for (int ff = 0; ff < yf; ff++) {
        wf[ii + (int64_t)ff*nx + (int64_t)kk*nxy] = 0.0f;
        wf[ii + (int64_t)(ny-1-ff)*nx + (int64_t)kk*nxy] = 0.0f; }
    /* squash super-large values down */
    float clip = 3.0f * qw_cliplevel(wf, nv, 0.5f);
    for (int64_t i = 0; i < nv; i++) if (wf[i] > clip) wf[i] = clip;
    /* median (r=2.25) over the positive mask */
    unsigned char *mmm = (unsigned char *)malloc((size_t)nv);
    if (!mmm) { free(wf); return NULL; }
    for (int64_t i = 0; i < nv; i++) mmm[i] = (wf[i] > 0.0f);
    float *wm = qw_median_sphere(wf, nx, ny, nz, mmm, 2.25);
    free(wf); if (!wm) { free(mmm); return NULL; }
    wf = wm;
    /* Gaussian blur sigma 4.5 voxels (niimath) — a failure would silently skip the required
     * smoothing, so treat it as an error (the caller checks for a NULL weight). */
    if (qw_gauss_blur_vox(wf, nx, ny, nz, 4.5) != 0) { free(wf); free(mmm); return NULL; }
    /* lower clip + keep largest cluster (only for true 3D) */
    float mx = 0.0f; for (int64_t i = 0; i < nv; i++) if (wf[i] > mx) mx = wf[i];
    clip = 0.05f * mx;
    float clip2 = 0.33f * qw_cliplevel(wf, nv, 0.33f);
    if (clip2 > clip) clip = clip2;
    if (nz > 2) {
        for (int64_t i = 0; i < nv; i++) mmm[i] = (wf[i] >= clip);
        qw_mask_clust(nx, ny, nz, mmm);
        qw_mask_erode(nx, ny, nz, mmm);
        qw_mask_clust(nx, ny, nz, mmm);
        for (int64_t i = 0; i < nv; i++) if (!mmm[i]) wf[i] = 0.0f;
    }
    free(mmm);
    /* normalize to [0,1] */
    mx = 0.0f; for (int64_t i = 0; i < nv; i++) if (wf[i] > mx) mx = wf[i];
    if (mx <= 0.0f) { free(wf); return NULL; }
    float inv = 1.0f / mx;
    for (int64_t i = 0; i < nv; i++) wf[i] *= inv;
    return wf;
}

/* ---- IndexWarp3D: the dense displacement field (mri_nwarp.c, public domain) --------
 * Voxel (i,j,k) on the base grid maps to SOURCE sampling location (i+xd, j+yd, k+zd) —
 * displacements in INDEX (voxel) units; identity warp = all zeros. Only xd/yd/zd are used
 * on the source->base (s2bim) path; AFNI's amat (affine part), hv/je/se (energy) and
 * cmat/geomstring (dataset I/O) fields are omitted here (unused for a plain warp). Private
 * to qwarp.c — future work (warp output, transform composition) will surface it deliberately. */
typedef struct {
    int nx, ny, nz;
    int64_t nvox;
    float *xd, *yd, *zd;   /* index-unit displacements, nvox each */
    float *je, *se;        /* per-voxel bulk/shear penalty energies (lazy; NULL until loaded) */
} qw_warp;

static void qw_warp_destroy(qw_warp *w) {
    if (!w) return;
    free(w->xd); free(w->yd); free(w->zd); free(w->je); free(w->se); free(w);
}
static qw_warp *qw_warp_create(int nx, int ny, int nz) {   /* IW3D_create: zero (identity) warp */
    qw_warp *w = (qw_warp *)calloc(1, sizeof *w);
    if (!w) return NULL;
    w->nx = nx; w->ny = ny; w->nz = nz; w->nvox = (int64_t)nx * ny * nz;
    w->xd = (float *)calloc((size_t)w->nvox, sizeof(float));
    w->yd = (float *)calloc((size_t)w->nvox, sizeof(float));
    w->zd = (float *)calloc((size_t)w->nvox, sizeof(float));
    if (!w->xd || !w->yd || !w->zd) { qw_warp_destroy(w); return NULL; }
    return w;   /* je/se stay NULL until IW3D_load_energy allocates them */
}
/* ---- deformation penalty energy (mri_nwarp.c hexahedron_energy / IW3D_load_energy /
 * HPEN_addup, public domain). Faithful port: bulk energy en.a = 0.333*(J-1/J)^2, shear+
 * vorticity en.b, per-voxel penalty = max(je-cut,0)^4 + max(se-cut,0)^4 (cut = 1). ---- */
#define QW_HPEN_CUT 1.0f

/* hexahedron_energy: strain from 8 corner displacements d[0..7] (000,100,010,110,001,101,
 * 011,111), each a float[3] (x,y,z). Writes {bulk, shear} into en[2]. */
static void qw_hexahedron_energy(const float d[8][3], float en[2]) {
#define QW_DA(p,q) (d[p][0]-d[q][0])
#define QW_DB(p,q) (d[p][1]-d[q][1])
#define QW_DC(p,q) (d[p][2]-d[q][2])
    float fxx = (QW_DA(1,0) + QW_DA(7,6)) * 0.5f + 1.0f;
    float fxy = (QW_DB(1,0) + QW_DB(7,6)) * 0.5f;
    float fxz = (QW_DC(1,0) + QW_DC(7,6)) * 0.5f;
    float fyx = (QW_DA(2,0) + QW_DA(7,5)) * 0.5f;
    float fyy = (QW_DB(2,0) + QW_DB(7,5)) * 0.5f + 1.0f;
    float fyz = (QW_DC(2,0) + QW_DC(7,5)) * 0.5f;
    float fzx = (QW_DA(4,0) + QW_DA(7,3)) * 0.5f;
    float fzy = (QW_DB(4,0) + QW_DB(7,3)) * 0.5f;
    float fzz = (QW_DC(4,0) + QW_DC(7,3)) * 0.5f + 1.0f;
#undef QW_DA
#undef QW_DB
#undef QW_DC
    float JJ = fxx*(fyy*fzz-fyz*fzy) + fyx*(fzy*fxz-fzz*fxy) + fzx*(fxy*fyz-fxz*fyy);
    if (JJ < 0.1f) JJ = 0.1f; else if (JJ > 10.0f) JJ = 10.0f;
    float II = fxx*fxx + fyy*fyy + fzz*fzz + fxy*fxy + fyx*fyx
             + fxz*fxz + fzx*fzx + fyz*fyz + fzy*fzy;
    float vx = fyz - fzy, vy = fxz - fzx, vz = fxy - fyx;
    float VV = 2.0f * (vx*vx + vy*vy + vz*vz);
    float jcb = cbrtf(JJ*JJ);
    II = (II + VV) / jcb - 3.0f; if (II < 0.0f) II = 0.0f;
    jcb = JJ - 1.0f/JJ;
    en[0] = 0.333f * jcb * jcb; en[1] = II;
}

/* IW3D_load_energy: fill w->je/se from the displacement field; return the summed penalty. */
static double qw_load_energy(qw_warp *w) {
    if (!w) return 0.0;
    int nx = w->nx, ny = w->ny, nz = w->nz, nxy = nx*ny; int64_t nxyz = w->nvox;
    const float *xda = w->xd, *yda = w->yd, *zda = w->zd;
    if (!w->je) w->je = (float *)calloc((size_t)nxyz, sizeof(float));
    if (!w->se) w->se = (float *)calloc((size_t)nxyz, sizeof(float));
    if (!w->je || !w->se) return 0.0;
    float *jea = w->je, *sea = w->se;
    double esum = 0.0;
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) reduction(+:esum) if(nxyz > 9999)
#endif
    for (int64_t qq = 0; qq < nxyz; qq++) {
        int ii = qq % nx, kk = (int)(qq / nxy), jj = (int)((qq - (int64_t)kk*nxy) / nx);
        int ip = ii+1, jp = jj+1, kp = kk+1;
        if (ip == nx) ip--; if (jp == ny) jp--; if (kp == nz) kp--;
        int64_t c[8];
        c[0] = qq;
        c[1] = ip + (int64_t)jj*nx + (int64_t)kk*nxy;
        c[2] = ii + (int64_t)jp*nx + (int64_t)kk*nxy;
        c[3] = ip + (int64_t)jp*nx + (int64_t)kk*nxy;
        c[4] = ii + (int64_t)jj*nx + (int64_t)kp*nxy;
        c[5] = ip + (int64_t)jj*nx + (int64_t)kp*nxy;
        c[6] = ii + (int64_t)jp*nx + (int64_t)kp*nxy;
        c[7] = ip + (int64_t)jp*nx + (int64_t)kp*nxy;
        float d[8][3];
        for (int m = 0; m < 8; m++) { d[m][0] = xda[c[m]]; d[m][1] = yda[c[m]]; d[m][2] = zda[c[m]]; }
        float en[2]; qw_hexahedron_energy(d, en);
        jea[qq] = en[0]; sea[qq] = en[1];
        float ev = en[0] - QW_HPEN_CUT; if (ev > 0.0f) esum += (double)(ev*ev)*(ev*ev);
        ev = en[1] - QW_HPEN_CUT; if (ev > 0.0f) esum += (double)(ev*ev)*(ev*ev);
    }
    return esum;
}

/* HPEN_addup: penalty sum over pre-computed je/se arrays (out-of-patch seed). */
static double qw_hpen_addup(int64_t n, const float *je, const float *se) {
    double esum = 0.0;
    for (int64_t i = 0; i < n; i++) {
        float ev = je[i] - QW_HPEN_CUT; if (ev > 0.0f) esum += (double)(ev*ev)*(ev*ev);
        ev = se[i] - QW_HPEN_CUT; if (ev > 0.0f) esum += (double)(ev*ev)*(ev*ev);
    }
    return esum;
}

/* ---- interpolation (mri_genalign_util.c, public domain) ------------------
 * Codes: QW_LINEAR (fast, used during optimization) and QW_WSINC5 (the single
 * final reslice). Both fill out-of-FOV points with QW_OUTVAL and CLAMP in-range
 * near-edge taps to [0,n-1] (edge replication) — faithful to GA_interp_*. */
#define QW_LINEAR  1
#define QW_WSINC5  2
#define QW_OUTVAL  0.0f
#define QW_FAR(i,j,k) far[(i) + (j)*nx + (int64_t)(k)*nxy]

static void qw_interp_linear(const float *far, int nx, int ny, int nz,
                             int np, const float *ip, const float *jp, const float *kp, float *vv) {
    int nxy = nx * ny, nx1 = nx - 1, ny1 = ny - 1, nz1 = nz - 1;
    float nxh = nx - 0.501f, nyh = ny - 0.501f, nzh = nz - 0.501f;
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) if(np > 4444)
#endif
    for (int pp = 0; pp < np; pp++) {
        float xx = ip[pp]; if (xx < -0.499f || xx > nxh) { vv[pp] = QW_OUTVAL; continue; }
        float yy = jp[pp]; if (yy < -0.499f || yy > nyh) { vv[pp] = QW_OUTVAL; continue; }
        float zz = kp[pp]; if (zz < -0.499f || zz > nzh) { vv[pp] = QW_OUTVAL; continue; }
        float ixf = floorf(xx), fx = xx - ixf, jyf = floorf(yy), fy = yy - jyf, kzf = floorf(zz), fz = zz - kzf;
        int i0 = (int)ixf, i1 = i0 + 1; if (i0 < 0) i0 = 0; else if (i0 > nx1) i0 = nx1; if (i1 < 0) i1 = 0; else if (i1 > nx1) i1 = nx1;
        int j0 = (int)jyf, j1 = j0 + 1; if (j0 < 0) j0 = 0; else if (j0 > ny1) j0 = ny1; if (j1 < 0) j1 = 0; else if (j1 > ny1) j1 = ny1;
        int k0 = (int)kzf, k1 = k0 + 1; if (k0 < 0) k0 = 0; else if (k0 > nz1) k0 = nz1; if (k1 < 0) k1 = 0; else if (k1 > nz1) k1 = nz1;
        float w0 = 1.0f - fx, w1 = fx;
#define QW_XINT(j,k) (w0*QW_FAR(i0,j,k) + w1*QW_FAR(i1,j,k))
        float a00 = QW_XINT(j0,k0), a10 = QW_XINT(j1,k0), a01 = QW_XINT(j0,k1), a11 = QW_XINT(j1,k1);
#undef QW_XINT
        w0 = 1.0f - fy; w1 = fy;
        float b0 = w0*a00 + w1*a10, b1 = w0*a01 + w1*a11;
        vv[pp] = (1.0f - fz) * b0 + fz * b1;
    }
}

/* WSINC5 default config (no AFNI_WSINC5_* env): IRAD=5 (10 taps/axis, 1000-tap cube),
 * WRAD=5.001, taper = M3 minimum-sidelobe window, WCUT=0; value = sum / (wx*wy*wz). */
#define QW_IRAD  5
#define QW_IRAD1 4
#define QW_WRAD  5.001f
static inline float qw_sinc(float x) {
    return (x > 0.01f) ? sinf(3.1415927f * x) / (3.1415927f * x) : 1.0f - 1.6449341f * x * x;
}
static inline float qw_taper_m3(float x) {   /* M3(x) minimum-sidelobe 3-term window */
    return 0.4243801f + 0.4973406f * cosf(3.1415927f * x) + 0.0782793f * cosf(3.1415927f * x * 2.0f);
}

static void qw_interp_wsinc5(const float *far, int nx, int ny, int nz,
                             int np, const float *ip, const float *jp, const float *kp, float *vv) {
    int nxy = nx * ny, nx1 = nx - 1, ny1 = ny - 1, nz1 = nz - 1;
    float nxh = nx - 0.501f, nyh = ny - 0.501f, nzh = nz - 0.501f;
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) if(np > 444)
#endif
    for (int pp = 0; pp < np; pp++) {
        float xx = ip[pp]; if (xx < -0.499f || xx > nxh) { vv[pp] = QW_OUTVAL; continue; }
        float yy = jp[pp]; if (yy < -0.499f || yy > nyh) { vv[pp] = QW_OUTVAL; continue; }
        float zz = kp[pp]; if (zz < -0.499f || zz > nzh) { vv[pp] = QW_OUTVAL; continue; }
        int ix = (int)floorf(xx); float fx = xx - ix;
        int jy = (int)floorf(yy); float fy = yy - jy;
        int kz = (int)floorf(zz); float fz = zz - kz;
        if (fabsf(fx) < 0.0001f && fabsf(fy) < 0.0001f && fabsf(fz) < 0.0001f) {
            int i = ix, j = jy, k = kz;
            if (i < 0) i = 0; else if (i > nx1) i = nx1;
            if (j < 0) j = 0; else if (j > ny1) j = ny1;
            if (k < 0) k = 0; else if (k > nz1) k = nz1;
            vv[pp] = QW_FAR(i, j, k); continue;
        }
        const int NT = 2 * QW_IRAD;   /* 10 */
        float wtt[2 * QW_IRAD]; int iqq[2 * QW_IRAD];
        float wsum = 0.0f;
        for (int qq = -QW_IRAD1; qq <= QW_IRAD; qq++) {
            float xw = fabsf(fx - qq), wt = qw_sinc(xw);
            xw /= QW_WRAD; if (xw > 0.0f) wt *= qw_taper_m3(xw);
            wtt[qq + QW_IRAD1] = wt; wsum += wt;
            int iq = ix + qq; if (iq < 0) iq = 0; else if (iq > nx1) iq = nx1; iqq[qq + QW_IRAD1] = iq;
        }
        float wfac = wsum;
        float fjk[2 * QW_IRAD][2 * QW_IRAD];
        for (int jj = -QW_IRAD1; jj <= QW_IRAD; jj++) {
            int jq = jy + jj; if (jq < 0) jq = 0; else if (jq > ny1) jq = ny1;
            for (int kk = -QW_IRAD1; kk <= QW_IRAD; kk++) {
                int kq = kz + kk; if (kq < 0) kq = 0; else if (kq > nz1) kq = nz1;
                float sum = 0.0f;
                for (int t = 0; t < NT; t++) sum += QW_FAR(iqq[t], jq, kq) * wtt[t];
                fjk[jj + QW_IRAD1][kk + QW_IRAD1] = sum;
            }
        }
        wsum = 0.0f;
        for (int qq = -QW_IRAD1; qq <= QW_IRAD; qq++) {
            float yw = fabsf(fy - qq), wt = qw_sinc(yw);
            yw /= QW_WRAD; if (yw > 0.0f) wt *= qw_taper_m3(yw);
            wtt[qq + QW_IRAD1] = wt; wsum += wt;
        }
        wfac *= wsum;
        float fk[2 * QW_IRAD];
        for (int kk = 0; kk < NT; kk++) {
            float sum = 0.0f;
            for (int jj = 0; jj < NT; jj++) sum += wtt[jj] * fjk[jj][kk];
            fk[kk] = sum;
        }
        wsum = 0.0f;
        for (int qq = -QW_IRAD1; qq <= QW_IRAD; qq++) {
            float zw = fabsf(fz - qq), wt = qw_sinc(zw);
            zw /= QW_WRAD; if (zw > 0.0f) wt *= qw_taper_m3(zw);
            wtt[qq + QW_IRAD1] = wt; wsum += wt;
        }
        wfac *= wsum;
        float sum = 0.0f;
        for (int kk = 0; kk < NT; kk++) sum += wtt[kk] * fk[kk];
        vv[pp] = sum / wfac;
    }
}

/* IW3D_warp_floatim / _into_floatim / THD_interp_floatim (mri_nwarp.c, PD): warp `src`
 * (on the warp's grid) by the displacement field -> a new float image. Base voxel (i,j,k)
 * samples the source at (i+xd*fac, j+yd*fac, k+zd*fac). High-order (WSINC5) output is
 * clamped to the source intensity range (sinc-overshoot suppression). Caller frees. */
static float *qw_warp_apply(const qw_warp *w, const float *src, int nx, int ny, int nz,
                            int code, float fac) {
    int64_t nv = (int64_t)nx * ny * nz;
    float *ip = (float *)malloc(sizeof(float) * (size_t)nv);
    float *jp = (float *)malloc(sizeof(float) * (size_t)nv);
    float *kp = (float *)malloc(sizeof(float) * (size_t)nv);
    float *out = (float *)malloc(sizeof(float) * (size_t)nv);
    if (!ip || !jp || !kp || !out) { free(ip); free(jp); free(kp); free(out); return NULL; }
    int nxy = nx * ny;
    for (int kk = 0; kk < nz; kk++) for (int jj = 0; jj < ny; jj++) for (int ii = 0; ii < nx; ii++) {
        int64_t ijk = ii + (int64_t)jj * nx + (int64_t)kk * nxy;
        ip[ijk] = ii + w->xd[ijk] * fac;
        jp[ijk] = jj + w->yd[ijk] * fac;
        kp[ijk] = kk + w->zd[ijk] * fac;
    }
    if (code == QW_WSINC5) qw_interp_wsinc5(src, nx, ny, nz, (int)nv, ip, jp, kp, out);
    else                   qw_interp_linear(src, nx, ny, nz, (int)nv, ip, jp, kp, out);
    free(ip); free(jp); free(kp);
    if (code == QW_WSINC5) {   /* clamp to source [min,max] (high-order overshoot) */
        float bot = src[0], top = src[0];
        for (int64_t i = 1; i < nv; i++) { if (src[i] < bot) bot = src[i]; else if (src[i] > top) top = src[i]; }
        for (int64_t i = 0; i < nv; i++) { if (out[i] < bot) out[i] = bot; else if (out[i] > top) out[i] = top; }
    }
    return out;
}

/* ---- INCOR similarity: incomplete (clipped-)Pearson (thd_incorrelate.c, PD) --------
 * The DEFAULT 3dQwarp cost. A running-sum accumulator (sx..sw); a fixed base is primed
 * once via _addto, then each warp candidate is scored by cloning + adding the varying
 * data (qw_incor_eval). PEARCLP soft-clips each coord to [cbot,ctop] (mapping a clipped
 * value to a destination just beyond the threshold) and down-weights by 1/(#clipped+1),
 * centering by the clip-window midpoint; the finalizer is atanh(r), clamped to +/-4,
 * POSITIVE = better (the caller negates for a minimizer). Plain Pearson (no clip) is the
 * fallback used when either image has negative values (clip bounds are undefined then). */
#define QW_INCOR_PEARSON 1
#define QW_INCOR_PEARCLP 3

typedef struct {
    int meth;
    double sx, sxx, sy, syy, sxy, sw;                 /* weighted running sums */
    double xcbot, xctop, ycbot, yctop;                /* clip thresholds (x=base, y=source) */
    double xdbot, xdtop, ydbot, ydtop;                /* clip destinations */
} qw_incor;

static double qw_myatanh(double x) {   /* MYatanh: atanh clamped to +/-4 for |x|>0.9993293 */
    return (x < -0.9993293) ? -4.0 : (x > 0.9993293) ? 4.0 : atanh(x);
}

static void qw_incor_addto(qw_incor *c, int n, const float *x, const float *y, const float *w) {
    double sx = c->sx, sxx = c->sxx, sy = c->sy, syy = c->syy, sxy = c->sxy, sw = c->sw;
    if (c->meth == QW_INCOR_PEARCLP) {
        double xcb = c->xcbot, xct = c->xctop, ycb = c->ycbot, yct = c->yctop;
        double xmid = 0.5 * (xcb + xct), ymid = 0.5 * (ycb + yct);
        double xdb = c->xdbot, xdt = c->xdtop, ydb = c->ydbot, ydt = c->ydtop;
#ifdef _OPENMP
        #pragma omp parallel for schedule(static) reduction(+:sx,sxx,sy,syy,sxy,sw) if(n > 9999)
#endif
        for (int ii = 0; ii < n; ii++) {
            double ww = w ? (double)w[ii] : 1.0;
            if (ww <= 0.0) continue;
            int cl = 1;
            double xx = (double)x[ii];
            if (xx <= xcb) { xx = xdb; cl++; } else if (xx >= xct) { xx = xdt; cl++; }
            double yy = (double)y[ii];
            if (yy <= ycb) { yy = ydb; cl++; } else if (yy >= yct) { yy = ydt; cl++; }
            ww /= cl; xx -= xmid; yy -= ymid;
            sx += xx*ww; sxx += xx*xx*ww; sy += yy*ww; syy += yy*yy*ww; sxy += xx*yy*ww; sw += ww;
        }
    } else {   /* plain Pearson: raw weighted sums, no clip/centering */
#ifdef _OPENMP
        #pragma omp parallel for schedule(static) reduction(+:sx,sxx,sy,syy,sxy,sw) if(n > 9999)
#endif
        for (int ii = 0; ii < n; ii++) {
            double ww = w ? (double)w[ii] : 1.0;
            if (w && ww <= 0.0) continue;
            double xx = (double)x[ii], yy = (double)y[ii];
            sx += xx*ww; sxx += xx*xx*ww; sy += yy*ww; syy += yy*yy*ww; sxy += xx*yy*ww; sw += ww;
        }
    }
    c->sx = sx; c->sxx = sxx; c->sy = sy; c->syy = syy; c->sxy = sxy; c->sw = sw;
}

static double qw_incor_finalize(const qw_incor *c) {
    if (c->sw <= 0.0) return 0.0;
    double swi = 1.0 / c->sw;
    double xv = c->sxx - c->sx * c->sx * swi;
    double yv = c->syy - c->sy * c->sy * swi;
    double xy = c->sxy - c->sx * c->sy * swi;
    if (xv <= 0.0 || yv <= 0.0) return 0.0;
    return qw_myatanh(xy / sqrt(xv * yv));
}

/* Score the varying (x,y,w) against a primed base WITHOUT modifying the base: clone (a plain
 * struct copy carries the sums + clip bounds), fold in the varying data, finalize. */
static double qw_incor_eval(const qw_incor *base, int n, const float *x, const float *y, const float *w) {
    qw_incor c = *base;
    qw_incor_addto(&c, n, x, y, w);
    return qw_incor_finalize(&c);
}

/* mri_quantile(alpha): the alpha-quantile of the good (finite) values, linear-interpolated. */
static float qw_quantile(const float *a, int n, double alpha) {
    if (n < 1) return 0.0f;
    float *s = (float *)malloc(sizeof(float) * n); if (!s) return 0.0f;
    memcpy(s, a, sizeof(float) * n);
    qsort(s, n, sizeof(float), qw_cmp_flt);   /* shared comparator (identical to the removed qw_cmp_flt2) */
    double pos = alpha * (n - 1); int lo = (int)pos; double fr = pos - lo;
    float v = (lo + 1 < n) ? (float)((1.0 - fr) * s[lo] + fr * s[lo + 1]) : s[n - 1];
    free(s); return v;
}

/* INCOR_clipate (thd_incorrelate.c): bottom = cliplevel(0.321), top = 98.7th pct capped at
 * 6.543*bottom; needs >=666 good values and a non-negative image, else degenerate (bot>=top). */
static void qw_clipate(const float *xar, int64_t n, float *cbot, float *ctop) {
    float *buf = (float *)malloc(sizeof(float) * (size_t)n); if (!buf) { *cbot = 1.0f; *ctop = 0.0f; return; }
    int nq = 0; float mn = 0.0f;
    for (int64_t i = 0; i < n; i++) if (xar[i] < 1.0e10f) { if (nq == 0 || xar[i] < mn) mn = xar[i]; buf[nq++] = xar[i]; }
    if (nq < 666 || mn < 0.0f) { *cbot = 1.0f; *ctop = 0.0f; free(buf); return; }
    float cb = qw_cliplevel(buf, nq, 0.321f);
    float ct = qw_quantile(buf, nq, 0.987);
    if (ct > 6.543f * cb) ct = 6.543f * cb;
    *cbot = cb; *ctop = ct; free(buf);
}

/* Fill a PEARCLP incor's 8 clip params from the (fixed) base x[] and source y[] samples,
 * reproducing 3dQwarp's mpar mapping (clip thresholds = quantile bounds; destinations sit a
 * little beyond). Returns 0 on success, nonzero if the clip is degenerate (caller -> Pearson). */
static int qw_incor_setup_clip(qw_incor *c, const float *x, const float *y, int64_t n) {
    float xcb, xct, ycb, yct;
    qw_clipate(x, n, &xcb, &xct);
    qw_clipate(y, n, &ycb, &yct);
    if (xcb >= xct || ycb >= yct) return 1;   /* degenerate (e.g. negative image) */
    float xmn = x[0], xmx = x[0], ymn = y[0], ymx = y[0];
    for (int64_t i = 1; i < n; i++) {
        if (x[i] < xmn) xmn = x[i]; else if (x[i] > xmx) xmx = x[i];
        if (y[i] < ymn) ymn = y[i]; else if (y[i] > ymx) ymx = y[i];
    }
    float d1, d2, dif;
    d2 = 0.05f * (xct - xcb);
    d1 = 0.5f * (xcb - xmn); dif = (d1 < d2) ? d1 : d2; c->xdbot = xcb - dif;
    d1 = 0.5f * (xmx - xct); dif = (d1 < d2) ? d1 : d2; c->xdtop = xct + dif;
    d2 = 0.05f * (yct - ycb);
    d1 = 0.5f * (ycb - ymn); dif = (d1 < d2) ? d1 : d2; c->ydbot = ycb - dif;
    d1 = 0.5f * (ymx - yct); dif = (d1 < d2) ? d1 : d2; c->ydtop = yct + dif;
    c->xcbot = xcb; c->xctop = xct; c->ycbot = ycb; c->yctop = yct;
    return 0;
}

/* ======================================================================
 * PATCH OPTIMIZER + WARP DRIVER (mri_nwarp.c, public domain) — the full
 * IW3D_warpomatic engine: level-0 basis escalation (SINCC->CUBIC_LITE->CUBIC
 * ->QUINTIC_LITE->QUINTIC), the shrinking-patch level schedule, alternating
 * sweeps, and the deformation penalty for levels >= 3. AFNI's ~150 "H..."
 * file-scope globals are centralized into one qw_ctx; a single process-global
 * pointer g_qw is read by the NEWUOA cost callback (serial-only, mirroring
 * coreg_fast.c's g_cf).
 *
 * DEVIATION (documented): AFNI toggles the NEWUOA constraint region between a
 * box and a ball (BOXOPT/BALLOPT). The standalone powell_newuoa.c drop-in has
 * no exposed setter for its box/ball mode (con_meth is a file-static defaulting
 * to box) and is byte-for-byte niimath-owned, so it cannot be edited here. We
 * therefore emulate AFNI's box/ball SELECTION for the per-basis max-displacement
 * scale (Hbasis_parmax) but always run the box-constrained optimizer. This is a
 * numerical (not byte) match, consistent with the plan's equivalence contract.
 * ==================================================================== */
#define QW_CUBIC        2     /* AFNI MRI_CUBIC        */
#define QW_QUINTIC      4     /* AFNI MRI_QUINTIC      */
#define QW_SINCC        311   /* AFNI MRI_SINCC        */
#define QW_CUBIC_LITE   398   /* AFNI MRI_CUBIC_LITE   */
#define QW_QUINTIC_LITE 399   /* AFNI MRI_QUINTIC_LITE */

#define QW_NGMIN    5         /* absolute smallest patch edge          */
#define QW_NGMIN_Q  7         /* smallest quintic patch edge           */
#define QW_HNGMIN   25        /* default minimum patch (3dQwarp Hngmin) */
#define QW_SHRINK   0.749999f /* inter-level patch shrink factor        */
#define QW_SC_BOX   1
#define QW_SC_BALL  2
#define QW_HPEN_FBASE 0.033333 /* base penalty factor (3dQwarp)         */
#define QW_HPEN_FIRST_LEV 3    /* first level to apply the penalty       */
#define QW_STOPCOST (-3.991f)  /* Pearson costs: stop if better than this */

#define QW_IS_QUINTIC(c) ((c)==QW_QUINTIC || (c)==QW_QUINTIC_LITE)

typedef struct {
    int nx, ny, nz, nxy; int64_t nxyz;          /* base grid dims */
    const float *basim, *srcim_blur;            /* base, blurred source (match images) */
    const float *wtim;                           /* weight image */
    unsigned char *bmask;                        /* base mask (weight>0) */
    double wbar;                                 /* mean nonzero weight */
    int imin, imax, jmin, jmax, kmin, kmax;      /* weight autobbox */
    qw_warp *aawarp;                             /* Haawarp: global warp being improved */
    float *aasrcim;                              /* Haasrcim: globally-warped blurred source */
    int negate;                                  /* negate the cost for the minimizer */
    qw_incor mpar;                               /* clip-param template (meth + 8 clips) */
    qw_incor incor_base;                         /* per-patch: primed with out-of-patch data */
    /* --- driver / global optimization state (AFNI H... globals) --- */
    double Hcost;                                /* current best cost (incl penalty) */
    int    Hforce;                               /* force even a light patch */
    int    fatal;                                /* set on a real allocation failure -> abort the run
                                                  * (distinct from a legitimate patch skip / NEWUOA fail) */
    float  Hfactor;                              /* max-displacement scale (1.0 here) */
    int    con;                                  /* emulated box/ball selection */
    int    ndone, nskipped;                      /* patch tallies (Hdone/Hskipped) */
    int    pen_use;                              /* Hpen_use for this level */
    double pen_fff;                              /* Hpen_fff: penalty factor this level */
    double pen_sum;                              /* Hpen_sum: out-of-patch penalty seed */
    /* --- current patch basis --- */
    int nbx, nby, nbz, nbxy;                     /* patch dims */
    float *b1d[3][3];                            /* [axis 0..2][func 0..2] 1D basis arrays */
    double dci[3];                               /* dxci/dyci/dzci scale per axis */
    int nfunc;                                   /* funcs per axis: cubic 2, quintic 3, sincc 1 */
    int npar_per_dim;                            /* params per dimension */
    int uxf[81], uyf[81], uzf[81];               /* used product decode -> axis func indices */
    float **bbb; int nbbb;                       /* npar_per_dim 3D product arrays, or NULL */
    qw_warp *hwarp, *ahwarp;                     /* patch incr H(x) and composed A(H(x)) */
    int need_ah;
    double basis_parmax;
    int npar, nparmap; int *parmap;
    float *par, *xpar, *ypar, *zpar;
    int dox, doy, doz;
    int ibot, itop, jbot, jtop, kbot, ktop;      /* inclusive patch range on the base grid */
    int nval; float *wval, *bval, *aawt;         /* patch: warped src, base, weight */
} qw_ctx;

static qw_ctx *g_qw = NULL;   /* serial-only: the NEWUOA callback reads this */

/* 1D basis func values at x in [-1,1] into f[0..nfunc-1] (HCwarp/HQwarp/HSCwarp_eval_basis). */
static void qw_basis1d(int code, float x, float *f) {
    float aa = fabsf(x);
    if (QW_IS_QUINTIC(code)) {
        if (aa >= 1.0f) { f[0] = f[1] = f[2] = 0.0f; return; }
        float bb = 1.0f - aa; bb = bb*bb*bb; float aq = aa*aa;
        f[0] = bb * ((6.0f*aq+3.0f)*aa + 1.0f);
        f[1] = bb * x * (3.0f*aa+1.0f) * 5.0625f;
        f[2] = aq * bb * 28.935f;
    } else if (code == QW_SINCC) {
        f[0] = (aa >= 1.0f) ? 0.0f : (0.5f * (1.0f + qw_sinc(aa) - qw_sinc(1.0f-aa)));
    } else {   /* cubic */
        if (aa >= 1.0f) { f[0] = f[1] = 0.0f; return; }
        float bb = 1.0f - aa; bb = bb*bb;
        f[0] = bb * (1.0f + 2.0f*aa);
        f[1] = bb * x * 6.75f;
    }
}

/* Fill the used-product list (npar_per_dim indices) + decode each into per-axis func indices.
 * Cubic products: index = zf + 2*yf + 4*xf (zf,yf,xf in {0,1}); LITE uses {0,1,2,4}.
 * Quintic products: index = zf + 3*yf + 9*xf (in {0,1,2}); LITE uses {0,1,2,3,4,6,9,10,12,18}. */
static void qw_basis_used(qw_ctx *H, int code) {
    static const int cub_lite[4]  = {0,1,2,4};
    static const int qui_lite[10] = {0,1,2,3,4,6,9,10,12,18};
    int npd; const int *lst = NULL; int quintic = QW_IS_QUINTIC(code);
    if (code == QW_CUBIC)            { npd = 8;  }
    else if (code == QW_CUBIC_LITE)  { npd = 4;  lst = cub_lite; }
    else if (code == QW_QUINTIC)     { npd = 27; }
    else if (code == QW_QUINTIC_LITE){ npd = 10; lst = qui_lite; }
    else                             { npd = 1;  }   /* sincc */
    H->npar_per_dim = npd;
    for (int m = 0; m < npd; m++) {
        int p = lst ? lst[m] : m;
        if (code == QW_SINCC) { H->uxf[m] = H->uyf[m] = H->uzf[m] = 0; }
        else if (quintic)     { H->uzf[m] = p % 3; H->uyf[m] = (p/3) % 3; H->uxf[m] = p / 9; }
        else                  { H->uzf[m] = p & 1; H->uyf[m] = (p>>1) & 1; H->uxf[m] = (p>>2) & 1; }
    }
}

/* HCwarp/HQwarp/HSCwarp_setup_basis (unified): 1D basis arrays per axis over [-1,1] via
 * COMPUTE_CAB (ILEFT=-0.5, IRGHT=n-0.5 -> cb=2/n, dci=n/2), optionally the 3D tensor-product
 * arrays (<=1GB), and fresh patch warps. An axis narrower than NGMIN gets a constant basis
 * (no displacement, dci=0), mirroring IW3D_munge_flags' NOxDEP handling. */
static void qw_setup_basis(qw_ctx *H, int code, int nx, int ny, int nz) {
    for (int a = 0; a < 3; a++) for (int f = 0; f < 3; f++) { free(H->b1d[a][f]); H->b1d[a][f] = NULL; }
    if (H->bbb) { for (int i = 0; i < H->nbbb; i++) free(H->bbb[i]); free(H->bbb); H->bbb = NULL; H->nbbb = 0; }
    qw_warp_destroy(H->hwarp); qw_warp_destroy(H->ahwarp); H->hwarp = H->ahwarp = NULL;

    H->nbx = nx; H->nby = ny; H->nbz = nz; H->nbxy = nx * ny;
    H->nfunc = QW_IS_QUINTIC(code) ? 3 : (code == QW_SINCC ? 1 : 2);
    qw_basis_used(H, code);

    int dims[3] = { nx, ny, nz };
    for (int a = 0; a < 3; a++) {
        int n = dims[a];
        int ok1d = 1;
        for (int f = 0; f < H->nfunc; f++) { H->b1d[a][f] = (float *)malloc(sizeof(float) * n); if (!H->b1d[a][f]) ok1d = 0; }
        if (!ok1d) return;   /* OOM: leave hwarp/ahwarp NULL so qw_improve_warp's -2 guard fires */
        if (n < QW_NGMIN) {   /* NOxDEP: constant basis, no displacement on this axis */
            for (int i = 0; i < n; i++) { H->b1d[a][0][i] = 1.0f; for (int f = 1; f < H->nfunc; f++) H->b1d[a][f][i] = 0.0f; }
            H->dci[a] = 0.0;
        } else {
            float cb = 2.0f / ((n - 0.5f) - (-0.5f)), ca = -1.0f - cb * (-0.5f);
            H->dci[a] = 1.0f / cb;
            for (int i = 0; i < n; i++) { float fv[3]; qw_basis1d(code, ca + i * cb, fv);
                for (int f = 0; f < H->nfunc; f++) H->b1d[a][f][i] = fv[f]; }
        }
    }

    int64_t npv = (int64_t)nx * ny * nz;
    int npd = H->npar_per_dim;
    if ((double)sizeof(float) * npd * npv / (1000.0 * 1048576.0) <= 1.0) {   /* HmaxmemG = 1GB */
        H->bbb = (float **)calloc((size_t)npd, sizeof(float *));   /* calloc: entries after a failed
                        * inner malloc stay NULL, so the cleanup free() below is always safe. */
        int ok = (H->bbb != NULL);
        for (int m = 0; ok && m < npd; m++) { H->bbb[m] = (float *)malloc(sizeof(float) * npv); if (!H->bbb[m]) ok = 0; }
        if (ok) {
            H->nbbb = npd;   /* MUST track the count — the free loops in qw_setup_basis/qw_ctx_free
                              * iterate i < nbbb; leaving it 0 leaks every patch's inner arrays. */
            int64_t qq = 0;
            for (int kk = 0; kk < nz; kk++) for (int jj = 0; jj < ny; jj++) for (int ii = 0; ii < nx; ii++, qq++)
                for (int m = 0; m < npd; m++)
                    H->bbb[m][qq] = H->b1d[2][H->uzf[m]][kk] * H->b1d[1][H->uyf[m]][jj] * H->b1d[0][H->uxf[m]][ii];
        } else { if (H->bbb) { for (int m = 0; m < npd; m++) free(H->bbb[m]); free(H->bbb); H->bbb = NULL; } }
    }
    H->hwarp  = qw_warp_create(nx, ny, nz);
    H->ahwarp = qw_warp_create(nx, ny, nz);
}

/* H?warp_eval: patch displacement at patch index qq (uses the 3D arrays if built, else the
 * on-the-fly tensor product from the 1D arrays — numerically identical). */
static inline void qw_heval(const qw_ctx *H, int qq, float *xx, float *yy, float *zz) {
    int n = H->npar_per_dim;
    const float *xp = H->xpar, *yp = H->ypar, *zp = H->zpar;
    float sx = 0.0f, sy = 0.0f, sz = 0.0f;
    if (H->bbb) {
        for (int m = 0; m < n; m++) { float b = H->bbb[m][qq]; sx += b*xp[m]; sy += b*yp[m]; sz += b*zp[m]; }
    } else {
        int ii = qq % H->nbx, kk = qq / H->nbxy, jj = (qq - kk * H->nbxy) / H->nbx;
        for (int m = 0; m < n; m++) {
            float b = H->b1d[2][H->uzf[m]][kk] * H->b1d[1][H->uyf[m]][jj] * H->b1d[0][H->uxf[m]][ii];
            sx += b*xp[m]; sy += b*yp[m]; sz += b*zp[m];
        }
    }
    *xx = H->dox ? (float)(H->dci[0] * sx) : 0.0f;
    *yy = H->doy ? (float)(H->dci[1] * sy) : 0.0f;
    *zz = H->doz ? (float)(H->dci[2] * sz) : 0.0f;
}

/* Hwarp_apply: for each patch voxel, evaluate the patch displacement, linearly interpolate the
 * GLOBAL warp Haawarp at the patch-warped location to compose (-> AHwarp), then linearly sample
 * the SOURCE at the composed location -> val. (Faithful to mri_nwarp.c Hwarp_apply: compose
 * interp is top-clamp-only QLIP since patches never touch the edge; source interp is full CLIP
 * after clamping loc to [-0.499, n-0.501].) */
static void qw_hwarp_apply(qw_ctx *H, float *val) {
    int nbx = H->nbx, nby = H->nby, nbz = H->nbz, nbxy = nbx * nby;
    int64_t nbxyz = (int64_t)nbxy * nbz;
    int nAx = H->aawarp->nx, nAy = H->aawarp->ny, nAz = H->aawarp->nz, nAxy = nAx * nAy;
    int nAx1 = nAx-1, nAy1 = nAy-1, nAz1 = nAz-1;
    float nAxh = nAx-0.501f, nAyh = nAy-0.501f, nAzh = nAz-0.501f;
    const float *Axd = H->aawarp->xd, *Ayd = H->aawarp->yd, *Azd = H->aawarp->zd;
    float *hxd = H->hwarp->xd, *hyd = H->hwarp->yd, *hzd = H->hwarp->zd;
    float *bxd = H->ahwarp->xd, *byd = H->ahwarp->yd, *bzd = H->ahwarp->zd;
    const float *sar = H->srcim_blur;   /* SRCIM = the ORIGINAL blurred source (sampled at the
                                         * composed AHwarp location; NOT the already-warped
                                         * aasrcim — sampling that would double-apply the warp) */
    unsigned char *bmask = H->bmask;
    int need_ah = H->need_ah;
#define QW_IJK(i,j,k) ((i) + (j)*nAx + (int64_t)(k)*nAxy)
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) if(nbxyz > 9999)
#endif
    for (int64_t qq = 0; qq < nbxyz; qq++) {
        int ii = qq % nbx, kk = qq / nbxy, jj = (qq - (int64_t)kk*nbxy) / nbx;
        int gi = ii + H->ibot, gj = jj + H->jbot, gk = kk + H->kbot;
        int need_val = bmask[QW_IJK(gi, gj, gk)] != 0;
        if (!need_val && !need_ah) { val[qq] = 0.0f; continue; }
        float hx, hy, hz; qw_heval(H, (int)qq, &hx, &hy, &hz);
        hxd[qq] = hx; hyd[qq] = hy; hzd[qq] = hz;
        /* --- compose: interp Haawarp displacement at (gi+hx, gj+hy, gk+hz) --- */
        float xq = gi + hx; int ix = (int)xq; float fx = xq - ix;
        float yq = gj + hy; int jy = (int)yq; float fy = yq - jy;
        float zq = gk + hz; int kz = (int)zq; float fz = zq - kz;
        /* AFNI QLIP only top-clamps (patch bottoms are floored at ibbb>=1 and the basis is ~0
         * at patch edges, so the composed location never reaches index 0). The lower clamp is a
         * defensive no-op that keeps the gather in-bounds if those patch-bound invariants ever
         * change — it never fires here, so output stays bit-identical. */
        int i0 = ix, i1 = ix+1; if (i0 < 0) i0 = 0; else if (i0 > nAx1) i0 = nAx1; if (i1 < 0) i1 = 0; else if (i1 > nAx1) i1 = nAx1;
        int j0 = jy, j1 = jy+1; if (j0 < 0) j0 = 0; else if (j0 > nAy1) j0 = nAy1; if (j1 < 0) j1 = 0; else if (j1 > nAy1) j1 = nAy1;
        int k0 = kz, k1 = kz+1; if (k0 < 0) k0 = 0; else if (k0 > nAz1) k0 = nAz1; if (k1 < 0) k1 = 0; else if (k1 > nAz1) k1 = nAz1;
        float w00 = 1.0f-fx, w1 = fx;
#define QW_AX(arr,j,k) (w00*arr[QW_IJK(i0,j,k)] + w1*arr[QW_IJK(i1,j,k)])
        float fx00 = QW_AX(Axd,j0,k0), fx10 = QW_AX(Axd,j1,k0), fx01 = QW_AX(Axd,j0,k1), fx11 = QW_AX(Axd,j1,k1);
        float gx00 = QW_AX(Ayd,j0,k0), gx10 = QW_AX(Ayd,j1,k0), gx01 = QW_AX(Ayd,j0,k1), gx11 = QW_AX(Ayd,j1,k1);
        float hx00 = QW_AX(Azd,j0,k0), hx10 = QW_AX(Azd,j1,k0), hx01 = QW_AX(Azd,j0,k1), hx11 = QW_AX(Azd,j1,k1);
#undef QW_AX
        float wy = 1.0f-fy;
        float fk0 = wy*fx00 + fy*fx10, fk1 = wy*fx01 + fy*fx11;
        float gk0 = wy*gx00 + fy*gx10, gk1 = wy*gx01 + fy*gx11;
        float hk0 = wy*hx00 + fy*hx10, hk1 = wy*hx01 + fy*hx11;
        float wz = 1.0f-fz;
        bxd[qq] = wz*fk0 + fz*fk1 + hx;    /* AHwarp = interp(Haawarp) + Hwarp */
        byd[qq] = wz*gk0 + fz*gk1 + hy;
        bzd[qq] = wz*hk0 + fz*hk1 + hz;
        if (!need_val) { val[qq] = 0.0f; continue; }
        /* --- sample source at the composed absolute location --- */
        xq = bxd[qq] + gi; yq = byd[qq] + gj; zq = bzd[qq] + gk;
        if (xq < -0.499f) xq = -0.499f; else if (xq > nAxh) xq = nAxh;
        if (yq < -0.499f) yq = -0.499f; else if (yq > nAyh) yq = nAyh;
        if (zq < -0.499f) zq = -0.499f; else if (zq > nAzh) zq = nAzh;
        ix = (int)floorf(xq); fx = xq - ix; jy = (int)floorf(yq); fy = yq - jy; kz = (int)floorf(zq); fz = zq - kz;
        i0 = ix; i1 = ix+1; if (i0 < 0) i0 = 0; else if (i0 > nAx1) i0 = nAx1; if (i1 < 0) i1 = 0; else if (i1 > nAx1) i1 = nAx1;
        j0 = jy; j1 = jy+1; if (j0 < 0) j0 = 0; else if (j0 > nAy1) j0 = nAy1; if (j1 < 0) j1 = 0; else if (j1 > nAy1) j1 = nAy1;
        k0 = kz; k1 = kz+1; if (k0 < 0) k0 = 0; else if (k0 > nAz1) k0 = nAz1; if (k1 < 0) k1 = 0; else if (k1 > nAz1) k1 = nAz1;
        w00 = 1.0f-fx; w1 = fx;
#define QW_SX(j,k) (w00*sar[QW_IJK(i0,j,k)] + w1*sar[QW_IJK(i1,j,k)])
        float s00 = QW_SX(j0,k0), s10 = QW_SX(j1,k0), s01 = QW_SX(j0,k1), s11 = QW_SX(j1,k1);
#undef QW_SX
        wy = 1.0f-fy;
        float sk0 = wy*s00 + fy*s10, sk1 = wy*s01 + fy*s11;
        val[qq] = (1.0f-fz)*sk0 + fz*sk1;
    }
#undef QW_IJK
}

/* IW3D_scalar_costfun: the NEWUOA callback. Scatter params -> H->par, warp, INCOR-evaluate,
 * negate for the minimizer, add the deformation penalty when enabled. Reads process-global g_qw. */
static double qw_scalar_costfun(int npar, double *dpar) {
    qw_ctx *H = g_qw;
    if (H->parmap) { for (int i = 0; i < H->npar; i++) H->par[i] = 0.0f;
                     for (int i = 0; i < npar; i++) H->par[H->parmap[i]] = (float)dpar[i]; }
    else           { for (int i = 0; i < H->npar; i++) H->par[i] = (float)dpar[i]; }
    qw_hwarp_apply(H, H->wval);
    double cost = qw_incor_eval(&H->incor_base, H->nval, H->bval, H->wval, H->aawt);
    if (H->negate) cost = -cost;
    /* NaN guard: map to a cost WORSE than any real one (post-negate cost is atanh(r) in [-4,4],
     * minimized) so the optimizer never prefers a NaN evaluation. Inputs are finite-checked at
     * entry, so this should not fire — it is pure defense and leaves the normal path unchanged. */
    if (!(cost == cost)) cost = 1.0e9;
    if (H->pen_use) {   /* HPEN_penalty: Hpen_sum (out-of-patch) + energy(AHwarp), scaled */
        double hsum = H->pen_sum + qw_load_energy(H->ahwarp);
        if (!H->ahwarp->je || !H->ahwarp->se) H->fatal = 1;   /* energy OOM: don't silently drop the penalty */
        if (hsum > 0.0) hsum = H->pen_fff * pow(hsum, 0.25);
        cost += hsum;
    }
    return cost;
}

extern int powell_newuoa_con(int ndim, double *x, double *xbot, double *xtop,
                             int nrand, double rstart, double rend, int maxcall,
                             double (*cost)(int, double *));
extern void powell_set_mfac(float mm, float aa);
extern void powell_get_mfac(float *mm, float *aa);

/* IW3D_improve_warp: optimize one patch with the given basis, compositing the result into
 * aawarp/aasrcim (and the je/se energy fields). Returns NEWUOA iterations, or <= 0 on skip/fail.
 * Faithful to mri_nwarp.c: light-patch and too-few-points skips, INCOR primed from the
 * out-of-patch data once per patch, optional penalty seed, final composite. */
static int qw_improve_warp(qw_ctx *H, int code, int ibot, int itop, int jbot, int jtop, int kbot, int ktop) {
    if (H->fatal) return -4;   /* a prior allocation failure aborted the run — do no more work */
    if (ibot < 0) ibot = 0; if (itop > H->nx-1) itop = H->nx-1;
    if (jbot < 0) jbot = 0; if (jtop > H->ny-1) jtop = H->ny-1;
    if (kbot < 0) kbot = 0; if (ktop > H->nz-1) ktop = H->nz-1;
    int nxh = itop-ibot+1, nyh = jtop-jbot+1, nzh = ktop-kbot+1;
    if (nxh < QW_NGMIN && nyh < QW_NGMIN && nzh < QW_NGMIN) { H->nskipped++; return 0; }
    H->ibot = ibot; H->itop = itop; H->jbot = jbot; H->jtop = jtop; H->kbot = kbot; H->ktop = ktop;
    int nval = nxh * nyh * nzh;

    /* ball constraint only for tiny patches (emulated: selects parmax, optimizer stays box) */
    int ballopt = (H->con == QW_SC_BALL);
    int nball = QW_IS_QUINTIC(code) ? 19 : 15;
    if (nxh < nball || nyh < nball || nzh < nball) { ballopt = 1; H->con = QW_SC_BALL; }

    /* is this region weighty enough to bother with? */
    double wsum = 0.0; int nwb = 0;
    for (int kk = kbot; kk <= ktop; kk++) for (int jj = jbot; jj <= jtop; jj++) for (int ii = ibot; ii <= itop; ii++) {
        int64_t q = ii + (int64_t)jj*H->nx + (int64_t)kk*H->nxy;
        if (H->bmask[q]) { wsum += H->wtim[q]; nwb++; }
    }
    if (!H->Hforce && (nwb < 0.333f*nval || wsum < 0.166f*nval*H->wbar)) { H->nskipped++; return 0; }

    /* basis + max-displacement scale (box vs ball value; * Hfactor) */
    double prad = 0.333, boxv, ballv; int npar;
    switch (code) {
        case QW_CUBIC:        boxv = 0.0280; ballv = 0.0790; npar = 24; break;
        case QW_QUINTIC:      boxv = 0.0099; ballv = 0.0611; npar = 81; break;
        case QW_QUINTIC_LITE: boxv = 0.0267; ballv = 0.1111; npar = 30; break;
        case QW_SINCC:        boxv = 0.1666; ballv = 0.2345; npar = 3;  break;
        default:              code = QW_CUBIC_LITE;
        case QW_CUBIC_LITE:   boxv = 0.0421; ballv = 0.1141; npar = 12; break;
    }
    H->basis_parmax = (ballopt ? ballv : boxv) * H->Hfactor;
    H->npar = npar; H->nparmap = npar; free(H->parmap); H->parmap = NULL;
    qw_setup_basis(H, code, nxh, nyh, nzh);
    if (!H->hwarp || !H->ahwarp) { H->fatal = 1; return -2; }   /* basis/warp OOM -> abort the run */

    if (nwb < 5 * H->nparmap) { H->nskipped++; return 0; }
    H->dox = H->doy = H->doz = 1;

    /* grow the reusable scratch buffers (realloc via a temp so an OOM failure doesn't leak the
     * original block — the ctx still owns the old pointer and frees it in qw_ctx_free). */
    float *tp;
    tp = (float *)realloc(H->par,  sizeof(float) * H->npar); if (!tp) { H->fatal = 1; return -4; } H->par  = tp;
    tp = (float *)realloc(H->wval, sizeof(float) * nval);    if (!tp) { H->fatal = 1; return -4; } H->wval = tp;
    tp = (float *)realloc(H->bval, sizeof(float) * nval);    if (!tp) { H->fatal = 1; return -4; } H->bval = tp;
    tp = (float *)realloc(H->aawt, sizeof(float) * nval);    if (!tp) { H->fatal = 1; return -4; } H->aawt = tp;
    for (int i = 0; i < H->npar; i++) H->par[i] = 0.0f;
    int npd = H->npar_per_dim;
    H->xpar = H->par; H->ypar = H->par + npd; H->zpar = H->par + 2*npd;
    H->nval = nval;

    /* prime INCOR with the OUT-of-patch data: working copy of the weight, extract patch
     * base/weight, zero the patch weight, accumulate the whole volume (in-patch weight 0). */
    float *wbfar = (float *)malloc(sizeof(float) * (size_t)H->nxyz);
    if (!wbfar) { H->fatal = 1; return -4; }
    memcpy(wbfar, H->wtim, sizeof(float) * (size_t)H->nxyz);
    int pp = 0; int base_const = 1; float b0v = 0.0f;
    for (int kk = kbot; kk <= ktop; kk++) for (int jj = jbot; jj <= jtop; jj++) for (int ii = ibot; ii <= itop; ii++, pp++) {
        int64_t q = ii + (int64_t)jj*H->nx + (int64_t)kk*H->nxy;
        H->bval[pp] = H->basim[q];
        H->aawt[pp] = wbfar[q];
        wbfar[q] = 0.0f;
        if (pp == 0) b0v = H->bval[0]; else if (H->bval[pp] != b0v) base_const = 0;
    }
    if (base_const) { free(wbfar); H->nskipped++; return 0; }   /* can't correlate a constant base */
    H->incor_base = H->mpar;
    qw_incor_addto(&H->incor_base, (int)H->nxyz, H->basim, H->aasrcim, wbfar);
    free(wbfar);

    /* penalty seed: zero the patch je/se in the global energy field, sum the rest */
    H->need_ah = H->pen_use;
    if (H->pen_use) {
        float *je = H->aawarp->je, *se = H->aawarp->se;
        if (je && se) {
            for (int kk = kbot; kk <= ktop; kk++) for (int jj = jbot; jj <= jtop; jj++) for (int ii = ibot; ii <= itop; ii++) {
                int64_t q = ii + (int64_t)jj*H->nx + (int64_t)kk*H->nxy; je[q] = se[q] = 0.0f;
            }
            H->pen_sum = qw_hpen_addup(H->nxyz, je, se);
        } else H->pen_sum = 0.0;
    }

    /* optimize (box-constrained NEWUOA; see file-header DEVIATION note on box/ball) */
    double *parvec = (double *)calloc(H->nparmap, sizeof(double));
    double *xbot = (double *)malloc(sizeof(double) * H->nparmap);
    double *xtop = (double *)malloc(sizeof(double) * H->nparmap);
    if (!parvec || !xbot || !xtop) { free(parvec); free(xbot); free(xtop); H->fatal = 1; return -4; }
    for (int i = 0; i < H->nparmap; i++) { parvec[i] = 0.0; xbot[i] = -H->basis_parmax; xtop[i] = H->basis_parmax; }
    /* NEWUOA's sampling factors are thread-local shared optimizer state. qwarp needs AFNI's
     * lean 1*n+2 policy, but an embedding caller must not inherit it after this patch. Likewise,
     * restore the callback context even though the documented API is serial/non-reentrant. */
    float saved_mfac, saved_afac;
    qw_ctx *saved_qw = g_qw;
    powell_get_mfac(&saved_mfac, &saved_afac);
    powell_set_mfac(1.001f, 2.001f);
    int itmax = 8 * H->nparmap + 31;
    g_qw = H;
    int iter = powell_newuoa_con(H->nparmap, parvec, xbot, xtop, 0, prad, 0.009 * prad, itmax, qw_scalar_costfun);
    if (iter == -1)   /* the one recoverable code: retry once with a larger end radius */
        iter = powell_newuoa_con(H->nparmap, parvec, xbot, xtop, 0, prad, 0.09 * prad, itmax, qw_scalar_costfun);
    powell_set_mfac(saved_mfac, saved_afac);
    /* After the retry, a NEGATIVE result is a NEWUOA hard failure (workspace alloc -7, bad args
     * -2..-6) — treat as fatal so an optimizer OOM cannot yield a silently-degraded warp. Zero is
     * a benign non-convergence / no-improvement result and stays a patch skip. */
    if (iter < 0)  { g_qw = saved_qw; free(parvec); free(xbot); free(xtop); H->fatal = 1; return -4; }
    if (iter == 0) { g_qw = saved_qw; free(parvec); free(xbot); free(xtop); H->nskipped++; return 0; }

    /* final eval at the optimum: fills wval + ahwarp displacements (and cost incl penalty) */
    H->need_ah = 1;
    H->Hcost = qw_scalar_costfun(H->nparmap, parvec);
    g_qw = saved_qw;
    (void)qw_load_energy(H->ahwarp);   /* fill ahwarp je/se for the composite */
    /* ahwarp's je/se are lazily allocated PER PATCH; a partial OOM (e.g. je ok, se NULL) would
     * NULL-deref in the composite below. Require BOTH before proceeding, else abort the run. */
    if (!H->ahwarp->je || !H->ahwarp->se) { free(parvec); free(xbot); free(xtop); H->fatal = 1; return -4; }
    /* composite AHwarp/Hwval back into the globals (displacements, warped source, energies) */
    float *sar = H->aasrcim;
    float *Axd = H->aawarp->xd, *Ayd = H->aawarp->yd, *Azd = H->aawarp->zd;
    float *Aje = H->aawarp->je, *Ase = H->aawarp->se;
    const float *bxd = H->ahwarp->xd, *byd = H->ahwarp->yd, *bzd = H->ahwarp->zd;
    const float *bje = H->ahwarp->je, *bse = H->ahwarp->se;
    pp = 0;
    for (int kk = kbot; kk <= ktop; kk++) for (int jj = jbot; jj <= jtop; jj++) for (int ii = ibot; ii <= itop; ii++, pp++) {
        int64_t q = ii + (int64_t)jj*H->nx + (int64_t)kk*H->nxy;
        sar[q] = H->wval[pp];
        Axd[q] = bxd[pp]; Ayd[q] = byd[pp]; Azd[q] = bzd[pp];
        if (Aje && Ase && bje && bse) { Aje[q] = bje[pp]; Ase[q] = bse[pp]; }
    }
    free(parvec); free(xbot); free(xtop);
    H->ndone++;
    return iter;
}

/* IW3D_setup_for_improvement (no-penalty subset): fill the context from base/blurred-source/
 * weight, build the mask + INCOR clip params (from the blurred source, = AFNI SRCIM), seed
 * aawarp=identity, aasrcim=blurred source, and load the initial energy field. Returns 0 on OK. */
static int qw_setup_for_improvement(qw_ctx *H, const float *base, const float *src_blur,
                                    const float *wt, int nx, int ny, int nz) {
    memset(H, 0, sizeof *H);
    H->nx = nx; H->ny = ny; H->nz = nz; H->nxy = nx * ny; H->nxyz = (int64_t)nx * ny * nz;
    H->basim = base; H->srcim_blur = src_blur; H->wtim = wt;
    H->con = QW_SC_BOX; H->Hfactor = 1.0f;
    H->bmask = (unsigned char *)malloc((size_t)H->nxyz);
    if (!H->bmask) return 1;
    double wsum = 0.0; int64_t nwb = 0;
    for (int64_t i = 0; i < H->nxyz; i++) { H->bmask[i] = (wt[i] > 0.0f); if (H->bmask[i]) { wsum += wt[i]; nwb++; } }
    H->wbar = (nwb > 0) ? wsum / nwb : 0.0;
    /* INCOR clip params from the masked base(x) vs blurred source(y) samples */
    float *xar = (float *)malloc(sizeof(float) * (size_t)(nwb ? nwb : 1));
    float *yar = (float *)malloc(sizeof(float) * (size_t)(nwb ? nwb : 1));
    if (!xar || !yar) { free(xar); free(yar); return 1; }
    int64_t k = 0;
    for (int64_t i = 0; i < H->nxyz; i++) if (H->bmask[i]) { xar[k] = base[i]; yar[k] = src_blur[i]; k++; }
    memset(&H->mpar, 0, sizeof H->mpar);
    H->mpar.meth = QW_INCOR_PEARCLP;
    if (qw_incor_setup_clip(&H->mpar, xar, yar, k) != 0)     /* degenerate -> plain Pearson */
        H->mpar.meth = QW_INCOR_PEARSON;
    free(xar); free(yar);
    H->negate = 1;
    H->aawarp = qw_warp_create(nx, ny, nz);          /* identity */
    H->aasrcim = (float *)malloc(sizeof(float) * (size_t)H->nxyz);
    if (!H->aawarp || !H->aasrcim) return 1;
    memcpy(H->aasrcim, src_blur, sizeof(float) * (size_t)H->nxyz);   /* identity-warped blurred src */
    (void)qw_load_energy(H->aawarp);   /* allocate + initialize je/se for the penalty */
    if (!H->aawarp->je || !H->aawarp->se) return 1;   /* je/se OOM: the penalty would silently be 0 */
    return 0;
}

static void qw_ctx_free(qw_ctx *H) {
    free(H->bmask); qw_warp_destroy(H->aawarp); free(H->aasrcim);
    qw_warp_destroy(H->hwarp); qw_warp_destroy(H->ahwarp);
    for (int a = 0; a < 3; a++) for (int f = 0; f < 3; f++) free(H->b1d[a][f]);
    if (H->bbb) { for (int i = 0; i < H->nbbb; i++) free(H->bbb[i]); free(H->bbb); }
    free(H->parmap); free(H->par); free(H->wval); free(H->bval); free(H->aawt);
}

/* IW3D_warpomatic: the level-driven optimization. Level 0 is one global patch with a fixed
 * basis escalation; then patches shrink by QW_SHRINK each level (forced odd, ~50% overlap,
 * alternating sweep parity) until they reach the minimum patch size. The penalty turns on at
 * level >= QW_HPEN_FIRST_LEV. On return the optimized warp is in H->aawarp. Faithful port of
 * mri_nwarp.c IW3D_warpomatic for the default 3dQwarp configuration (cubic-lite ordinary
 * levels, quintic-lite at level 0, penalty on, box-emulated constraints). */
static void qw_warpomatic(qw_ctx *H) {
    int Hnx = H->nx, Hny = H->ny, Hnz = H->nz;
    int imin = H->imin, imax = H->imax, jmin = H->jmin, jmax = H->jmax, kmin = H->kmin, kmax = H->kmax;
    int verb = (getenv("QW_VERB") != NULL);   /* zero release overhead; diagnostics only */
    H->ndone = H->nskipped = 0; H->Hforce = 0; H->Hfactor = 1.0f;

    /* ---- level 0: one global patch, basis escalation, no penalty ---- */
    int xwid = (imax-imin)/8, ywid = (jmax-jmin)/8, zwid = (kmax-kmin)/8;
    int ibbb = (1 > imin-xwid) ? 1 : imin-xwid, jbbb = (1 > jmin-ywid) ? 1 : jmin-ywid, kbbb = (1 > kmin-zwid) ? 1 : kmin-zwid;
    int ittt = (Hnx-2 < imax+xwid) ? Hnx-2 : imax+xwid, jttt = (Hny-2 < jmax+ywid) ? Hny-2 : jmax+ywid, kttt = (Hnz-2 < kmax+zwid) ? Hnz-2 : kmax+zwid;
    int xwid0 = ittt-ibbb+1, ywid0 = jttt-jbbb+1, zwid0 = kttt-kbbb+1;

    H->Hforce = 1; H->pen_use = 0;
    H->con = QW_SC_BOX;
    (void)qw_improve_warp(H, QW_SINCC,        ibbb, ittt, jbbb, jttt, kbbb, kttt);
    (void)qw_improve_warp(H, QW_CUBIC_LITE,   ibbb, ittt, jbbb, jttt, kbbb, kttt);
    H->con = QW_SC_BALL;
    (void)qw_improve_warp(H, QW_CUBIC,        ibbb, ittt, jbbb, jttt, kbbb, kttt);
    H->con = QW_SC_BALL;
    (void)qw_improve_warp(H, QW_QUINTIC_LITE, ibbb, ittt, jbbb, jttt, kbbb, kttt);
    H->con = QW_SC_BALL;
    (void)qw_improve_warp(H, QW_QUINTIC,      ibbb, ittt, jbbb, jttt, kbbb, kttt);
    if (H->fatal) return;   /* allocation failure at level 0 -> abort (qwarp_run writes no output) */
    if (verb) fprintf(stderr, "[qw] lev0 patch %d..%d %d..%d %d..%d done=%d skip=%d Hcost=%.5f\n",
                      ibbb, ittt, jbbb, jttt, kbbb, kttt, H->ndone, H->nskipped, H->Hcost);

    H->Hforce = 0; H->con = QW_SC_BOX;

    int ngmin = QW_HNGMIN; if (ngmin < QW_NGMIN) ngmin = QW_NGMIN; else if (ngmin % 2 == 0) ngmin--;
    if (ngmin >= Hnx && ngmin >= Hny && ngmin >= Hnz) return;

    /* ---- iterate down to finer and finer patches ---- */
    int levs = 1, leve = 666, levdone = 0;
    for (int lev = levs; lev <= leve && !levdone; lev++) {
        float flev = powf((float)(lev - levs + 1), 0.333f);
        H->pen_fff = QW_HPEN_FBASE * (flev < 3.21f ? flev : 3.21f);
        H->pen_use = (H->pen_fff > 0.0) && (lev >= QW_HPEN_FIRST_LEV);
        if (lev == QW_HPEN_FIRST_LEV) H->pen_fff *= 0.5;

        flev = powf(QW_SHRINK, (float)lev);
        xwid = (int)(xwid0*flev); if (xwid % 2 == 0) xwid++;
        ywid = (int)(ywid0*flev); if (ywid % 2 == 0) ywid++;
        zwid = (int)(zwid0*flev); if (zwid % 2 == 0) zwid++;

        int dox = (xwid >= ngmin), doy = (ywid >= ngmin), doz = (zwid >= ngmin);
        if (!dox && !doy && !doz) break;
        if (xwid < ngmin) xwid = (Hnx < ngmin) ? Hnx : ngmin;
        if (ywid < ngmin) ywid = (Hny < ngmin) ? Hny : ngmin;
        if (zwid < ngmin) zwid = (Hnz < ngmin) ? Hnz : ngmin;

        /* jump straight to the minimum patch if the next shrink would undershoot it */
        float g1 = xwid/(float)ngmin, g2 = ywid/(float)ngmin, g3 = zwid/(float)ngmin;
        if (g2 > g1) g1 = g2; if (g3 > g1) g1 = g3;
        if (g1 > 1.0f && g1*QW_SHRINK <= 1.00001f) {
            if (xwid > ngmin) xwid = ngmin; if (ywid > ngmin) ywid = ngmin; if (zwid > ngmin) zwid = ngmin;
            levdone = 1;
        } else {
            int mx = (xwid > ywid) ? xwid : ywid; if (zwid > mx) mx = zwid; levdone = (mx == ngmin);
        }

        int xdel = (xwid-1)/2; if (xdel == 0) xdel = 1;
        int ydel = (ywid-1)/2; if (ydel == 0) ydel = 1;
        int zdel = (zwid-1)/2; if (zdel == 0) zdel = 1;
        int diii = xdel, djjj = ydel, dkkk = zdel;

        ibbb = imin-xdel/4-1; if (ibbb <= 0) ibbb = 1;
        jbbb = jmin-ydel/4-1; if (jbbb <= 0) jbbb = 1;
        kbbb = kmin-zdel/4-1; if (kbbb <= 0) kbbb = 1;
        ittt = imax+xdel/4+1; if (ittt >= Hnx-1) ittt = Hnx-2;
        jttt = jmax+ydel/4+1; if (jttt >= Hny-1) jttt = Hny-2;
        kttt = kmax+zdel/4+1; if (kttt >= Hnz-1) kttt = Hnz-2;

        (void)qw_load_energy(H->aawarp);   /* refresh energy field for this level's penalty */
        {   /* re-warp the source from scratch to reset compose-interpolation drift (AFNI does
             * this each level: Haasrcim = IW3D_warp_floatim(Haawarp, SRCIM, LINEAR)) */
            float *fresh = qw_warp_apply(H->aawarp, H->srcim_blur, Hnx, Hny, Hnz, QW_LINEAR, 1.0f);
            if (!fresh) { H->fatal = 1; return; }   /* OOM: don't silently continue with a stale reslice */
            free(H->aasrcim); H->aasrcim = fresh;
        }

        H->ndone = 0; H->nskipped = 0;
        int lev_odd = (lev % 2 == 1);

        if (lev_odd) {   /* bottom-to-top sweep, ijk order */
            for (int kdon = 0, kbot = kbbb; !kdon; kbot += dkkk) {
                int ktop = kbot+zwid-1;
                if (ktop >= kttt)          { ktop = kttt; kbot = ktop+1-zwid; kdon = 1; }
                else if (ktop >= kttt-zwid/4) { ktop = kttt; kdon = 1; }
                for (int jdon = 0, jbot = jbbb; !jdon; jbot += djjj) {
                    int jtop = jbot+ywid-1;
                    if (jtop >= jttt)          { jtop = jttt; jbot = jtop+1-ywid; jdon = 1; }
                    else if (jtop >= jttt-ywid/4) { jtop = jttt; jdon = 1; }
                    for (int idon = 0, ibot = ibbb; !idon; ibot += diii) {
                        int itop = ibot+xwid-1;
                        if (itop >= ittt)          { itop = ittt; ibot = itop+1-xwid; idon = 1; }
                        else if (itop >= ittt-xwid/4) { itop = ittt; idon = 1; }
                        (void)qw_improve_warp(H, QW_CUBIC_LITE, ibot, itop, jbot, jtop, kbot, ktop);
                        if (H->fatal) return;   /* allocation failure -> abort */
                        if (H->Hcost < QW_STOPCOST) return;
                    }
                }
            }
        } else {   /* top-to-bottom sweep, kji order */
            for (int idon = 0, itop = ittt; !idon; itop -= diii) {
                int ibot = itop+1-xwid;
                if (ibot <= ibbb)          { ibot = ibbb; itop = ibot+xwid-1; idon = 1; }
                else if (ibot <= ibbb+xwid/4) { ibot = ibbb; idon = 1; }
                for (int jdon = 0, jtop = jttt; !jdon; jtop -= djjj) {
                    int jbot = jtop+1-ywid;
                    if (jbot <= jbbb)          { jbot = jbbb; jtop = jbot+ywid-1; jdon = 1; }
                    else if (jbot <= jbbb+ywid/4) { jbot = jbbb; jdon = 1; }
                    for (int kdon = 0, ktop = kttt; !kdon; ktop -= dkkk) {
                        int kbot = ktop+1-zwid;
                        if (kbot <= kbbb)          { kbot = kbbb; ktop = kbot+zwid-1; kdon = 1; }
                        else if (kbot <= kbbb+zwid/4) { kbot = kbbb; kdon = 1; }
                        (void)qw_improve_warp(H, QW_CUBIC_LITE, ibot, itop, jbot, jtop, kbot, ktop);
                        if (H->fatal) return;   /* allocation failure -> abort */
                        if (H->Hcost < QW_STOPCOST) return;
                    }
                }
            }
        }

        /* if nothing was optimized this level, force one patch centered on the autobox */
        if (H->ndone == 0) {
            int ibot = (imin+imax-xwid)/2; if (ibot < 0) ibot = 0;
            int jbot = (jmin+jmax-ywid)/2; if (jbot < 0) jbot = 0;
            int kbot = (kmin+kmax-zwid)/2; if (kbot < 0) kbot = 0;
            int itop = ibot+xwid-1; if (itop >= Hnx) itop = Hnx-1;
            int jtop = jbot+ywid-1; if (jtop >= Hny) jtop = Hny-1;
            int ktop = kbot+zwid-1; if (ktop >= Hnz) ktop = Hnz-1;
            H->Hforce = 1;
            (void)qw_improve_warp(H, QW_CUBIC_LITE, ibot, itop, jbot, jtop, kbot, ktop);
            H->Hforce = 0;
        }
        if (verb) fprintf(stderr, "[qw] lev=%d patch=%dx%dx%d pen=%s(%.4f) done=%d skip=%d Hcost=%.5f\n",
                          lev, xwid, ywid, zwid, H->pen_use ? "on" : "off", H->pen_fff,
                          H->ndone, H->nskipped, H->Hcost);
    }
}

#ifdef QWARP_CHECKPOINT
/* Developer-only: dump a float volume on `ref`'s grid to /tmp/qwchk_<tag>.nii.gz for
 * numerical comparison against captured AFNI intermediates. Compiled out of release. */
static void qw_checkpoint(const nifti_image *ref, const float *data,
                          int nx, int ny, int nz, const char *tag) {
    nifti_image *im = (nifti_image *)malloc(sizeof *im); if (!im) return;
    *im = *ref;
    im->fname = NULL; im->iname = NULL; im->num_ext = 0; im->ext_list = NULL;
    im->ndim = 3; im->dim[0] = 3;
    im->nx = nx; im->ny = ny; im->nz = nz; im->dim[1] = nx; im->dim[2] = ny; im->dim[3] = nz;
    im->nt = im->nu = im->nv = im->nw = 1; im->dim[4] = im->dim[5] = im->dim[6] = im->dim[7] = 1;
    im->nvox = (int64_t)nx * ny * nz; im->datatype = DT_FLOAT32; im->nbyper = 4;
    im->scl_slope = 0; im->scl_inter = 0;
    size_t nb = (size_t)im->nvox * 4; im->data = malloc(nb);
    if (im->data) {
        memcpy(im->data, data, nb);
        char path[600]; snprintf(path, sizeof path, "/tmp/qwchk_%s", tag);
        nifti_set_filenames(im, path, 0, 0); nifti_image_write(im);
    }
    nifti_image_free(im);
}
#define QW_CHECKPOINT(ref, data, nx, ny, nz, tag) qw_checkpoint(ref, data, nx, ny, nz, tag)
#else
#define QW_CHECKPOINT(ref, data, nx, ny, nz, tag) ((void)0)
#endif

/*--- public entry ---------------------------------------------------------*/

int qwarp_run(const nifti_image *moving, const nifti_image *stationary,
              nifti_image **result) {
    if (result) *result = NULL;
    if (!result) return 1;
    if (qw_grid_compat(moving, stationary)) return 1;

    /* ===================== the full ported qwarp pipeline =====================
     * extract -> data-dependent zeropad -> blur source (FWHM 3 vox) -> weightize
     * the PADDED base -> warpomatic (level-0 basis escalation + shrinking patches
     * + penalty) -> single WSINC5 reslice of the UNBLURRED padded source through
     * the optimized warp -> crop back to the stationary grid. Failure is atomic:
     * *result is written only after the complete warp + reslice + crop succeed. */
    int onx = (int)stationary->nx, ony = (int)stationary->ny, onz = (int)stationary->nz;
    int rc = 1;
    float *base = qw_extract_float(stationary);   /* base = stationary (unpadded) */
    float *src  = qw_extract_float(moving);        /* source = moving  (unpadded) */
    float *basep = NULL, *srcp = NULL, *srcbp = NULL, *wtp = NULL, *warpedp = NULL, *cropped = NULL;
    qw_ctx H; int have_H = 0;
    if (!base || !src) goto pipe_done;
    /* Input policy (decision 3): reject non-finite voxels up front rather than let NaN/Inf
     * propagate through blur/weight/correlation/interp into a plausible-looking bad output. */
    int64_t onv = (int64_t)onx * ony * onz;
    if (!qw_all_finite(base, onv) || !qw_all_finite(src, onv)) {
        fprintf(stderr, "qwarp: input contains non-finite (NaN/Inf) voxels\n"); goto pipe_done; }

    int xm, xp, ym, yp, zm, zp;
    qw_compute_pad(base, onx, ony, onz, &xm, &xp, &ym, &yp, &zm, &zp);
    int pnx = onx + xm + xp, pny = ony + ym + yp, pnz = onz + zm + zp;
    /* Several inner routines take `int` voxel counts (interpolators, INCOR); a padded volume
     * past INT_MAX would truncate to a bogus count and silently produce wrong output. Reject it
     * BEFORE the padding allocations (mirrors the fast engine's nvox guard). ~2.1 Gvoxel is far
     * beyond any real brain grid. */
    int64_t pnv = (int64_t)pnx * pny * pnz;
    if (pnv > (int64_t)INT_MAX || !al_float_nvox_fits((uint64_t)pnv)) {
        fprintf(stderr, "qwarp: padded volume %dx%dx%d exceeds the "
        "supported voxel count\n", pnx, pny, pnz); goto pipe_done; }
    int cpnx, cpny, cpnz;
    basep = qw_zeropad(base, onx, ony, onz, xm, xp, ym, yp, zm, zp, &cpnx, &cpny, &cpnz);
    srcp  = qw_zeropad(src,  onx, ony, onz, xm, xp, ym, yp, zm, zp, &cpnx, &cpny, &cpnz);
    if (!basep || !srcp || cpnx != pnx || cpny != pny || cpnz != pnz) goto pipe_done;

    srcbp = (float *)malloc(sizeof(float) * (size_t)pnv);   /* blurred padded source */
    if (!srcbp) goto pipe_done;
    memcpy(srcbp, srcp, sizeof(float) * (size_t)pnv);
    /* the source blur IS the `-blur 0 3` step — a silent failure would change the algorithm */
    if (qw_gauss_blur_vox(srcbp, pnx, pny, pnz, QW_FWHM_TO_SIGMA(QW_SRC_BLUR_FWHM)) != 0) {
        fprintf(stderr, "qwarp: source blur failed\n"); goto pipe_done; }

    wtp = qw_weightize(basep, pnx, pny, pnz);   /* weight from the PADDED base */
    if (!wtp) goto pipe_done;

    if (qw_setup_for_improvement(&H, basep, srcbp, wtp, pnx, pny, pnz) != 0) { qw_ctx_free(&H); goto pipe_done; }
    have_H = 1;
    qw_autobbox(wtp, pnx, pny, pnz, &H.imin, &H.imax, &H.jmin, &H.jmax, &H.kmin, &H.kmax);
    if (H.imax < H.imin || H.jmax < H.jmin || H.kmax < H.kmin) goto pipe_done;  /* empty weight */

#ifdef QWARP_CHECKPOINT
    { qw_incor c0 = H.mpar; qw_incor_addto(&c0, (int)pnv, basep, H.aasrcim, wtp);
      fprintf(stderr, "[qwarp chk] padded %dx%dx%d bbox i[%d..%d] j[%d..%d] k[%d..%d] corr@identity=%.5f\n",
              pnx, pny, pnz, H.imin, H.imax, H.jmin, H.jmax, H.kmin, H.kmax, qw_incor_finalize(&c0)); }
#endif

    qw_warpomatic(&H);
    if (H.fatal) { fprintf(stderr, "qwarp: allocation failure during optimization\n"); goto pipe_done; }

#ifdef QWARP_CHECKPOINT
    { qw_incor c1 = H.mpar; qw_incor_addto(&c1, (int)pnv, basep, H.aasrcim, wtp);
      fprintf(stderr, "[qwarp chk] warpomatic done: patches=%d skipped=%d corr(blurred)=%.5f Hcost=%.5f\n",
              H.ndone, H.nskipped, qw_incor_finalize(&c1), H.Hcost); }
#endif

    /* single final reslice: warp the ORIGINAL (unblurred) padded source with WSINC5 */
    warpedp = qw_warp_apply(H.aawarp, srcp, pnx, pny, pnz, QW_WSINC5, 1.0f);
    if (!warpedp) goto pipe_done;
    QW_CHECKPOINT(stationary, warpedp, pnx, pny, pnz, "warpedpad");

    /* crop back to the stationary grid (negative pads undo the zeropad) */
    int cnx, cny, cnz;
    cropped = qw_zeropad(warpedp, pnx, pny, pnz, -xm, -xp, -ym, -yp, -zm, -zp, &cnx, &cny, &cnz);
    if (!cropped || cnx != onx || cny != ony || cnz != onz) goto pipe_done;

    {   /* build the output on the STATIONARY grid (dims + spatial metadata), float32 */
        nifti_image *out = (nifti_image *)malloc(sizeof *out);
        if (!out) { fprintf(stderr, "qwarp: failed to allocate output header\n"); goto pipe_done; }
        *out = *stationary;
        out->fname = NULL; out->iname = NULL;
        out->num_ext = 0;  out->ext_list = NULL;
        out->data = cropped; cropped = NULL;   /* transfer ownership */
        out->datatype = DT_FLOAT32; out->nbyper = 4;
        out->scl_slope = 0.0; out->scl_inter = 0.0;
        out->cal_min = 0.0;   out->cal_max = 0.0;
        *result = out;
        rc = 0;
    }

pipe_done:
    if (have_H) qw_ctx_free(&H);
    free(base); free(src); free(basep); free(srcp); free(srcbp); free(wtp); free(warpedp); free(cropped);
    return rc;
}
