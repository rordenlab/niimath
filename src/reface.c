/*----------------------------------------------------------------------------
 * reface.c — face-replacement (anonymization) compositing stage.
 *
 * CLEAN-ROOM reimplementation of the AFNI afni_refacer2.csh "reface" recipe.
 * See reface.h for the provenance/licensing boundary. No AFNI source is copied;
 * the Gaussian blur is niimath's own BSD/PD nifti_smooth_gauss_f32 used as a
 * normalized convolution, NOT a port of AFNI's 3dBlurInMask.
 *
 * Serial and deterministic: `ifac` uses a fixed-order scalar accumulation and
 * the separable Gaussian parallelizes only across independent rows, so the whole
 * operation is bit-identical at -p 1 vs -p N.
 *--------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include "reface.h"
#include "core32.h"   /* nifti_smooth_gauss_f32 */
#include "al_size_guard.h"

/* FWHM -> Gaussian sigma: FWHM = 2*sqrt(2*ln2)*sigma. */
#define REFACE_FWHM_MM        2.666f
#define REFACE_FWHM_TO_SIGMA  2.3548200450309493f  /* 2*sqrt(2*ln2) */

/* Finite guard usable under -ffast-math -fno-finite-math-only (isfinite() is
   unreliable there): true only for a real magnitude in [-FLT_MAX, FLT_MAX].
   fabsf(NaN) > FLT_MAX is false, so NaN and +/-Inf both fail this test. */
static int reface_finitef(float v) { return fabsf(v) <= FLT_MAX; }
static int reface_finited(double v) { return fabs(v) <= DBL_MAX; }

/* Overflow-checked voxel count for a public API boundary: a C-API caller can pass an
   image whose dims multiply past SIZE_MAX (wrapping an allocation and driving loops out
   of bounds). Writes nx*ny*nz to *out and returns 0 only if positive and non-wrapping
   (including the later *sizeof(float) malloc). Portable (no __builtin), matching mc_mul.
   Also caps at INT_MAX voxels: the Gaussian backend (nifti_smooth_gauss_f32) takes int
   counts and rejects a larger volume, so fail immediately and cheaply here rather than
   scanning the input and attempting three huge float allocations first (mirrors the fast
   engine's / qwarp's INT_MAX voxel policy). Real brain grids are nowhere near this. */
static int reface_nvox(int nx, int ny, int nz, size_t *out) {
    if (nx < 1 || ny < 1 || nz < 1) return 1;
    size_t a = (size_t)nx, b = (size_t)ny, c = (size_t)nz;
    if (b > SIZE_MAX / a) return 1;
    size_t ab = a * b;
    if (c > SIZE_MAX / ab) return 1;
    size_t n = ab * c;
    if (n > (size_t)INT_MAX || !al_float_nvox_fits((uint64_t)n))
        return 1;   /* int backend + float-buffer byte limit; fail cheap, pre-alloc */
    *out = n;
    return 0;
}

void reface_isola(float *v, int nx, int ny, int nz) {
    size_t nvox;
    if (!v || reface_nvox(nx, ny, nz, &nvox)) return;
    /* Snapshot the original nonzero pattern so a removal cannot cascade into a
       neighbor's isolation test (matches AFNI's single-pass -isola). One byte
       per voxel; on OOM skip the cleanup (denoiser only, never required). */
    unsigned char *was_nz = (unsigned char *)malloc(nvox);
    if (!was_nz) return;
    for (size_t i = 0; i < nvox; i++) was_nz[i] = (v[i] != 0.0f) ? 1u : 0u;
    size_t nxy = (size_t)nx * ny;
    for (int z = 0; z < nz; z++) {
        for (int y = 0; y < ny; y++) {
            for (int x = 0; x < nx; x++) {
                size_t idx = (size_t)z * nxy + (size_t)y * nx + x;
                if (!was_nz[idx]) continue;
                int isolated = 1;
                for (int dz = -1; dz <= 1 && isolated; dz++) {
                    int zz = z + dz; if (zz < 0 || zz >= nz) continue;
                    for (int dy = -1; dy <= 1 && isolated; dy++) {
                        int yy = y + dy; if (yy < 0 || yy >= ny) continue;
                        for (int dx = -1; dx <= 1; dx++) {
                            if (dx == 0 && dy == 0 && dz == 0) continue;
                            int xx = x + dx; if (xx < 0 || xx >= nx) continue;
                            if (was_nz[(size_t)zz * nxy + (size_t)yy * nx + xx]) { isolated = 0; break; }
                        }
                    }
                }
                if (isolated) v[idx] = 0.0f;
            }
        }
    }
    free(was_nz);
}

int reface_apply(const nifti_image *subject, const nifti_image *shell,
                 nifti_image **result) {
    if (result) *result = NULL;
    if (!result || !subject || !shell || !subject->data || !shell->data) {
        fprintf(stderr, "reface: NULL image\n"); return 1;
    }
    if (subject->datatype != DT_FLOAT32 || shell->datatype != DT_FLOAT32) {
        fprintf(stderr, "reface: inputs must be float32\n"); return 1;
    }
    int nx = subject->nx, ny = subject->ny, nz = subject->nz;
    size_t nvox;
    if (reface_nvox(nx, ny, nz, &nvox)) {   /* overflow-checked (public API boundary) */
        fprintf(stderr, "reface: bad or too-large subject dims\n"); return 1;
    }
    if (subject->nvox != nvox) { fprintf(stderr, "reface: subject must be single-volume 3D\n"); return 1; }
    if (shell->nx != nx || shell->ny != ny || shell->nz != nz ||
        shell->nvox != nvox) {
        fprintf(stderr, "reface: shell/subject grids differ\n"); return 1;
    }
    const float *sub = (const float *)subject->data;
    const float *sh  = (const float *)shell->data;

    /* --- ifac: brightness match (AFNI 3dBrickStat -non-zero -mean) ---------
       ibar = mean subject over {shell>0 && subject!=0}; mbar = mean shell over
       {shell>0}. Fixed-order double accumulation -> deterministic. A non-finite shell
       voxel (a malformed C-API input) is rejected rather than silently taking the
       preserve branch (NaN>0 and NaN<0 are both false). */
    double isum = 0.0, msum = 0.0;
    size_t icnt = 0, mcnt = 0;
    for (size_t i = 0; i < nvox; i++) {
        float c = sh[i];
        if (!reface_finitef(c)) { fprintf(stderr, "reface: non-finite shell voxel\n"); return 1; }
        if (c > 0.0f) {
            msum += (double)c; mcnt++;
            float a = sub[i];
            if (a != 0.0f) { isum += (double)a; icnt++; }
        }
    }
    if (mcnt == 0 || icnt == 0) {
        fprintf(stderr, "reface: shell has no positive region overlapping subject data (cannot scale)\n");
        return 1;
    }
    double ibar = isum / (double)icnt;
    double mbar = msum / (double)mcnt;
    /* Magnitude guards, NOT isfinite(): under -ffast-math -fno-finite-math-only the
       compiler may fold isfinite() to true (see reface_finitef). A NaN subject voxel
       (NaN != 0.0f passes the accumulation filter) would then poison ibar and slip past
       an isfinite() check; the magnitude guard fails closed on it. */
    if (!reface_finited(ibar) || !reface_finited(mbar) || mbar <= 0.0) {
        fprintf(stderr, "reface: degenerate brightness match (ibar=%g mbar=%g)\n", ibar, mbar);
        return 1;
    }
    double ifac = ibar / mbar;
    if (!reface_finited(ifac)) { fprintf(stderr, "reface: non-finite ifac\n"); return 1; }

    /* --- composite + masked-blur buffers ---------------------------------- */
    float *out = (float *)malloc(nvox * sizeof(float));  /* composite, becomes final */
    float *num = (float *)malloc(nvox * sizeof(float));  /* blur(out * mask)          */
    float *den = (float *)malloc(nvox * sizeof(float));  /* blur(mask)                */
    if (!out || !num || !den) {
        free(out); free(num); free(den);
        fprintf(stderr, "reface: out of memory\n"); return 1;
    }
    for (size_t i = 0; i < nvox; i++) {
        float c = sh[i], o;
        if (c > 0.0f)      o = (float)((double)c * ifac);   /* insert scaled shell   */
        else if (c < 0.0f) o = 0.0f;                        /* outer shell -> zero   */
        else { float a = sub[i]; o = (a > 0.0f) ? a : 0.0f; } /* preserve subject    */
        out[i] = o;
        int inside = (c > 0.0f);
        num[i] = inside ? o : 0.0f;
        den[i] = inside ? 1.0f : 0.0f;
    }

    /* Normalized Gaussian convolution inside the shell>0 mask (AFNI 3dBlurInMask
       -preserve equivalent): final = blur(out*mask)/blur(mask) where mask, else
       out. Uses the subject pixdims; a nonpositive pixdim falls back to 1 mm. */
    float dx = (subject->dx > 0.0f) ? subject->dx : 1.0f;
    float dy = (subject->dy > 0.0f) ? subject->dy : 1.0f;
    float dz = (subject->dz > 0.0f) ? subject->dz : 1.0f;
    float sig = REFACE_FWHM_MM / REFACE_FWHM_TO_SIGMA;
    if (nifti_smooth_gauss_f32(num, nx, ny, nz, 1, dx, dy, dz, sig, sig, sig, -6.0f) ||
        nifti_smooth_gauss_f32(den, nx, ny, nz, 1, dx, dy, dz, sig, sig, sig, -6.0f)) {
        free(out); free(num); free(den);
        fprintf(stderr, "reface: masked blur failed\n"); return 1;
    }
    for (size_t i = 0; i < nvox; i++) {
        if (sh[i] > 0.0f) {
            float d = den[i];
            if (d > 1e-6f) {
                float val = num[i] / d;
                if (reface_finitef(val)) out[i] = val;
            }
        }
    }
    free(num); free(den);

    /* Final finiteness guard (fail closed rather than emit a corrupt scan). */
    for (size_t i = 0; i < nvox; i++) {
        if (!reface_finitef(out[i])) {
            free(out);
            fprintf(stderr, "reface: non-finite output voxel\n"); return 1;
        }
    }

    /* Wrap `out` in a float32 image on the subject grid. Copy the subject header
       and swap in our data; NULL out the owned name/extension pointers so
       nifti_image_free() on the result frees only its own buffer. subject is
       already plain float32 single-volume, so datatype/nbyper/dims are correct. */
    nifti_image *r = (nifti_image *)malloc(sizeof *r);
    if (!r) { free(out); fprintf(stderr, "reface: out of memory\n"); return 1; }
    *r = *subject;
    r->fname = NULL; r->iname = NULL;
    r->num_ext = 0;  r->ext_list = NULL;
    r->data = out;
    r->scl_slope = 0.0f; r->scl_inter = 0.0f;
    r->cal_min = 0.0f;   r->cal_max = 0.0f;
    *result = r;
    return 0;
}
