// resample.h - one-shot resampling of a 4D BOLD series (-resample), for fastprep.
//
// Every output voxel is pulled through, in order: an optional greedy nonlinear chain
// (warp then affine, greedy's `-r warp.nii.gz aff.mat` order), an optional rigid/affine
// world-mm pull matrix (allineate -savemat "fixed_to_moving"), the input's own voxel grid, and a
// per-volume -moco correction row.  The input is then sampled ONCE with a cubic B-spline
// (prefiltered, scipy map_coordinates order=3 equivalent; 0 outside the field of view).
// Built with -moco (HAVE_MOCO); it reuses its per-volume geometry.

#ifndef RESAMPLE_H
#define RESAMPLE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "core.h"

// Susceptibility voxel-shift map on the INPUT grid (float32, voxels; sign already folded in as
// fMRIPrep's vsm = Hz * readout * (-1 for a '-' PE direction)), added along `axis` after motion.
// A 4D map (one volume per frame, e.g. MEDIC) is sampled in each frame's own native position.
typedef struct {
	const nifti_image *map;
	int axis;       // 0, 1, 2 = i, j, k
	int jacobian;   // multiply by 1 + d(vsm)/d(axis)
} resample_vsm;

typedef struct {
	const nifti_image *grid;   // output geometry (header only)
	double res;                // > 0: isotropic voxel size replacing the grid's, same field of view
	const double *hmc;         // nt rows of six -moco correction parameters, or NULL
	const double *ftm;         // 4x4 row-major world pull matrix (output/T1w -> input world), or NULL
	const nifti_image *warp;   // greedy displacement field (LPS mm, 5D, 3 components), or NULL
	const double *aff;         // 4x4 row-major greedy RAS affine (fixed -> moving); with warp
	const resample_vsm *vsm;   // or NULL
	const nifti_image *bwarp;  // greedy field in BOLD world (undistorted -> distorted, LPS mm), after ftm; or NULL
} resample_opts;

// Writes the resampled float32 series to `out` (format per FSLOUTPUTTYPE/gz); `nim` (float32,
// 3D or 4D) is left untouched.  Returns 0 on success.
int nii_resample(nifti_image *nim, const resample_opts *o, const char *out, gzModes gz);

// In-memory helpers (caller frees).  resample_hmc: the series motion-corrected onto its own grid
// (nt x nvox; hmc, vsm and bwarp may be NULL).  resample_volume: 3D `src` pulled onto `grid` through a world pull matrix (grid
// world -> src world, 4x4 row-major; NULL = identity).  nii_read_float32: read + convert to float32.
float *resample_hmc(const nifti_image *nim, const double *hmc, const resample_vsm *vsm, const nifti_image *bwarp);
float *resample_volume(const nifti_image *src, const nifti_image *grid, const double *pull);
nifti_image *nii_read_float32(const char *fn);

#ifdef __cplusplus
}
#endif

#endif // RESAMPLE_H
