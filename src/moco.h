// moco.h - rigid-body motion correction for 4D datasets (-moco)
//
// Clean-room implementation of the algorithm published in
//   Cox RW & Jesmanowicz A, "Real-Time 3D Image Registration for Functional MRI",
//   Magn Reson Med 42:1014-1018 (1999).
// The paper supplies the objective function, the repeated-linearization optimizer, the
// four-shear factorization and its Appendix proof, and the translation folding.  AFNI's
// 3dvolreg implementation is GPL-2 and was NOT read, translated or paraphrased; it was used
// only as a black-box oracle.  Every convention the paper does not fix was measured through
// the public 3dvolreg/3drotate executables and is recorded in the moco_bench repository's test/moco_reference_manifest.md.
//
// Guarded by HAVE_MOCO.

#ifndef MOCO_H
#define MOCO_H

#ifdef __cplusplus
extern "C" {
#endif

#include "core.h"

// -moco [-ref <n|image>] [-1Dfile <path>]: ordinary chain operation.  Registers every sub-brick
// of a 4D image to a reference and replaces `nim` with the corrected series.  `nim` must be
// float32 and 4D with nt > 1.  When `par_path` is non-NULL the six motion parameters per
// sub-brick are written there in AFNI `-1Dfile` format (roll pitch yaw dS dL dP).  The caller
// must supply a path ending in `.1D`; that is enforced at parse time so the parameter file can
// never name a supported NIfTI output and be silently replaced by the image writer.
//
// The reference is chosen by exactly one of:
//   * `ref_vol` >= 0 with `ref_file` NULL -- sub-brick `ref_vol` of the input series (0 is the
//     historical default).  That sub-brick is copied through unchanged and its parameter row is
//     all zeros, as it always has been for volume 0.
//   * `ref_file` non-NULL -- an EXTERNAL reference image (`ref_vol` is ignored).  It must be on
//     the same voxel grid as the input: identical nx/ny/nz and a voxel-to-world transform
//     agreeing within 0.001 mm, the same gate --qc and --medic use.  -moco registers inside the
//     input's own voxel grid, so a reference on any other grid cannot be honoured without a
//     resampling step whose conventions are outside the measured contract; it is rejected with a
//     diagnostic rather than silently resliced.  A 4D reference contributes its volume 0.  With
//     an external reference NO sub-brick is copied through: every one of the nt volumes is
//     registered and gets a non-zero parameter row.
//
// Failure behaviour, stated precisely because a chain operation cannot make the LATER image
// write atomic by itself:
//   * on any failure inside this call `nim` is left untouched, no parameter file is created or
//     modified (it is written through an exclusive sibling temporary and renamed only on
//     success), and a diagnostic is printed;
//   * the parameter file is published BEFORE the corrected image, which niimath's output stage
//     writes after this call returns.  A later image-write failure therefore leaves a valid
//     parameter file with no image.  Aliasing `par_path` with the input or the output image is
//     rejected at parse time.
// Returns 0 on success, non-zero on failure.
int nii_moco(nifti_image *nim, const char *par_path, int ref_vol, const char *ref_file);

// -moco -relative: MEASUREMENT ONLY.  Fits every sub-brick t (1..nt-1) onto sub-brick t-1 and
// writes those parameters; `nim` is NOT modified and NO image is written -- on the command line
// `-relative` is a flag and the trailing positional that would otherwise name the output image is
// `rel_path` (the caller enforces the `.1D` suffix and suppresses the image write).  The base of each
// pair is the ORIGINAL predecessor, never a corrected one, which keeps the pairs independent (so
// the loop is parallel and thread-count invariant) and makes the numbers raw frame-to-frame motion
// rather than residual drift after correction.  Row 0 is all zeros: volume 0 has no predecessor.
//
// Two files are published, both through exclusive sibling temporaries renamed on success:
//   * `path` -- text, nt rows of six %12.8f fields (roll pitch yaw dS dL dP), the same convention
//     and column order as -1Dfile but at the precision the estimator actually carries;
//   * `path` + ".bin" -- the same numbers as nt*6 raw little-endian float64, row-major, no header,
//     so a reader is np.fromfile(path + ".bin").reshape(-1, 6).
// Neither is published if writing fails; only a failure BETWEEN the two renames can leave the text
// file without its companion.  Because the estimator is shared with nii_moco (moco_base_setup +
// moco_fit_one), the two modes cannot drift apart.
//
// Rebuilding the weight, the six derivative images and the normal equations for every pair makes
// this roughly an order of magnitude more work per volume than the ordinary mode.
// Returns 0 on success, non-zero on failure.
int nii_moco_relative(nifti_image *nim, const char *rel_path);

#ifdef __cplusplus
}
#endif

#endif // MOCO_H
