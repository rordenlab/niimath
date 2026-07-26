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

// -moco [-1Dfile <path>]: ordinary chain operation.  Registers every sub-brick of a 4D image
// to sub-brick 0 and replaces `nim` with the corrected series.  `nim` must be float32 and 4D
// with nt > 1.  When `par_path` is non-NULL the six motion parameters per sub-brick are written
// there in AFNI `-1Dfile` format (roll pitch yaw dS dL dP).  The caller must supply a path ending
// in `.1D`; that is enforced at parse time so the parameter file can never name a supported NIfTI
// output and be silently replaced by the image writer.
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
int nii_moco(nifti_image *nim, const char *par_path);

#ifdef __cplusplus
}
#endif

#endif // MOCO_H
