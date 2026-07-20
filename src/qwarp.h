/*----------------------------------------------------------------------------
 * qwarp.h — minimal nonlinear (3dQwarp) registration stage: public API.
 *
 * ATTRIBUTED PUBLIC-DOMAIN PORT. This module is a direct, faithful port of the
 * nonlinear-warp code from AFNI (public domain, NIMH/NIH), by Robert W. Cox
 * ("Zhark"):
 *     3dQwarp.c, mri_nwarp.c, thd_incorrelate.c, thd_cliplevel.c,
 *     thd_automask.c, edt_volpad.c, mri_genalign_util.c
 * pinned at AFNI revision 506e48403 (AFNI_26.1.00-6-g506e48403).
 * The derivative-free optimizer is M.J.D. Powell's NEWUOA (powell_newuoa.c),
 * the same public-domain routine AFNI itself uses.
 *
 * This is NOT clean-room code (unlike the fast engine in coreg_fast.c, which is
 * deliberately designed only from published papers to avoid AFNI's GPL parts).
 * Direct porting is limited to those public-domain files. AFNI's GPL/Medical
 * College of Wisconsin exceptions and niimath/src/GPL/ are explicitly excluded;
 * see AGENTS.md "Clean-room provenance" for the precise boundary.
 *
 * Scope (initial): the single AFNI operation
 *     3dQwarp -blur 0 3 -source <moving> -base <stationary> -prefix <output>
 * i.e. the nonlinear stage only. The moving image must already be unifized,
 * skull-stripped, affine-aligned, and sampled on the stationary grid.
 *--------------------------------------------------------------------------*/
#ifndef QWARP_H
#define QWARP_H

#include "nifti_io.h"

/* Nonlinearly warp `moving` onto `stationary` (equivalent to the AFNI call above).
 *
 * Contract:
 *   - `moving` and `stationary` must share a grid: equal 3D dims, finite positive
 *     voxel sizes, a usable (sform- or qform-coded) and mutually-agreeing
 *     voxel->world transform, and single-volume (3D) data. Neither is mutated.
 *   - On success returns 0 and sets *result to a NEWLY ALLOCATED nifti_image on
 *     the stationary grid (stationary's dims + spatial metadata), datatype
 *     float32, holding the warped source. The CALLER owns *result and must free
 *     it with nifti_image_free().
 *   - On any failure returns nonzero and sets *result to NULL — no partial state,
 *     no allocation leaked (failure is atomic).
 *   - Serial and non-reentrant: the NEWUOA callback uses process-global context.
 *     Do not overlap qwarp_run() calls from multiple host threads.
 *
 * The dense index-space displacement field used internally is private to qwarp.c
 * (not exposed) — future work (warp output, transform composition) will surface
 * it deliberately.
 */
int qwarp_run(const nifti_image *moving, const nifti_image *stationary,
              nifti_image **result);

#endif /* QWARP_H */
