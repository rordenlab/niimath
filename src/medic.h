// medic.h - MEDIC (Multi-Echo DIstortion Correction) for niimath
//
// Clean-room emulation of the MEDIC workflow described in Van et al., Imaging Neuroscience 4
// (2026), doi:10.1162/IMAG.a.1262. No Warpkit implementation, test, build product or debug symbol
// was read; every convention not fixed by the paper was measured through the public executables
// and is recorded in test/medic_reference_manifest.md.
//
// Phase unwrapping is the MIT ROMEO port in romeo.c, reached through its in-memory frame API
// (romeo_unwrap_frame); MEDIC does not duplicate ROMEO and does not shell out per frame.
//
// Guarded by HAVE_MEDIC, which requires HAVE_ROMEO.

#ifndef MEDIC_H
#define MEDIC_H

#ifdef __cplusplus
extern "C" {
#endif

#include "core.h"

// Entry point for "niimath --medic ...", a terminal subcommand like --dtifit and --qc: it parses
// its own argv and returns EXIT_SUCCESS/EXIT_FAILURE.  Writes <prefix>_fieldmaps_native,
// <prefix>_fieldmaps and <prefix>_displacementmaps.
int nii_medic(int argc, char *argv[]);

// -unwarp <displacement-map> <axis>: ordinary chain operation.  Resamples the working image
// through a scalar displacement map in millimetres along one phase-encoding axis.  `nim` must be
// float32; it is replaced in place.  Returns 0 on success.
int medic_unwarp(nifti_image *nim, const char *mapfile, const char *axis);

#ifdef __cplusplus
}
#endif

#endif // MEDIC_H
