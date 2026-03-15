#ifndef ALLINEATE_H
#define ALLINEATE_H

/* Affine image registration using local Pearson correlation
   Adapted from AFNI's 3dAllineate by RW Cox (public domain)
   Supports lpc+ZZ, lpa+ZZ, Hellinger, and Pearson cost functions
   with twopass coarse-to-fine optimization */

#include "nifti_io.h"
#include <string.h>

/* Cost function codes (user-facing subset) */
#define AL_COST_LPC      0  /* lpc+ZZ (cross-modal, e.g. fMRI→T1) */
#define AL_COST_LPA      1  /* lpa+ZZ (cross-modal, default for deface) */
#define AL_COST_HELLINGER 2  /* Hellinger (default for allineate, matches AFNI) */
#define AL_COST_PEARSON   3  /* Global Pearson correlation (within-modality, fast) */

/* Center-of-mass modes */
#define AL_CMASS_NONE     0  /* No center-of-mass alignment */
#define AL_CMASS_YES      1  /* Use center-of-mass for initial shift */

/* Options for nii_allineate and nii_deface */
typedef struct {
    int cost;    /* AL_COST_* code (default: AL_COST_HELLINGER) */
    int cmass;   /* AL_CMASS_* code (default: AL_CMASS_NONE) */
} al_opts;

/* Initialize options to defaults */
static inline al_opts al_opts_default(void) {
    al_opts o;
    o.cost = AL_COST_HELLINGER;
    o.cmass = AL_CMASS_NONE;
    return o;
}

/* Parse cost function name string into AL_COST_* code.
   Returns 0 on success, 1 on unrecognized name. */
static inline int al_parse_cost(const char *name, int *cost_out) {
    if (!strcmp(name, "lpc"))                        { *cost_out = AL_COST_LPC; return 0; }
    if (!strcmp(name, "lpa"))                        { *cost_out = AL_COST_LPA; return 0; }
    if (!strcmp(name, "hel"))                        { *cost_out = AL_COST_HELLINGER; return 0; }
    if (!strcmp(name, "ls") || !strcmp(name, "pearson")) { *cost_out = AL_COST_PEARSON; return 0; }
    return 1;
}

/* Parse sub-arguments (-cost, -cmass, -nocmass) from argv.
   *ac points to the last positional arg consumed; on return it points to the
   last sub-argument consumed. Returns 0 on success, 1 on error (with message). */
static inline int al_parse_subopts(int *ac, int argc, char **argv, al_opts *opts,
                                   const char *cmd_name) {
    while (*ac + 1 < argc && argv[*ac + 1][0] == '-') {
        (*ac)++;
        if (!strcmp(argv[*ac], "-cmass")) {
            opts->cmass = AL_CMASS_YES;
        } else if (!strcmp(argv[*ac], "-nocmass")) {
            opts->cmass = AL_CMASS_NONE;
        } else if (!strcmp(argv[*ac], "-cost")) {
            (*ac)++;
            if (*ac >= argc) {
                fprintf(stderr, "%s -cost requires a cost function name\n", cmd_name);
                return 1;
            }
            if (al_parse_cost(argv[*ac], &opts->cost)) {
                fprintf(stderr, "Unknown cost function '%s' (use: lpc, lpa, hel, ls)\n", argv[*ac]);
                return 1;
            }
        } else {
            /* Not a recognized sub-argument, back up */
            (*ac)--;
            break;
        }
    }
    return 0;
}

/* Register source image to base image grid using affine (12 DOF) alignment.
   source: the moving image (will be modified in-place: data replaced, dims updated)
   base: the stationary/reference image
   opts: registration options (cost function, cmass, etc.)
   Returns 0 on success, nonzero on error. */
int nii_allineate(nifti_image *source, nifti_image *base, al_opts opts);

/* Deface: register template to input, warp mask to input space, zero masked voxels.
   input: the image to deface (modified in-place, stays in its own space)
   tmpl: template image (moving image for registration)
   mask: mask in template space (non-zero = keep, zero = remove face)
   opts: registration options (cost function, cmass)
   Returns 0 on success, nonzero on error. */
int nii_deface(nifti_image *input, nifti_image *tmpl, nifti_image *mask, al_opts opts);

#endif /* ALLINEATE_H */
