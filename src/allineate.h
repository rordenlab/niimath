#ifndef ALLINEATE_H
#define ALLINEATE_H

/* Affine (12 DOF) image registration
   Adapted from AFNI's 3dAllineate by RW Cox (public domain)
   Supports Hellinger (default), lpc, lpa, and Pearson (ls) cost functions
   with twopass coarse-to-fine optimization and TOHD blok local correlation.
   Compile with -DAL_LPC_MICHO to use lpc+ZZ combined cost (slower, adds
   Hellinger/MI/NMI/CrA helper costs with +ZZ pure-lpc refinal pass). */

#include "nifti_io.h"
#include <string.h>

/* Cost function codes (user-facing subset) */
#define AL_COST_LPC      0  /* lpc (cross-modal, e.g. fMRI→T1) */
#define AL_COST_LPA      1  /* lpa (cross-modal, default for deface) */
#define AL_COST_HELLINGER 2  /* Hellinger (default for allineate) */
#define AL_COST_PEARSON   3  /* Global Pearson correlation (within-modality, fast) */

/* Center-of-mass modes */
#define AL_CMASS_NONE     0  /* No center-of-mass alignment */
#define AL_CMASS_YES      1  /* Use center-of-mass for initial shift */

/* Warp type codes (number of free parameters) */
#define AL_WARP_SHIFT_ONLY          3  /* shift_only / sho: 3 DOF */
#define AL_WARP_SHIFT_ROTATE        6  /* shift_rotate / shr: 6 DOF */
#define AL_WARP_SHIFT_ROTATE_SCALE  9  /* shift_rotate_scale / srs: 9 DOF */
#define AL_WARP_AFFINE_GENERAL     12  /* affine_general / aff: 12 DOF (default) */

/* Interpolation codes (shared between -interp and -final) */
#define AL_INTERP_NN       0  /* Nearest-neighbor */
#define AL_INTERP_LINEAR   1  /* Trilinear (default for -interp, matching AFNI) */
#define AL_INTERP_CUBIC    3  /* Tricubic */
#define AL_INTERP_DEFAULT -1  /* Use mode-appropriate default (cubic for allineate, linear for deface) */

/* Options for nii_allineate and nii_deface */
typedef struct {
    int cost;              /* AL_COST_* code (default: AL_COST_HELLINGER) */
    int cmass;             /* AL_CMASS_* code (default: AL_CMASS_NONE) */
    int source_automask;   /* if nonzero, fill outside of source automask with noise */
    int interp;            /* AL_INTERP_* for fine-pass matching (default: LINEAR) */
    int final_interp;      /* AL_INTERP_* for output reslicing (default: AL_INTERP_DEFAULT) */
    int warp;              /* AL_WARP_* DOF count (default: AL_WARP_AFFINE_GENERAL = 12) */
    const char *skullstrip; /* brain mask file for -skullstrip mode (NULL = normal registration) */
} al_opts;

/* Initialize options to defaults */
static inline al_opts al_opts_default(void) {
    al_opts o;
    o.cost = AL_COST_HELLINGER;
    o.cmass = AL_CMASS_NONE;
    o.source_automask = 0;
    o.interp = AL_INTERP_LINEAR;
    o.final_interp = AL_INTERP_DEFAULT;
    o.warp = AL_WARP_AFFINE_GENERAL;
    o.skullstrip = NULL;
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

/* Parse interpolation name string into AL_INTERP_* code.
   Returns 0 on success, 1 on unrecognized name. */
static inline int al_parse_interp(const char *name, int *code_out) {
    if (!strcmp(name, "NN") || !strcmp(name, "nearest") ||
        !strcmp(name, "nearestneighbour") || !strcmp(name, "nearestneighbor"))
        { *code_out = AL_INTERP_NN; return 0; }
    if (!strcmp(name, "linear") || !strcmp(name, "trilinear"))
        { *code_out = AL_INTERP_LINEAR; return 0; }
    if (!strcmp(name, "cubic") || !strcmp(name, "tricubic"))
        { *code_out = AL_INTERP_CUBIC; return 0; }
    return 1;
}

/* Parse warp type name into AL_WARP_* code.
   Returns 0 on success, 1 on unrecognized name. */
static inline int al_parse_warp(const char *name, int *warp_out) {
    if (!strcmp(name, "shift_only") || !strcmp(name, "sho"))
        { *warp_out = AL_WARP_SHIFT_ONLY; return 0; }
    if (!strcmp(name, "shift_rotate") || !strcmp(name, "shr"))
        { *warp_out = AL_WARP_SHIFT_ROTATE; return 0; }
    if (!strcmp(name, "shift_rotate_scale") || !strcmp(name, "srs"))
        { *warp_out = AL_WARP_SHIFT_ROTATE_SCALE; return 0; }
    if (!strcmp(name, "affine_general") || !strcmp(name, "aff"))
        { *warp_out = AL_WARP_AFFINE_GENERAL; return 0; }
    return 1;
}

/* Return human-readable name for an AL_WARP_* code. */
static inline const char *al_warp_name(int warp) {
    switch (warp) {
        case AL_WARP_SHIFT_ONLY:         return "shift_only";
        case AL_WARP_SHIFT_ROTATE:       return "shift_rotate";
        case AL_WARP_SHIFT_ROTATE_SCALE: return "shift_rotate_scale";
        default:                         return "affine_general";
    }
}

/* Parse sub-arguments from argv.
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
        } else if (!strcmp(argv[*ac], "-source_automask")) {
            opts->source_automask = 1;
        } else if (!strcmp(argv[*ac], "-nearest") || !strcmp(argv[*ac], "-NN")) {
            opts->final_interp = AL_INTERP_NN;
        } else if (!strcmp(argv[*ac], "-linear") || !strcmp(argv[*ac], "-trilinear")) {
            opts->final_interp = AL_INTERP_LINEAR;
        } else if (!strcmp(argv[*ac], "-cubic") || !strcmp(argv[*ac], "-tricubic")) {
            opts->final_interp = AL_INTERP_CUBIC;
        } else if (!strcmp(argv[*ac], "-skullstrip")) {
            (*ac)++;
            if (*ac >= argc) {
                fprintf(stderr, "%s -skullstrip requires a brain mask filename\n", cmd_name);
                return 1;
            }
            opts->skullstrip = argv[*ac];
        } else if (!strcmp(argv[*ac], "-warp")) {
            (*ac)++;
            if (*ac >= argc) {
                fprintf(stderr, "%s -warp requires a type (sho, shr, srs, aff)\n", cmd_name);
                return 1;
            }
            if (al_parse_warp(argv[*ac], &opts->warp)) {
                fprintf(stderr, "Unknown warp '%s' (use: sho, shr, srs, aff)\n", argv[*ac]);
                return 1;
            }
        } else if (!strcmp(argv[*ac], "-interp")) {
            (*ac)++;
            if (*ac >= argc) {
                fprintf(stderr, "%s -interp requires an interpolation name\n", cmd_name);
                return 1;
            }
            if (al_parse_interp(argv[*ac], &opts->interp)) {
                fprintf(stderr, "Unknown interp '%s' (use: NN, linear, cubic)\n", argv[*ac]);
                return 1;
            }
        } else if (!strcmp(argv[*ac], "-final")) {
            (*ac)++;
            if (*ac >= argc) {
                fprintf(stderr, "%s -final requires an interpolation name\n", cmd_name);
                return 1;
            }
            if (al_parse_interp(argv[*ac], &opts->final_interp)) {
                fprintf(stderr, "Unknown final interp '%s' (use: NN, linear, cubic)\n", argv[*ac]);
                return 1;
            }
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
   opts: registration options (cost function, cmass, interpolation, etc.)
   final_interp default: cubic.
   Returns 0 on success, nonzero on error. */
int nii_allineate(nifti_image *source, nifti_image *base, al_opts opts);

/* Deface/skullstrip: register template to input, warp mask to input space,
   set voxels where warped mask < 0.5 to the input's minimum value.
   input: the image to modify (modified in-place, stays in its own space)
   tmpl: template image (moving image for registration)
   mask: mask in template space (non-zero = keep)
   opts: registration options (cost function, cmass)
   final_interp default: linear (to avoid ringing in the mask).
   Returns 0 on success, nonzero on error. */
int nii_deface(nifti_image *input, nifti_image *tmpl, nifti_image *mask, al_opts opts);

#endif /* ALLINEATE_H */
