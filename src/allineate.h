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
#define AL_COST_LPA      1  /* lpa (cross-modal) */
#define AL_COST_HELLINGER 2  /* Hellinger — the ORDINARY engine's cost (-cost hel). The CLI
                                default cost is the fast engine (see AL_ENGINE_* below); these
                                AL_COST_* codes select the ordinary AFNI-style engine. */
#define AL_COST_PEARSON   3  /* Global Pearson correlation (within-modality) */
#define AL_COST_NMI       4  /* Normalized mutual information, H(x,y)/(H(x)+H(y)) — an
                                entropy-ratio alternative to Hellinger for images with a
                                huge shared background/zero bin (masked or skull-stripped
                                data), where Hellinger's sqrt(joint) term is dominated by
                                that one concentrated cell. Minimized (1 == independent). */

/* Center-of-mass modes */
#define AL_CMASS_NONE     0  /* No center-of-mass alignment */
#define AL_CMASS_YES      1  /* Use center-of-mass for initial shift */

/* Registration engine selector (al_opts.fast). The fast engine (coreg_fast.c) is a distinct
   estimator that is the CLI DEFAULT (via `-cost fast`); the value also carries which fast cost
   to use. The CLI default lives in the dispatch, not here — `al_opts_default()` leaves
   `.fast = 0` and the command dispatch sets AL_ENGINE_FAST_X when no -cost is given. */
#define AL_ENGINE_ALLINEATE 0  /* ordinary AFNI-style allineate engine (nii_allineate; -cost hel/lpc/lpa/ls) */
#define AL_ENGINE_FAST_CR   1  /* fast engine, correlation-ratio cost  (-cost fastcr) */
#define AL_ENGINE_FAST_HEL  2  /* fast engine, Hellinger cost          (-cost fasthel) */
#define AL_ENGINE_FAST_X    3  /* fast engine, HEL/CR coarse compete  (-cost fast/-cost fastx; default) */

/* al_opts.cli_set bits — which override options the user explicitly passed. Split
   matching-interp (-interp) from output-interp (-final/-nearest/-linear/-cubic): the
   fast engine honors -final for its one output reslice but ignores matching -interp. */
#define AL_CLI_COST   0x1u
#define AL_CLI_WARP   0x2u
#define AL_CLI_INTERP 0x4u   /* -interp (fine-pass MATCHING interpolation) */
#define AL_CLI_FINAL  0x8u   /* -final / -nearest / -linear / -cubic (OUTPUT interpolation) */
#define AL_CLI_CMASS  0x10u  /* -cmass / -nocmass */

/* al_parse_subopts capability bits — the set of option GROUPS a command actually implements.
   Each command passes its capabilities so an option it does not honor is rejected AT PARSE
   TIME (before images load) rather than silently ignored (the `-deface` silent-no-op bug
   class). Intra-command mode restrictions (e.g. the fast engine within -allineate rejecting
   -warp/-zoom — it DOES honor the -com/-sym header seeds) are separate and enforced at dispatch. */
#define AL_CAP_TUNING  0x1u   /* -cost <normal>/-warp/-interp/-cmass/-nocmass/-source_automask/-dark_automask */
#define AL_CAP_FINAL   0x2u   /* -final / -nearest / -linear / -cubic (output interpolation) */
#define AL_CAP_FAST    0x4u   /* -cost fast / fastx / fasthel / fastcr (fast engine) */
#define AL_CAP_MASTER  0x8u   /* -master <grid> */
#define AL_CAP_MATRIX  0x10u  /* -savemat / -applymat */
#define AL_CAP_SEED    0x20u  /* -com / -sym / -symd / -symb / -nosagseed / -zoom */
#define AL_CAP_FILL    0x40u  /* -fill zero/nan/auto (out-of-FOV output fill; not for -deface) */
#define AL_CAP_WEIGHT  0x80u  /* -weight <img> (AFNI 3dAllineate graded base-space weight; not for -deface) */
#define AL_CAP_ALL     (AL_CAP_TUNING|AL_CAP_FINAL|AL_CAP_FAST|AL_CAP_MASTER|AL_CAP_MATRIX|AL_CAP_SEED|AL_CAP_FILL|AL_CAP_WEIGHT)

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

/* Out-of-FOV fill for the -allineate OUTPUT reslice (al_opts.fillmode). Resolved to a
   concrete float by al_resolve_fillv(). AUTO is the default: it keeps the historical 0
   fill for positive-only images (e.g. MRI, background 0 — output stays byte-identical),
   but for images whose darkest voxel is negative (CT/X-ray Hounsfield units, where air ≈
   -1000) it fills with that darkest value so out-of-FOV reads as air, not soft tissue. */
#define AL_FILL_AUTO  0  /* 0 unless the source minimum is < 0, then the source minimum (default) */
#define AL_FILL_ZERO  1  /* always 0 */
#define AL_FILL_NAN   2  /* NaN */

/* Options for nii_allineate and nii_deface */
typedef struct {
    int cost;              /* AL_COST_* code (default: AL_COST_HELLINGER) */
    int cmass;             /* AL_CMASS_* code (default: AL_CMASS_NONE) */
    int source_automask;   /* if nonzero, fill outside of source automask with noise */
    int dark_automask;     /* -dark_automask: drop matched pairs where the base or warped
                              source value is at that image's darkest value (background/pad) */
    int interp;            /* AL_INTERP_* for fine-pass matching (default: LINEAR) */
    int final_interp;      /* AL_INTERP_* for output reslicing (default: AL_INTERP_DEFAULT) */
    int fillmode;          /* AL_FILL_* out-of-FOV output fill (default: AL_FILL_AUTO) */
    int warp;              /* AL_WARP_* DOF count (default: AL_WARP_AFFINE_GENERAL = 12) */
    /* --- CLI-workflow fields below: interpreted by the host's pre-/post-processing chain
       (niimath's nifti_allineate_wrap in coreFLT.c, and the standalone allineate main.c),
       NOT by nii_allineate()/nii_deface(); a direct library caller setting those gets a
       silent no-op. EXCEPTION: `zoom` (bottom of the struct) IS read by nii_allineate()/
       al_register() as the scale-relaxation flag. These are all wired to CLI flags in both
       hosts (savemat/applymat/com/sym[/-symd/-symb via sym_deoblique]/sagseed[-nosagseed]/
       zoom as -allineate sub-options; robustfov is a standalone niimath op; skullstrip = the
       standalone's flag, provided in niimath by -deface with a brain mask). Splitting the
       estimator options from CLI state is a possible future refactor. --- */
    const char *skullstrip; /* CLI-only: -skullstrip brain mask (NULL = normal registration) */
    const char *weight;     /* CLI-only: -weight base(stationary)-space GRADED weight, AFNI
                               3dAllineate style, for BOTH engines. Must share the base grid (dims +
                               world frame). Normalized to [0,1] (divide by max) and used per base
                               voxel — replaces the manufactured autoweight (ordinary engine) and is
                               NOT binarized (graded, even for box costs). A voxel weighted 0 is
                               excluded; keep the out-of-ROI head attenuated (nonzero) to anchor
                               global scale (a fully-zeroed exterior risks a scale collapse, as in
                               AFNI). NULL = autoweight (ordinary) / unweighted (fast). */
    double robustfov;       /* CLI-only: -robustfov crop mm applied to the moving image (0 = off) */
    const char *savemat;    /* CLI-only: -savemat path; main.c saves the fitted affine as JSON
                               (via nii_last_affine()). NULL = don't save. */
    const char *applymat;   /* CLI-only: -applymat path; apply a saved affine JSON to reslice
                               the moving image onto the stationary grid (no registration). */
    const char *master;     /* CLI-only: -master output-grid image. Register at the stationary
                               resolution, but reslice the result onto THIS grid instead (must
                               share the stationary world frame, e.g. a higher-res template).
                               NULL = output on the stationary grid. */
    int com;                /* CLI-only: -com set origin to brightness center of mass (1 = on) */
    int sym;                /* CLI-only: -sym/-symd midsagittal alignment (1 = enabled, 0 = disabled) */
    int sym_deoblique;      /* CLI-only: 0 = -sym; 1 = -symd (snap the frame axis-aligned
                               before the mirror fit); 2 = -symb (auto-compete both) */
    int sagseed;            /* CLI-only: -sagseed in-MSP rigid seed after -sym (default 1;
                               -nosagseed disables). Only acts when -sym runs with a template. */
    unsigned cli_set;       /* CLI-only: bitmask of explicitly-passed override options
                               (AL_CLI_*), so a distinct engine can reject options it
                               cannot honor (e.g. the fast engine rejects -warp/-interp; it
                               honors -cmass/-nocmass via AL_CLI_CMASS). */
    int fast;               /* CLI-only: registration engine selector (AL_ENGINE_*). Interpreted
                               by main.c (dispatches to coreg_fast_estimate); nii_allineate()
                               ignores it. */
    int zoom;               /* CLI-only: -zoom — flags abnormal size (e.g. infant vs adult
                               template): adds a global isotropic scale DOF to the -sagseed
                               seed fit AND widens the main affine's scale range to [0.5,2.0]
                               (from 0.711..1.406). Only the scale limits relax; all other
                               regularization and default (adult) behavior are unchanged. */
} al_opts;

/* Initialize options to defaults */
static inline al_opts al_opts_default(void) {
    al_opts o;
    o.cost = AL_COST_HELLINGER;
    o.cmass = AL_CMASS_NONE;
    o.source_automask = 0;
    o.dark_automask = 0;
    o.interp = AL_INTERP_LINEAR;
    o.final_interp = AL_INTERP_DEFAULT;
    o.fillmode = AL_FILL_AUTO;
    o.warp = AL_WARP_AFFINE_GENERAL;
    o.skullstrip = NULL;
    o.weight = NULL;
    o.robustfov = 0.0;
    o.savemat = NULL;
    o.applymat = NULL;
    o.master = NULL;
    o.com = 0;
    o.sym = 0;
    o.sym_deoblique = 0;
    o.sagseed = 1;
    o.cli_set = 0;
    o.fast = 0;
    o.zoom = 0;
    return o;
}

/* Parse cost function name string into AL_COST_* code.
   Returns 0 on success, 1 on unrecognized name. */
static inline int al_parse_cost(const char *name, int *cost_out) {
    if (!strcmp(name, "lpc"))                        { *cost_out = AL_COST_LPC; return 0; }
    if (!strcmp(name, "lpa"))                        { *cost_out = AL_COST_LPA; return 0; }
    if (!strcmp(name, "hel"))                        { *cost_out = AL_COST_HELLINGER; return 0; }
    if (!strcmp(name, "nmi"))                        { *cost_out = AL_COST_NMI; return 0; }
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

/* Parse out-of-FOV fill mode name into AL_FILL_* code.
   Returns 0 on success, 1 on unrecognized name. */
static inline int al_parse_fill(const char *name, int *fill_out) {
    if (!strcmp(name, "auto"))                        { *fill_out = AL_FILL_AUTO; return 0; }
    if (!strcmp(name, "zero") || !strcmp(name, "0"))  { *fill_out = AL_FILL_ZERO; return 0; }
    if (!strcmp(name, "nan") || !strcmp(name, "NaN")) { *fill_out = AL_FILL_NAN;  return 0; }
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

/* Canonical (re-parseable, matches al_parse_cost) name for an AL_COST_* code. */
static inline const char *al_cost_name(int cost) {
    switch (cost) {
        case AL_COST_LPC:      return "lpc";
        case AL_COST_LPA:      return "lpa";
        case AL_COST_PEARSON:  return "ls";
        case AL_COST_NMI:      return "nmi";
        default:               return "hel";
    }
}

/* Parse trailing `-sub` options for a `-allineate`/`-deface` command (e.g.
   `-allineate base.nii -cost fast -master hires.nii`). Consumes options while the
   NEXT token starts with '-' and is recognized; stops (backing up) at the first
   unrecognized token so the outer niimath CLI can continue. Returns 0 on success,
   1 on a malformed sub-option. Records which overrides the user passed in
   `opts->cli_set` so a distinct engine (the fast path) can reject what it cannot honor.
   `caps` is the command's AL_CAP_* capability set: an option outside it is rejected at
   parse time (`-allineate` passes AL_CAP_ALL; `-deface` passes
   AL_CAP_TUNING|AL_CAP_FINAL|AL_CAP_FAST — deface supports the fast engine and defaults to it),
   so a command can never silently accept an option it does not implement. */
static inline int al_parse_subopts(int *ac, int argc, char **argv, al_opts *opts,
                                   const char *cmd_name, unsigned caps) {
    /* Reject an option the command does not implement (parse-time, before images load). */
#define AL_NEED(bit, opt) do { if (!(caps & (bit))) { \
        fprintf(stderr, "%s does not support %s\n", cmd_name, (opt)); return 1; } } while (0)
    while (*ac + 1 < argc && argv[*ac + 1][0] == '-') {
        (*ac)++;
        if (!strcmp(argv[*ac], "-cmass")) {
            AL_NEED(AL_CAP_TUNING, "-cmass");
            opts->cmass = AL_CMASS_YES; opts->cli_set |= AL_CLI_CMASS;
        } else if (!strcmp(argv[*ac], "-nocmass")) {
            AL_NEED(AL_CAP_TUNING, "-nocmass");
            opts->cmass = AL_CMASS_NONE; opts->cli_set |= AL_CLI_CMASS;
        } else if (!strcmp(argv[*ac], "-source_automask")) {
            AL_NEED(AL_CAP_TUNING, "-source_automask");
            opts->source_automask = 1;
        } else if (!strcmp(argv[*ac], "-dark_automask")) {
            AL_NEED(AL_CAP_TUNING, "-dark_automask");
            opts->dark_automask = 1;
        } else if (!strcmp(argv[*ac], "-nearest") || !strcmp(argv[*ac], "-NN")) {
            AL_NEED(AL_CAP_FINAL, "-nearest");
            opts->final_interp = AL_INTERP_NN; opts->cli_set |= AL_CLI_FINAL;
        } else if (!strcmp(argv[*ac], "-linear") || !strcmp(argv[*ac], "-trilinear")) {
            AL_NEED(AL_CAP_FINAL, "-linear");
            opts->final_interp = AL_INTERP_LINEAR; opts->cli_set |= AL_CLI_FINAL;
        } else if (!strcmp(argv[*ac], "-cubic") || !strcmp(argv[*ac], "-tricubic")) {
            AL_NEED(AL_CAP_FINAL, "-cubic");
            opts->final_interp = AL_INTERP_CUBIC; opts->cli_set |= AL_CLI_FINAL;
        } else if (!strcmp(argv[*ac], "-warp")) {
            AL_NEED(AL_CAP_TUNING, "-warp");
            (*ac)++;
            if (*ac >= argc) {
                fprintf(stderr, "%s -warp requires a type (sho, shr, srs, aff)\n", cmd_name);
                return 1;
            }
            if (al_parse_warp(argv[*ac], &opts->warp)) {
                fprintf(stderr, "Unknown warp '%s' (use: sho, shr, srs, aff)\n", argv[*ac]);
                return 1;
            }
            opts->cli_set |= AL_CLI_WARP;
        } else if (!strcmp(argv[*ac], "-interp")) {
            AL_NEED(AL_CAP_TUNING, "-interp");
            (*ac)++;
            if (*ac >= argc) {
                fprintf(stderr, "%s -interp requires an interpolation name\n", cmd_name);
                return 1;
            }
            if (al_parse_interp(argv[*ac], &opts->interp)) {
                fprintf(stderr, "Unknown interp '%s' (use: NN, linear, cubic)\n", argv[*ac]);
                return 1;
            }
            opts->cli_set |= AL_CLI_INTERP;
        } else if (!strcmp(argv[*ac], "-final")) {
            AL_NEED(AL_CAP_FINAL, "-final");
            (*ac)++;
            if (*ac >= argc) {
                fprintf(stderr, "%s -final requires an interpolation name\n", cmd_name);
                return 1;
            }
            if (al_parse_interp(argv[*ac], &opts->final_interp)) {
                fprintf(stderr, "Unknown final interp '%s' (use: NN, linear, cubic)\n", argv[*ac]);
                return 1;
            }
            opts->cli_set |= AL_CLI_FINAL;
        } else if (!strcmp(argv[*ac], "-fill")) {
            AL_NEED(AL_CAP_FILL, "-fill");
            (*ac)++;
            if (*ac >= argc) {
                fprintf(stderr, "%s -fill requires a mode (auto, zero, nan)\n", cmd_name);
                return 1;
            }
            if (al_parse_fill(argv[*ac], &opts->fillmode)) {
                fprintf(stderr, "Unknown fill '%s' (use: auto, zero, nan)\n", argv[*ac]);
                return 1;
            }
        } else if (!strcmp(argv[*ac], "-master")) {
            AL_NEED(AL_CAP_MASTER, "-master");
            (*ac)++;
            if (*ac >= argc) {
                fprintf(stderr, "%s -master requires an output-grid image filename\n", cmd_name);
                return 1;
            }
            opts->master = argv[*ac];
        } else if (!strcmp(argv[*ac], "-weight")) {
            AL_NEED(AL_CAP_WEIGHT, "-weight");
            (*ac)++;
            if (*ac >= argc) {
                fprintf(stderr, "%s -weight requires a weight-image filename (base/fixed-space)\n", cmd_name);
                return 1;
            }
            opts->weight = argv[*ac];
        } else if (!strcmp(argv[*ac], "-savemat")) {
            AL_NEED(AL_CAP_MATRIX, "-savemat");
            (*ac)++;
            if (*ac >= argc) {
                fprintf(stderr, "%s -savemat requires an output filename (.json)\n", cmd_name);
                return 1;
            }
            opts->savemat = argv[*ac];
        } else if (!strcmp(argv[*ac], "-applymat")) {
            AL_NEED(AL_CAP_MATRIX, "-applymat");
            (*ac)++;
            if (*ac >= argc) {
                fprintf(stderr, "%s -applymat requires a matrix filename (.json from -savemat)\n", cmd_name);
                return 1;
            }
            opts->applymat = argv[*ac];
        } else if (!strcmp(argv[*ac], "-com")) {
            AL_NEED(AL_CAP_SEED, "-com");
            opts->com = 1;
        } else if (!strcmp(argv[*ac], "-sym")) {
            AL_NEED(AL_CAP_SEED, "-sym");
            opts->sym = 1; opts->sym_deoblique = 0;   /* last-option-wins: -symd -sym runs plain -sym */
        } else if (!strcmp(argv[*ac], "-symd")) {
            AL_NEED(AL_CAP_SEED, "-symd");
            opts->sym = 1; opts->sym_deoblique = 1;
        } else if (!strcmp(argv[*ac], "-symb")) {
            AL_NEED(AL_CAP_SEED, "-symb");
            opts->sym = 1; opts->sym_deoblique = 2;   /* auto-compete: best of -sym / -symd */
        } else if (!strcmp(argv[*ac], "-nosagseed")) {
            AL_NEED(AL_CAP_SEED, "-nosagseed");
            opts->sagseed = 0;
        } else if (!strcmp(argv[*ac], "-zoom")) {
            AL_NEED(AL_CAP_SEED, "-zoom");
            opts->zoom = 1;
        } else if (!strcmp(argv[*ac], "-cost")) {
            (*ac)++;
            if (*ac >= argc) {
                fprintf(stderr, "%s -cost requires a cost function name\n", cmd_name);
                return 1;
            }
            /* Four names select the FAST engine rather than an allineate cost:
               `fast`/`fastx` = HEL/CR coarse competition followed by HEL fine stages,
               `fasthel` = Hellinger only, `fastcr` = correlation-ratio only.
               `-cost` is last-one-wins: a normal cost after `-cost fast` clears the fast
               engine selection (and vice-versa) so the final `-cost` always decides. */
            if (!strcmp(argv[*ac], "fast")) {
                AL_NEED(AL_CAP_FAST, "-cost fast");
                opts->fast = AL_ENGINE_FAST_X;
                opts->cli_set &= ~AL_CLI_COST;   /* fast overrides an earlier ordinary -cost (symmetric last-one-wins) */
            } else if (!strcmp(argv[*ac], "fasthel")) {
                AL_NEED(AL_CAP_FAST, "-cost fasthel");
                opts->fast = AL_ENGINE_FAST_HEL;
                opts->cli_set &= ~AL_CLI_COST;
            } else if (!strcmp(argv[*ac], "fastcr")) {
                AL_NEED(AL_CAP_FAST, "-cost fastcr");
                opts->fast = AL_ENGINE_FAST_CR;
                opts->cli_set &= ~AL_CLI_COST;
            } else if (!strcmp(argv[*ac], "fastx")) {
                AL_NEED(AL_CAP_FAST, "-cost fastx");
                opts->fast = AL_ENGINE_FAST_X;
                opts->cli_set &= ~AL_CLI_COST;
            } else if (al_parse_cost(argv[*ac], &opts->cost)) {
                fprintf(stderr, "Unknown cost function '%s' (use: fast, fastx, fasthel, fastcr, lpc, lpa, hel, nmi, ls)\n",
                        argv[*ac]);
                return 1;
            } else {
                AL_NEED(AL_CAP_TUNING, "-cost");
                opts->fast = 0;   /* normal cost overrides an earlier fast-engine selector */
                opts->cli_set |= AL_CLI_COST;
            }
        } else {
            /* Not a recognized sub-argument, back up */
            (*ac)--;
            break;
        }
    }
    return 0;
#undef AL_NEED
}

/* The registration/deface API uses process-global and thread-local workspaces.
   It is safe for serial CLI-style use, but not for concurrent calls. */

/* Choose an image's index->world (mm) transform with NIfTI sform/qform precedence
   and validity policy (prefer sform when sform_code >= qform_code, else qform; fall
   back to whichever form is usable; a degenerate/bogus preferred form is skipped).
   Writes *out and returns 0 on success; returns 1 when neither form is usable.
   Pure helper (no global state) shared with the fast coreg path (coreg_fast.c).
   Exposed non-static in Phase 2 of the fast-coreg work with NO logic change —
   exact parity holds (a no-logic-change exposure). */
int al_image_xform(const nifti_image *nim, mat44 *out);

/* Index->world (mm) transform with the single no-form fallback policy: al_image_xform's
   coded sform/qform selection, else a pixdim-centered frame. The one policy used by every
   registration/geometry entry point (al_register, coreg_fast, -sym, -com, -sagseed, apply)
   so both-codes-zero inputs always yield a usable frame. Never fails. `who` non-NULL logs
   the fallback; NULL is silent. */
void al_image_xform_or_pixdim(const nifti_image *nim, mat44 *out, const char *who);

/* Register source image to base image grid using affine (12 DOF) alignment.
   source: the moving image (will be modified in-place: data replaced, dims updated)
   base: the stationary/reference image
   opts: registration options (cost function, cmass, interpolation, etc.)
   final_interp default: cubic.
   Returns 0 on success, nonzero on error. */
int nii_allineate(nifti_image *source, nifti_image *base, al_opts opts);

/* Estimate-only registration: fit `source` to `base` and write the world-mm
   FIXED(base)->MOVING(source) "pull" affine into *fixed_to_moving WITHOUT reslicing or
   mutating either image (the caller applies it, e.g. via nii_apply_affine). This is the
   estimate/apply seam: -master estimates once here then applies once onto its output
   grid (no reslice-then-discard), and -savemat can serialize the matrix directly.
   Returns 0 on success, nonzero on error. Serial-only. */
int nii_allineate_estimate(nifti_image *source, nifti_image *base, al_opts opts,
                           mat44 *fixed_to_moving);

/* Copy the most recent nii_allineate() fit's world-space FIXED(base)->MOVING(source)
   affine (mm) into *out. Returns 0 if a fit has run, nonzero otherwise. The matrix
   is the "pull"/resampling transform (fixed-grid world coords -> moving world coords)
   and reflects the moving image as passed to registration (after any -com/-sym fold).
   Serial-only (reads a process-global, like the rest of the engine). */
int nii_last_affine(mat44 *out);

/* Resolve an AL_FILL_* mode to the concrete out-of-FOV fill value for a float32
 * source of `n` voxels: AL_FILL_ZERO -> 0, AL_FILL_NAN -> NaN, AL_FILL_AUTO -> 0
 * unless the source minimum is negative (CT/HU air), then that minimum. Uses a
 * finite (non-NaN/Inf) guard on the scan. Shared BSD; used by -allineate's output
 * reslice so the fill matches wherever the reslice happens (fused warp or nii_apply_affine). */
float al_resolve_fillv(int fillmode, const float *data, size_t n);

/* al_resolve_fillv for a whole image of ANY datatype (extracts float internally so AUTO
 * detects negative Hounsfield air in a raw int16 CT). Returns 0 on extraction failure. */
float al_image_fillv(int fillmode, nifti_image *nim);

/* Reslice `source` onto `base`'s grid using an explicit base-index -> source-index
 * affine `gam` (0-based NIfTI voxel indices). Replaces source->data with the
 * resliced float volume and adopts base dims + sform/qform. interp: AL_INTERP_*;
 * out-of-FOV voxels set to `fillv` (0 or NaN). Returns 0 on success. Shared BSD
 * interpolation (used by -allineate output and the GPL -spm_coreg reslice). */
int nii_reslice_affine(nifti_image *source, const nifti_image *base,
                       mat44 gam, int interp, float fillv);

/* Apply a saved world-space FIXED->MOVING affine (`-savemat` `fixed_to_moving`) to
   reslice `input` (an image in the prior registration's MOVING space) onto `target`'s
   grid — `target` may have any resolution/FOV/origin sharing the fixed world frame
   (the matrix is world-mm, not grid-specific). input->data is replaced with the
   resliced volume on target's grid. interp: AL_INTERP_*; out-of-FOV -> fillv.
   Returns 0 on success, nonzero on error. */
int nii_apply_affine(nifti_image *input, const nifti_image *target,
                     mat44 fixed_to_moving, int interp, float fillv);

/* Template-free midsagittal-plane (MSP) alignment: register `nim` to its world-X
   mirror, take the half of the recovered rigid transform, and either resample the
   data symmetric about world X=0 (reslice != 0, standalone use) or fold the
   correction into the header as a registration seed (reslice == 0, pre-step use).
   If `C_out` is non-NULL it receives the 4x4 world-space correction. The mirror fit
   uses the `ls`/Pearson cost. `deoblique` (nonzero, `-symd`) first snaps the frame
   to axis-aligned (treats the voxel grid as anatomical) so an obliquely-acquired but
   grid-symmetric head is not rotated onto the oblique world frame. `deoblique`:
   0 = `-sym` (image's world frame); 1 = `-symd` (snap the frame axis-aligned first);
   2 = `-symb` (fit both frames; keep the de-obliqued frame only if its correction
   rotation is smaller by more than AL_SYMB_ROT_TOL_DEG, else keep the original world
   frame — a tie/dead-band bias, and only rotation is compared, not fit cost).
   `dark_automask` (nonzero, `-dark_automask`) drops background/pad matched pairs (at
   the image minimum) from the mirror-fit cost. Uses the process-global registration
   workspaces, so it is serial-only (see nii_allineate).
   Returns 0 on success, nonzero on error. */
int nii_symmetry(nifti_image *nim, mat44 *C_out, int reslice, int deoblique, int dark_automask);

/* -sagseed: in-MSP rigid seed, the complement of -sym. Precondition: -sym has
   folded its MSP correction into `nim`'s header (midsagittal plane at world X=0)
   and `tmpl` is an MSP-aligned template. Runs a 3-DOF-constrained fit freeing only
   the MSP-preserving isometries {y-shift, z-shift, pitch} — the rigid DOF -sym is
   blind to — and folds the correction into nim's header as a full-rigid seed for a
   subsequent nii_allineate. Uses the user's cost/interp/source_automask/cmass from
   `opts`. Serial-only (process-global workspaces). Returns 0 on success. */
int nii_sagseed(nifti_image *nim, nifti_image *tmpl, al_opts opts);

/* -com: set the image origin to its brightness center of mass. Computes the
   intensity-weighted centroid (positive voxels), maps it to world coordinates via
   the selected index->world transform, and folds a pure translation into the
   header so the centroid sits at world (0,0,0). Header-only (no reslice); a cheap
   origin reset, intended to run early (right after -robustfov). Template-free:
   falls back to a pixdim-centered frame when the input has no usable sform/qform.
   Returns 0 on success, nonzero on error. */
int nii_center_of_mass(nifti_image *nim);

/* Deface: register INPUT to TEMPLATE (the well-posed direction, base =
   template, same as -allineate), INVERT the transform, warp the template-space mask
   onto the input's native grid, and set voxels where warped mask < 0.5 to the input's
   minimum value. (Registering template->input instead mislocates the mask.)
   input: the image to modify (modified in-place, stays in its own native space)
   tmpl: template image — the registration base/fixed (input is the moving image)
   mask: mask in template space (>=0.5 = keep, <0.5 = remove). CONSUMED: its data is
         freed and replaced with the resliced mask on the input grid (regridded to
         input geometry). Caller still owns the nifti_image and must free it.
   opts: registration options (cost function, cmass)
   final_interp default: linear (to avoid ringing in the mask).
   Returns 0 on success, nonzero on error. */
int nii_deface(nifti_image *input, nifti_image *tmpl, nifti_image *mask, al_opts opts);

/* Zero (to image minimum) input voxels where the input-grid `warped_mask` < 0.5.
   Ensures input is float32. BSD; shared by -deface and the GPL -spm_deface.
   Returns #voxels masked, or -1 on error. */
long nii_apply_deface_mask(nifti_image *input, const float *warped_mask);

#endif /* ALLINEATE_H */
