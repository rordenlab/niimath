#ifndef COREG_FAST_H
#define COREG_FAST_H

/* Fast affine coregistration path (`-cost fast`/`-cost fastx`/`-cost fasthel`/
 * `-cost fastcr`) — an SPM/FLIRT-inspired, but
 * independently implemented, multiresolution 12-DOF affine estimator for 3D adult
 * human brain images. It is a distinct registration tier — selected via `-cost fast`,
 * which is the CLI default, but a SEPARATE estimator that does NOT change nii_allineate()'s
 * API or numerical behavior: nii_allineate() remains the ordinary AFNI-style engine, reached
 * with `-cost hel`/`lpc`/`lpa`/`ls`. (The CLI default lives in the dispatch, not here.)
 *
 * Contract:
 *  - coreg_fast_estimate() does NOT mutate its inputs. It returns a world-mm
 *    FIXED(base)->MOVING(source) "pull" transform (same convention as -savemat's
 *    fixed_to_moving). The caller applies it once (nii_apply_affine / a final
 *    reslice) and writes -savemat directly from the result.
 *  - Serial-only at the host-call level (OpenMP is used only inside owned loops).
 *
 * See AGENTS.md for the release contract and clean-room provenance.
 */

#include "nifti_io.h"

/* Cost functions for the fast path. */
#define CF_COST_LS   0   /* weighted Pearson/LS (diagnostic; same-modality) */
#define CF_COST_CR   1   /* directional Correlation Ratio ('-cost fastcr'): UNSTABLE
                            when modalities differ (spurious shrink/roll on cross-modal). */
#define CF_COST_HEL  2   /* Hellinger/MI family ('-cost fasthel'). Robust on
                            real cross-modal images (T2w->T1) and comparably fast. Only caveat: MI
                            is imprecise on smooth low-detail volumes (it can't sharply recover the
                            §7 synthetic Gaussian phantoms) — a property of MI on smooth data, not
                            a defect on real scans; the synthetic capture suites use '-cost fastcr'. */
#define CF_COST_HEL_CR 3 /* default / '-cost fast' / '-cost fastx': independently fit the 8 mm
                            rigid coarse stage from supplied-affine/COM starts with HEL and CR.
                            Whole-head bases select by HEL dependence*overlap and continue once.
                            Hard-zeroed bases additionally carry scale-bracketed HEL and CR-seeded
                            strategies through the 2 mm 7-DOF arbiter, then polish only the winner.
                            This is a mode selector, never passed to cf_cost_eval directly. */

typedef struct {
    int cost;          /* CF_COST_* (default CF_COST_HEL_CR) */
    int coarse_search; /* 1 = bounded deterministic orientation/scale search after
                          initialization selection (default); 0 = local descent only */
    int max_dof;       /* highest DOF to fit: 6/7/9/12 (default 12) */
    int use_cmass;     /* 1 = auto-select supplied-affine vs COM-recentered initialization
                          from initial dependence*overlap (default); 0 = use the
                          supplied frame only (-nocmass, or an already-applied -com) */
    int verbose;       /* nonzero: per-stage logging to stderr */
    const nifti_image *weight; /* optional (-weight): an AFNI 3dAllineate-style graded weight in
                          fixed(base) space with dims identical to `fixed`. Normalized to [0,1] over
                          the whole grid (÷ max) and applied GRADED per fixed sample (a voxel weighted
                          0 is excluded, ~1 dominates) — NOT the old soft-focus floor. Applied ONLY at
                          the FINEST pyramid level (2 mm) so the 8 mm rigid coarse and the 4 mm
                          global-scale bracket stay whole-head (a brain-concentrated weight must not
                          re-select global scale — the FLIRT Fig-1 shrink basin). Keep the out-of-ROI
                          head attenuated (nonzero) to anchor scale. Honored by the HEL and CR costs;
                          rejected with CF_COST_LS (which does not read it). NULL (default) leaves the
                          fit byte-for-byte unchanged. Not owned by the estimator (caller frees). */
    /* NOTE: the final reslice/interpolation is applied by the CALLER (main.c), not the
       estimator, so there is no interp field here — coreg_fast_estimate() only returns
       the world affine and does not touch image data. */
} coreg_fast_opts;

typedef struct {
    mat44  fixed_to_moving;  /* world-mm FIXED->MOVING pull transform */
    double final_cost;       /* cost at the returned transform (lower is better) */
    int    evaluations;      /* total cost evaluations across all stages */
    int    levels_completed; /* number of pyramid levels optimized */
    double registration_ms;  /* estimation wall time (excludes the final reslice) */
    int    resolved_cost;    /* CF_COST_* actually used after validation */
    int    resolved_dof;     /* highest DOF actually fitted (6/7/9/12) */
} coreg_fast_result;

static inline coreg_fast_opts coreg_fast_opts_default(void) {
    coreg_fast_opts o;
    o.cost = CF_COST_HEL_CR; /* Mixed capture is the default; hosts map '-cost fasthel' to HEL
                                and '-cost fastcr' to CR explicitly. */
    o.coarse_search = 1;
    o.max_dof = 12;
    o.use_cmass = 1;   /* auto-select supplied-affine vs COM initialization */
    o.verbose = 0;
    o.weight = NULL;   /* no region weighting unless -weight supplies a fixed-space image */
    return o;
}

/* Estimate the fast-path affine. Inputs are not modified. On success writes
 * *result and returns 0. On failure returns nonzero and leaves *result unchanged.
 * Serial-only at the host-call level (fastx uses callback TLS for its internal coarse jobs). */
int coreg_fast_estimate(const nifti_image *moving, const nifti_image *fixed,
                        const coreg_fast_opts *opts, coreg_fast_result *result);

#endif /* COREG_FAST_H */
