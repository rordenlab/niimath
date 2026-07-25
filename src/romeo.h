// romeo.h - ROMEO phase unwrapping for niimath (-romeo)
//
// Faithful C port of ROMEO.jl (https://github.com/korbinian90/ROMEO.jl) plus the
// MriResearchTools.jl helpers its command-line app depends on. Both upstream projects are
// MIT-licensed; see src/romeo.LICENSE for the preserved copyright/permission notices and
// romeo.c for the pinned reference versions and the numeric-type audit.
//
// Guarded by HAVE_ROMEO. romeo.c MUST be compiled with strict floating point
// (-fno-fast-math -ffp-contract=off): ROMEO's growth order is decided by 8-bit integer edge
// weights, so a single reassociated/contracted expression can move a weight across a bin
// boundary and shift a whole connected region by 2*pi.

#ifndef ROMEO_H
#define ROMEO_H

#ifdef __cplusplus
extern "C" {
#endif

#include "core.h"

#define ROMEO_MAX_TE 1024

// -w selection
enum {
	RM_W_ROMEO = 0, // resolved to romeo3 (with magnitude) or romeo4 (without)
	RM_W_ROMEO2,
	RM_W_ROMEO3,
	RM_W_ROMEO4,
	RM_W_ROMEO6,
	RM_W_FLAGS // explicit bit flags in opts.flags[]
};

// -k selection
enum {
	RM_MASK_ROBUST = 0, // default
	RM_MASK_NONE,
	RM_MASK_QUALITY,
	RM_MASK_FILE
};

typedef struct {
	int weights_sel;         // RM_W_*
	int flags[6];            // used when weights_sel == RM_W_FLAGS
	double TEs[ROMEO_MAX_TE];
	int nTE;                 // 0 = not supplied
	int te_epi;              // -t epi (all echoes share one TE)
	int mask_sel;            // RM_MASK_*
	const char *mask_file;   // RM_MASK_FILE
	double qmask_thresh;     // -k qualitymask <thr>
	int template_echo;       // 1-based, default 1
	int correctglobal;       // -g
	int individual;          // -i
	int verbose;             // -v
	int write_quality;       // -q
	int write_quality_all;   // -Q
	int no_mask_out;         // niimath-only: suppress the <base>_mask side output
	int no_phase_rescale;    // -no-phase-rescale
	int maxseeds;            // -max-seeds (only 1 supported)
	double temporal_uncertain; // -temporal-uncertain-unwrapping
	double wrap_addition;    // -wrap-addition (only 0 supported)
	const char *dumpdir;     // hidden -romeo-dump <dir>: raw parity dumps for test/romeo_compare.py
} romeo_opts;

romeo_opts romeo_opts_default(void);

// Consume recognized ROMEO sub-options starting at argv[*pac]; stops (without consuming) at the
// first unrecognized token so later niimath chain operations stay visible. Returns 0 on success.
int romeo_parse_subopts(int *pac, int argc, char *argv[], romeo_opts *o, const char *cmd);

// Run ROMEO on the working image. `nim` must be float32 and 3D or 4D (echoes on dim 4).
// `magfile` is the magnitude filename or NULL (the literal "none" resolves to NULL upstream).
// `phasefile` is the ORIGINAL input filename, needed by the readphase rescale branch that
// inspects unscaled stored values; pass NULL when unavailable (stdin), which makes a required
// rescale fail loudly rather than silently.  `is_first_op` records whether -romeo is the first
// computational operation (phase rescaling requires it).  `ihdr` is currently UNUSED: it is kept
// because the deferred `-fix-ge-phase` (plan M7) needs the ORIGINAL stored datatype to choose
// between the integer and float branches of MriResearchTools' fix_ge_phase!.  On success
// nim->data holds the unwrapped phase and side outputs have been written.  Returns 0 on success.
int romeo_run(nifti_image *nim, const char *magfile, const char *phasefile,
	const in_hdr *ihdr, int is_first_op, const romeo_opts *o, gzModes gzMode);

#ifdef __cplusplus
}
#endif

#endif // ROMEO_H
