# Plan: port ROMEO phase unwrapping to niimath (`-romeo`, `romeo.c`, `HAVE_ROMEO`)

## STATUS (2026-07-25)

**M0–M6, M8 and M9 are DONE and merged** (`romeo` branch). `-romeo` matches the real ROMEO CLI app (RomeoApp 4.5) on all three supplied validation cases: unwrapped phase EQUAL at `--compare 1e-7`, mask EQUAL at `--compare 0`. 602/602 oracle checks pass across 4 real + 11 synthetic cases, 10 weight selections and 6 B0 weighting modes, including the four unwrapping variants, the seed scalars, a ≤4 ULP gate on the pre-rescale weights, a byte-exact numeric-primitive table (the only coverage of the Payne-Hanek branch of `rem2pi`), and ROMEO.jl's own `test/features.jl` property. Thread parity holds at `-p 1/2/4/8`; UBSan and `leaks` are clean; the WASI reactor agrees with the native binary exactly.

**M8 is DONE.** `-B`/`-B0-phase-weighting` compute a B0 field map in Hz with all six weighting modes and the `exp(-TE/20)` synthetic-magnitude fallback, verified against the pinned Julia for every mode with and without a magnitude (`602/602` oracle checks). Two documented deviations: B0 is computed WITHOUT MCPC-3D-S (equivalent to `romeo --compute-B0 --phase-offset-correction off`, announced on stderr rather than silently implied), and where upstream's magnitude-free SNR map collapses to a single value niimath writes that same constant across the working grid.

**M7 remains DEFERRED.** Every option it covers is rejected at parse time with a specific message rather than silently ignored: `-u`, `-e`, `-threshold`, `-w bestpath`, `-max-seeds > 1`, `-merge-regions`, `-correct-regions`, `-wrap-addition != 0`, `-fix-ge-phase`. The experimental multi-seed/region-merging group is the part the plan itself flags as possibly not byte-reproducible (Julia `Dict` iteration order); the rest is small if anyone asks.

Three corrections to the plan as written, all confirmed empirically and recorded in `src/romeo.c`:

1. **§3.2 / quantile.** The pinned environment resolves the *registry* package `Statistics` v1.11.1, whose `_quantile` computes `aleph = n*p + m`. The copy bundled in Julia's stdlib tree (which the plan was written against) uses `fma(n, p, m)`. They differ in the last bits and move `maxmag` in the 11th digit.
2. **§3.2 / `rem2pi`.** libc `remainder(x, 2π_double)` is not merely "not specified to match" — it is *systematically* different, because Julia's `rem2pi` reduces against an infinitely precise 2π. Julia's `rem_pio2_kernel` (Cody-Waite + Payne-Hanek) is ported literally, with a portable 128-bit emulation so MSVC and wasm produce the same reduction.
3. **§3.3 / strict FP.** The plan asked for the exposure to be measured. It was, three ways: `-ffp-contract=off` (shipped) passes the full parity suite 422/422; FMA contraction is bit-identical on the validation volume but fails 9 checks on the full corpus; the repository-wide `-ffast-math` changes 360/797088 weight bytes and leaves 66 voxels off by a full 2π. The mechanism is not rounding — reassociation pushes a weight past 1.0 and `rescale()`'s `0 ≤ w ≤ 1` guard then returns bin 0, deleting the edge from the graph. Compute is 0.02 s under either policy (0.07 s wall including gzip output), so contraction buys nothing. **Any future FP-policy change must be re-measured on the full suite, not on one volume** — that single-volume shortcut is exactly what produced a wrong "FMA is clean" claim in the first place.

## Audit outcome and residual items (for the next session)

A three-agent audit (security/bugs, refactor, docs) ran over the branch. Everything it found that mattered is FIXED and committed; what follows is what was deliberately left.

**Fixed:** an MSVC build break (`strtok_r`, which would have failed the Windows release job); a `-template` `long`→`int` truncation giving SIGSEGV/SIGBUS and a silent wrong-echo run; a heap-use-after-free in `rm_read_f32`'s unsupported-datatype error path; undefined behaviour in the Payne-Hanek reduction (`idx << 6` on a negative int); `rm_pq_enqueue` failures being dropped (fail-open OOM → a region silently left wrapped, exit 0); unchecked `snprintf` truncation in `-romeo-dump` that could alias two dumps onto one filename; `-v` accepted but doing nothing; and several documentation errors, including one claim of mine that was flatly wrong (see §3 above).

**Deliberately NOT done, with reasons:**

- Hoisting the temporal-uncertain scratch buffers out of their loop (three duplicated free-blocks). Mechanical, low risk, but touches the one loop where a missed `free` becomes a leak — worth doing only alongside another edit in that function.
- Replacing the six per-target `$(if $(ROMEOOBJS),…)` Makefile lines with a `define`/`$(call)`. The audit judged the cure a wash: still six lines, plus a `$(call)`-inside-a-recipe escaping gotcha, against a Makefile whose house style is explicit-and-repetitive throughout.
- `rm_sample_capacity` over-allocates (`max(len², n)` because of the rare all-non-finite fallback), so `rm_robustmask` transiently costs ~134 MB on a 256³ magnitude. Pure allocation change, no effect on values; only worth doing if someone hits it.
- Unifying `rm_select_kth_f`/`rm_select_kth_d` behind a macro: numerically safe but saves 16 lines at the cost of macro-obscured debugging in the file you most want to read literally.

**Verified clean, so do not re-audit speculatively:** the 1-based↔0-based index conversions (ASan over 104 volume shapes × 4 option sets), leaks/double-free on all `goto done` paths (0 leaks over 6 success + 9 error paths), NaN/Inf reaching an index or allocation size, quickselect on non-finite input (400k fuzz trials), the `coreFLT.c` dispatch, and the 128-bit Payne-Hanek emulation — validated against 3000-bit mpmath over ~3000 values spanning the full double exponent range, worst error 0.4998 ULP.

One item in §3.1 is worth re-reading before touching the weights: `phaselinearity`'s return type is *data-dependent* (Float32 for the interior product, Float64 on its `isnan → 0.5` and boundary-`0.9` branches), and `unwrapedge!`'s `d = 0` is an **Int**, which silently selects a Float32 rather than Float64 subtraction inside `unwrapvoxel`.


## 1. Goal

Add an optional, self-contained C module (`src/romeo.c` / `src/romeo.h`, guarded by `HAVE_ROMEO`) that is a **faithful port** of [ROMEO.jl](https://github.com/korbinian90/ROMEO.jl) plus the small set of `MriResearchTools.jl` functions its command-line app depends on (`robustmask`, `readphase` rescaling, `gaussiansmooth3d` box filtering, `calculateB0_unwrapped`). Both upstream projects are MIT-licensed.

Target CLI (ordinary chain op, so further niimath ops may follow):

```bash
niimath phase0.nii.gz -romeo mag0.nii.gz -t 16.8 out
niimath phase.nii    -romeo mag.nii     -t '[16.8,38.56]' out
```

`out.nii.gz` is the unwrapped phase; with the default `robustmask`, `out_mask.nii.gz` is the mask that constrained unwrapping. This mirrors the Julia reference calls that produced `/Users/chris/src/ROMEO.jl/romeo/results{,0,1}/`.

**Fidelity is the primary requirement.** ROMEO is a minimum-spanning-tree region-growing algorithm whose growth order is decided by 8-bit integer edge weights. A one-ULP float difference that flips a weight from bin 137 to bin 138 can change the MST traversal order and shift an entire connected region by exactly 2π. The plan therefore requires exact integer intermediates plus a narrow end-to-end tolerance that is far below one wrap.

The minimum definition of done is the default Julia CLI path used by the three supplied validation cases: phase rescaling, default weight resolution, `robustmask`, single-echo spatial unwrapping, multi-echo temporal unwrapping, mask output, optional compilation, and native/WASM integration. The less common CLI paths (bestpath, individual/uncertain unwrapping, multiple regions, B0) remain required for full app compatibility, but are deliberately scheduled after that core parity gate.

## 2. What the Julia code actually is

Total surface to port is small — ~1,225 lines of Julia across ROMEO.jl, plus ~250 relevant lines of MriResearchTools. Expect roughly 1,800–2,400 lines of C.

### 2.1 Reference version (pin before writing an oracle)

Parity means parity with one immutable environment, not with whatever `Pkg` resolves later. The environment inspected for this plan is:

- Julia **1.12.3**
- ROMEO **1.4.0**, git commit `60d83fbb69669560d227c252dcd844afbb6648e1`
- MriResearchTools **3.5.0**, manifest tree `8ad52837d6686e0f8c0f2f7be45c36f3a98cb000`
- ImageMorphology **0.4.6**, NIfTI **0.6.2**, StatsBase **0.34.12**
- `/Users/chris/src/ROMEO.jl/Manifest.toml` SHA-256 `f6718247536608251ddef803f9674656cdb20b3141e204f250b90e6a6c368065`

The checkout currently has a modified `Project.toml`, so the commit alone does not reconstruct the app environment. M0 must copy the reference `Project.toml` and `Manifest.toml` (or record their hashes and refuse a mismatch) before generating golden data. Every oracle artifact records the Julia version, ROMEO commit, manifest hash, command, and input SHA-256.

| Julia source | Contents | C destination |
| --- | --- | --- |
| `src/utility.jl` | `γ` (single-wrap fold), `getdimoffsets` | `romeo.c` static inlines |
| `src/priorityqueue.jl` | `PQueue`: bucket queue, LIFO within bin | `romeo.c` (`rq_*`) |
| `src/weights.jl` | 6 weight terms, `calculateweights_romeo`, `rescale`, bestpath weights | `romeo.c` (`rm_weights_*`) |
| `src/seed.jl` | seed queue, `findseed!`, `seedcorrection!` | `romeo.c` (`rm_seed_*`) |
| `src/algorithm.jl` | `grow_region_unwrap!` MST loop, edge indexing, `unwrapvoxel` | `romeo.c` (`rm_grow_region`) |
| `src/region_handling.jl` | `merge_regions!`, `correct_regions!` (experimental, `maxseeds>1`) | `romeo.c` (`rm_merge_regions`) |
| `src/unwrapping.jl` | 3D/4D `unwrap!`, temporal unwrapping, `temporal_uncertain_unwrapping!`, `unwrap_individual!` | `romeo.c` (`rm_unwrap3d`, `rm_unwrap4d`) |
| `src/voxelquality.jl` | `voxelquality` quality map | `romeo.c` (`rm_voxelquality`) |
| `ext/RomeoApp/caller.jl` | app orchestration, mask selection, thresholding, B0 | `romeo.c` (`nii_romeo_run`) |
| `ext/RomeoApp/argparse.jl` | option defaults and `romeo`→`romeo3`/`romeo4` resolution | `coreFLT.c` dispatch + `romeo_parse_subopts` |
| MRT `masking.jl` | `robustmask`, `fill_holes` (`imfill`, 6-connectivity) | `romeo.c` (`rm_robustmask`) |
| MRT `smoothing.jl` | `gaussiansmooth3d` / `boxfilterline!` (running average) | `romeo.c` (`rm_boxfilter_line`) |
| MRT `utility.jl` | `sample` (block subsample), `approxextrema`, `estimatequantile` | `romeo.c` (`rm_sample`, `rm_quantile`) |
| MRT `niftihandling.jl` | `readphase` rescale-to-[-π,π] logic, `fix_ge_phase!` | `romeo.c` (`rm_phase_rescale`) |
| MRT `romeofunctions.jl` | `calculateB0_unwrapped`, `get_B0_snr` | `romeo.c` (`rm_compute_b0`) |

### 2.2 What we can reuse from niimath, and what we deliberately cannot

Reusable as-is: `nifti_image_read`, `set_input_hdr` / `nifti_image_change_datatype`, `nifti_save`, `nii_calloc` / `nii_malloc` / `nii_mul_size` (the `nim->data` ownership invariant applies), and `max_displacement_mm` for phase-vs-auxiliary world-frame sanity **warnings**. The `coreFLT.c` wrapper can call its local `nii_reject_oversize_aux` before loading magnitude/mask/weights; do not make `romeo.c` depend on that static helper. Copy the shape of `al_parse_subopts`'s "stop at first unrecognized token" parser pattern, not the allineate dependency itself.

**Not reusable, despite superficial similarity — this is a correctness trap, not a style preference:**

- `nifti_smooth_gauss` (coreFLT.c) is a true truncated-Gaussian convolution. `robustmask` needs MriResearchTools' `gaussiansmooth3d(mask; nbox=1, boxsizes=[[5],[5],[5]])`, which is *n* passes of a **running-average box filter** with a specific asymmetric edge normalization (`line[i] = lsum / (r + i)` at the leading edge). Substituting the Gaussian changes the mask.
- `nifti_robust_range` / `nifti_otsu` are 1000-bin / 256-bin histogram approximations. `robustmask`'s threshold uses exact type-7 quantiles (0.05, 0.15, 0.8, 0.99) over a **deterministically subsampled** vector. Different numbers, different mask.
- `bwlabel.c` gives connected components but is coupled to the optional mesh inventory and does not directly implement `imfill`'s "turn true components whose inclusive size is in `[1, length/20]` false" semantics. A local 6-connected component pass keeps ROMEO independent of mesh and is easier to validate exactly.

Wherever the Julia and the niimath primitive genuinely agree (NIfTI I/O, allocation, oversize gating), reuse niimath. Where they differ numerically, port the Julia.

## 3. Fidelity hazards (the part that decides whether this works)

These were verified against the pinned local environment in §2.1, not assumed.

### 3.1 Mixed Float32/Float64 promotion

Julia's promotion rules are not uniform across the six weight terms. Verified in the REPL against this exact code:

| Expression | Julia result type |
| --- | --- |
| `phasecoherence(P,i,j) = 1 - abs(γ(P[i]-P[j])/π)` | **Float32** (Irrational π promotes to Float32) |
| `phasegradientcoherence(...)` | **Float64** (`TEs` are Float64) |
| `phaselinearity(P,i,j,k)` | **Float32** normally; **Float64** on its `isnan → 0.5` branch |
| `phaselinearity(P,i,j)` | **Float32** for the interior product; **Float64** for the boundary fallback `0.9` (or a widened NaN fallback) |
| `magcoherence(small,big) = (small/big)^2` | **Float32** |
| `magweight`, `magweight2` | **Float64** |
| `weight` accumulator in `getweight` | **Float64** (`weight = 1.0`) |
| `0.1 + 0.9x` factors | **Float64**, with a Float32 `x` widened |
| `maxmag = quantile(mag[isfinite], 0.95)` | **Float64** even for a Float32 input array |
| `unwrapvoxel(new, old) = new - 2pi*round((new-old)/2pi)` | **Float64** intermediate, stored back into a Float32 array |

So `romeo.c` must compute each term in the matching width (`float` vs `double`) and widen at exactly the same points. Getting this wrong is the single most likely cause of a near-miss. **Action: a "numeric type audit" table is a checked deliverable of M2** — one row per expression, Julia type confirmed by `typeof(...)` in the REPL, C type in the port.

### 3.2 Rounding mode

Julia's `round(Int, x)` is **round-half-to-even** (confirmed: `0.5→0`, `1.5→2`, `2.5→2`). C's `lround`/`roundf` are half-away-from-zero. Use `nearbyint`/`rint` with the default `FE_TONEAREST`, never `round`/`lround`. This affects `rescale()` (which sets the 8-bit weight bin) and `unwrapvoxel()` (which sets the number of 2π wraps) — both load-bearing.

`rem2pi(x, RoundNearest)` in Julia does specialized range reduction; C's `remainder(x, 2π_double)` is not specified to match it bit-for-bit. M2 must exhaustively compare the actual arguments seen in validation and add synthetic values around ±π and half-integer wrap boundaries. Do not assume the arguments stay below 4π: `phaselinearity` and already-unwrapped input can exceed that range.

### 3.3 Whole-program `-ffast-math` (the biggest build-level risk)

`src/Makefile:140` does `CFLAGS += $(AFASTMATH)`; `src/CMakeLists.txt:186` and `notarize.sh` and the wasm targets do the same. `-ffast-math` permits reassociation, flushes denormals, and enables `-ffp-contract=fast` (FMA). Any of these can move a weight across a `rescale` bin boundary.

**Decision: `romeo.c` is compiled as a prebuilt object with strict FP** — `-fno-fast-math -ffp-contract=off` (and `-fno-lto` if link-time optimization is shown to erase the per-TU contract). These flags must occur **after** the repository-wide `CFLAGS`, which contain `-ffast-math`. It follows the existing prebuilt-object pattern but has its **own** `.romeo_obj_flags` signature: ROMEO must rebuild correctly when `CNAME`, optimization/debug/sanitizer flags, `OMP`, or `ROMEO` changes even when `AL=0`.

```make
ROMEOFLAGS= -DHAVE_ROMEO
ROMEOOBJS= romeo.o
ROMEO_STRICT_FP= -fno-fast-math -ffp-contract=off
```

M2 includes an objective measurement of how much this actually matters (build both ways, count differing weight bytes). If the fast-math build turns out to be byte-identical on the validation data we still keep the strict build — the whole point is that the guarantee should not depend on a codegen accident — but we will know the real exposure.

### 3.4 Priority-queue tie-breaking is part of the algorithm

`PQueue` is a 256-bucket queue. `enqueue!` does `push!` (append) and `dequeue!` does `pop!` (**take from the end**) — i.e. **LIFO within a bin**. `q.min` is lowered by `enqueue!` and advanced past empty bins by `dequeue!`. Any C implementation that uses FIFO per bin, or a heap, will produce a different (still valid, but different) spanning tree and hence different 2π assignments in ambiguous regions. Port the stack semantics literally.

Likewise the seed queue is built by `enumerate(sum(weights; dims=1))` in ascending linear index, so within a bin the **highest linear index is dequeued first**. Replicate.

### 3.5 Linear-index bounds checks and the zeroed last planes

Julia's `checkbounds(Bool, A, i)` on a linear index only checks `1 ≤ i ≤ length(A)`, so a neighbour lookup can cross a row or plane while remaining linearly in range. ROMEO relies on `calculateweights` zeroing `weights[1,end,:,:]`, `weights[2,:,end,:]`, `weights[3,:,:,end]` so invalid directed edges never enter the queue. The same linear-index behavior also appears in `phaselinearity` (`h`/`k`) and `unwrapedge!` (`oo`), where it is not equivalent to a fresh Cartesian bounds check. Port the linear formulas literally; add tiny non-cubic fixtures that expose every boundary.

### 3.6 The mask leaks by one voxel, on purpose

In `calculateweights_romeo` the guard is `if mask[I] && checkbounds(Bool, wrapped, J)` — only the **left** voxel of each positive-axis edge must be in the mask. A voxel immediately outside a positive-axis mask face can therefore be reached and unwrapped; it cannot continue growing because its outgoing weights are zero. This directional one-voxel fringe is observable. Replicate; do not tighten to `mask[I] && mask[J]`.

### 3.7 Deterministic subsampling

`sample(I; n=1e5)` takes `len = ceil(sqrt(n))` blocks of `len` contiguous elements, with block starts at `round.(Int, range(firstindex-1, lastindex-len; length=len))`, then filters non-finite, and falls back to filtering the whole array if the result is empty. The validation volumes are 76×76×46 = 265,696 voxels > 1e5, so **subsampling is active** for `approxextrema` and for every `robustmask` quantile. This is deterministic but idiosyncratic; port `range`'s endpoint-inclusive linear spacing and Julia's round-half-to-even exactly.

`quantile` is Julia's default **type 7** (linear interpolation between order statistics `h = (n-1)p`). Implement with two comparator-free quickselects (the WASM `qsort` prohibition applies) plus interpolation — not a sort, not a histogram.

### 3.8 Iteration and reduction order

Julia is column-major, which matches NIfTI's x-fastest voxel layout here, but only if C loops preserve Julia's `LinearIndices` order. This matters for queue insertion, `sum(weights; dims=1)`, mask statistics, medians, and all tie cases. Bestpath also constructs its unique linear neighbor-offset list from `Iterators.product(-1:1,-1:1,-1:1)` in a defined first-occurrence order; port that order rather than replacing it with a geometrically equivalent set.

External masks are read by the Julia app as `niread(file).raw .!= 0`, so mask truth is based on stored voxels and ignores `scl_slope`/`scl_inter`. Magnitude and phase, in contrast, use scaled values. Add a scaled-mask fixture to keep this distinction explicit.

### 3.9 Defaults that differ from the Julia library defaults

The CLI (`argparse.jl`) overrides several library defaults, and the validation runs used the CLI. From `results0/settings_romeo.txt`:

- `weights: romeo` resolves to **`romeo3`** when a magnitude is supplied (flags 1,2,4 = phasecoherence, phasegradientcoherence, magcoherence) and to **`romeo4`** (flags 1–4) when it is not.
- `updateflags` then disables flags 4–6 with no magnitude, and disables flag 2 with no `phase2`/`TEs` pair. **So the single-echo-with-magnitude validation case runs with only flags 1 and 4 active.**
- `temporal-uncertain-unwrapping` defaults to **0** from the CLI (the library default is 0.5).
- `mask` defaults to **`robustmask`**, so the reference runs *did* mask, and *did* write `mask.nii`.
- `max-seeds=1`, so `merge_regions`/`correct_regions` and the multi-seed threshold logic are dead code at default settings — real, but not on the critical path.

## 4. Design

### 4.1 Files

- `src/romeo.h` — public surface: an options struct, `romeo_opts_default()`, the sub-option parser, and one runner taking the working `nifti_image`, optional magnitude filename, original phase metadata needed by `readphase`, and `gzModes`. Keep ownership explicit: success replaces only `nim->data`/shape fields that change; failure writes no main output and leaves cleanup safe.
- `src/romeo.c` — everything else, all `static`. No dependency on `coreFLT.c` internals; only `core.h`, `nifti_io.h`, `<math.h>`.
- Dispatch hook in `coreFLT.c` next to `-reface`/`-qwarp`, wrapped in `#ifdef HAVE_ROMEO` and `#ifdef DT32` (float32 only, like `-reface`; `-dt double` prints an explicit refusal). In the current core, `fin` is the **input** filename; the output filename is already installed in `nim->fname`. Derive side outputs with `nifti_save(side_image, "_mask", gzMode)` (and analogous postfixes), not from `fin`.
- A universal `#else` dispatch stub recognizes `-romeo` in builds without `HAVE_ROMEO` and reports how to enable it, rather than falling through to “unsupported operation.”
- Help text in `niimath.c`.

### 4.2 CLI

```
niimath <phase> -romeo <mag|none> [options] <out>
```

The magnitude is a **required positional token**; pass the literal `none` for magnitude-free unwrapping. (This avoids the ambiguity of an optional positional followed by dashed options.) The main parser has already removed the trailing output from the operation range. A local `romeo_parse_subopts` consumes recognized ROMEO options and stops, backing up, at the first unrecognized token so later niimath chain operations remain visible.

| Option | ROMEO equivalent | Notes |
| --- | --- | --- |
| `-t <TEs>` | `--echo-times` | accepts `16.8`, `16.8,38.56`, `[16.8,38.56]`, `epi`, `epi 5.3`. Bracket form needs shell quoting; document it. |
| `-k <spec> [thr]` | `--mask` | `nomask` \| `robustmask` \| `qualitymask [thr]` \| `<file>` |
| `-u` | `--mask-unwrapped` | |
| `-w <spec>` | `--weights` | `romeo`\|`romeo2`\|`romeo3`\|`romeo4`\|`romeo6`\|`bestpath`\|bit flags e.g. `1010`\|external 4D weights file |
| `-e <spec>` | `--unwrap-echoes` | integer list / `:`; **no `eval`** — accept `:`, `n`, `a,b,c`, `a:b`, `a:s:b` only |
| `-g`, `-q`, `-Q`, `-i`, `-v` | `--correct-global`, `--write-quality`, `--write-quality-all`, `--individual-unwrapping`, `--verbose` | |
| `-template <n>`, `-threshold <x>`, `-max-seeds <n>`, `-merge-regions`, `-correct-regions`, `-wrap-addition <x>`, `-temporal-uncertain-unwrapping [x]` | same | |
| `-no-phase-rescale`, `-fix-ge-phase` | same | |
| `-B [name]` | `--compute-B0` | no name → `<base>_B0`; an explicit name is resolved relative to the main output directory, matching the Julia app |
| `-B0-phase-weighting <mode>` | `--B0-phase-weighting` | `phase_snr`\|`phase_var`\|`average`\|`TEs`\|`mag`\|`simulated_mag` |
| `-no-mask-out` | *(niimath-only)* | suppress the `<base>_mask` side output |

Side outputs normally use `nifti_save` postfixes on the already assigned output filename: `<base>_mask`, `<base>_quality`, `<base>_quality_<1..6>`, `<base>_regions`, `<base>_B0`, `<base>_B0_snr`. An explicit `-B <name>` overrides only the B0 stem as described above. They honor `FSLOUTPUTTYPE`/`-gz`, and every nonzero save status propagates. A mask is written only when a mask actually exists (default magnitude-backed `robustmask`, `qualitymask`, or an external mask); `nomask` has no invented all-ones side output.

Reject malformed syntax during parsing and data-dependent incompatibilities after dimensions are known. Match the Julia caller's echo/TE rules exactly, including echo selection before the multi-echo length check. At minimum reject an unknown mask/weight spec, invalid echo grammar, out-of-range template/echo, external weights with the wrong layout, and `-B` without echo times.

### 4.3 Loading, rescaling, and 4D

niimath has already applied `scl_slope`/`scl_inter` by the time the op runs, whereas `readphase` inspects both the scaled and the raw arrays. Preserve the original `in_hdr` and reproduce the three-way branch:

1. scaled range within 0.1 of 2π → no rescale (this is the case for the validation data, whose phase is exactly ±π);
2. else raw range within 0.1 of 2π → use the raw values (header slope/inter reset);
3. else linearly rescale to [-π, π].

Do not blindly invert the scaled float with `(scaled-inter)/slope`: that can lose stored integer values and is meaningless after an earlier chain op. For the direct/default path, obtain exact raw values before conversion or reread the original file without applying scaling; stdin needs an explicit snapshot because it cannot be reread. To keep semantics unambiguous, phase rescaling requires `-romeo` to be the first computational op. If it follows an earlier mutating op, fail with a targeted message unless `-no-phase-rescale` was supplied; further operations after `-romeo` remain supported.

`-romeo` accepts 3D or 4D (echoes on dim 4) and rejects >4D. Huge (>INT_MAX voxel) images are rejected — do **not** add `-romeo` to `kHugeSafeOps`; auxiliaries go through `nii_reject_oversize_aux`.

The magnitude must match the phase in `nx,ny,nz` and have at least `max(echoes)` volumes (ROMEO errors otherwise). An external mask must match the 3D dimensions and, as noted in §3.8, uses raw stored nonzero values. A world-transform mismatch between phase and magnitude/mask is a **warning** via `max_displacement_mm` — ROMEO does not check it, so making it fatal would be a behavioral divergence.

### 4.4 Build wiring

Adding a feature means touching all of the duplicated inventories (known-issue #2):

- `src/Makefile` — `ROMEOFLAGS=-DHAVE_ROMEO`, `ROMEOOBJS=romeo.o`, a `romeo.o` rule with `ROMEO_STRICT_FP` and its own flag-signature/dependency file, `ROMEO=0` to disable, plus every full-feature recipe (`all`, `native-noomp`, `static`, `verbose`, `debug`, `ubsan`, `sanitize`). These recipes must link the object, not also compile `romeo.c` directly. Target-specific debug/sanitizer flags must reach the object and participate in its signature; the Apple-Clang ASan recipe must compile the ROMEO object without OpenMP/libomp, just like the rest of that target.
- `src/CMakeLists.txt` — `option(ENABLE_ROMEO ... ON)` → `-DHAVE_ROMEO` + `ADDITIONAL_SRCS`, with per-source strict-FP options on GCC/Clang and an explicit `/fp:precise` policy on MSVC. Verify that release IPO/LTO preserves the source-level FP contract.
- `SuperBuild/SuperBuild.cmake` — declare and forward `ENABLE_ROMEO`.
- `src/notarize.sh` — compile a strict-FP ROMEO object separately for each architecture and link it into that slice; do not add `romeo.c` to the fast-math source line.
- `wasm`, `wasm-wasi`, and `wasm-emcc-core` — compile dedicated strict-FP `romeo_*_wasm.o` objects and link them; do not place `romeo.c` in a whole-program fast-math source list. `ROMEO=0` must remove the definition, object, and tests cleanly.

`-romeo` is off the `tiny`/`nano` builds by default. WASM includes it in the full feature build; ROMEO's strict object coexists with allineate's fast-math objects because the FP policies are per object.

## 5. Milestones

Each milestone has an objective, mechanically checkable exit criterion. Nothing is "done" on inspection.

### M0 — Oracle harness (no C yet)

Write `test/romeo_oracle.jl`, run with `JULIA_NUM_THREADS=1` under the exact environment in §2.1, that reruns each validation case and dumps every intermediate so C can be judged without modifying upstream source:

- `weights` → raw UInt8 `.bin` for byte comparison and 4D Float32 NIfTI, `permutedims` to `(nx,ny,nz,3)`, values 0–255 stored exactly
- `qmap` (voxelquality) → 3D Float32
- pre-mask intermediates: raw threshold, the post-`> threshold` binary mask, the post-first-smoothing mask, the post-`fill_holes` mask, the final mask
- `visited` (region labels) → 3D Float32
- `unwrapped` → 3D/4D Float32 (already produced as `results*/unwrapped.nii`)
- a canonical text manifest containing provenance/input hashes plus exact hexadecimal floats for `maxmag`, the four robustmask quantiles, `high_intensity`, `noise`, `threshold`, sampled raw weights, seed voxel linear index, and `new_seed_thresh`

The harness also emits small synthetic fixtures for half-even ties, non-cubic linear-index boundaries, the upstream `[30,7,30,0]` / `[30,119,0,0]` weight cases, NaN/Inf handling, the all-weights-zero failure, a stored-zero mask with nonzero scaling/intercept, degenerate magnitude values, disconnected masks, and exact ±π/±2π cases. These become the stable regression corpus; the supplied 76×76×46 data remains the realistic end-to-end corpus.

**Exit criterion:** `bash test/romeo_oracle.sh /Users/chris/src/ROMEO.jl` verifies the pinned environment, regenerates all three validation cases plus synthetic cases into `test/romeo_ref/`, and a second run produces identical hashes. Use uncompressed `.nii` and raw `.bin` for deterministic artifacts; do not require gzip container bytes to match.

### M1 — Scaffolding: build, CLI, pass-through

`romeo.c`/`romeo.h` exist, `HAVE_ROMEO` wires through every build inventory, `-romeo` parses every option in §4.2, loads and validates auxiliaries, applies the `readphase` rescale logic, and writes the (still wrapped) phase to `<out>`.

**Exit criteria:**
- `make`, `make ROMEO=0`, `make OMP=0 ubsan`, `make debug`, top-level CMake with `ENABLE_ROMEO=ON/OFF`, and `make wasm` all build clean with no new warnings.
- `make` → `make OMP=0 all` relinks with no `make -B` (flag-stamp works).
- `niimath phase0.nii.gz -romeo mag0.nii.gz -t 16.8 out` → `out.nii.gz` compares equal to the oracle's rescaled phase at threshold 0. Synthetic raw/scaled-header branches and `-fix-ge-phase` also match their oracle arrays exactly or at the stated float tolerance.
- Malformed inputs (unknown mask/weight spec, invalid TE/echo selection, 5D input, oversized auxiliary, `-dt double`, rescaling after an earlier mutating op) each fail with a specific message and nonzero status. A `ROMEO=0` binary gives the targeted feature-disabled message.

### M2 — Weights (the fidelity gate)

Port `γ`, all six weight terms, `getweight`, `updateflags`, `parsekwargs` (including the `mag .* mask` step and the Float64 `maxmag` quantile), `rescale`, the last-plane zeroing, and the `romeo`/`romeo2`/`romeo3`/`romeo4`/`romeo6`/flags selection. Add a documented test-only dump option (for example `-romeo-dump-weights <file>`) that is excluded from normal help.

**Exit criteria:**
- The numeric-type audit table (§3.1) is committed as a comment block in `romeo.c`, one row per expression and branch, each Julia type confirmed under the pinned environment. In particular, `phaselinearity` is data-dependent: ordinary interior results are Float32, while its literal fallback branches return Float64.
- Raw dumped UInt8 weights are byte-identical to the M0 oracle for all three validation cases, for every one of `romeo`, `romeo2`, `romeo3`, `romeo4`, `romeo6`, and at least four explicit flag combinations including the all-six case. The NIfTI view also passes `--compare 0`.
- The pre-`rescale` Float64 weights agree with the oracle's sampled edges to ≤ 4 ULP.
- Synthetic boundary cases establish the intended `rem2pi`/`remainder` equivalence. If libc `remainder` fails, port the required Julia reduction instead of widening the output tolerance.
- **Fast-math exposure measured and recorded:** build `romeo.o` with and without `-ffast-math -ffp-contract=fast`, report the count of differing weight bytes out of 3·N for each case. Ship the strict build regardless.

### M3 — Voxel quality map

Port `voxelquality` (the `type=Float32, rescale=identity` path plus the neighbour-sum and `/6`), and the 4D `calculateweights` overload.

**Exit criterion:** `<base>_quality` compares to the oracle `qmap` at `--compare 1e-6`, for both the single-echo and the two-echo case, and for each of the six individual-flag quality maps (`-Q`).

### M4 — `robustmask`

Port `sample`, type-7 `quantile`, the threshold estimation (0.05/0.15/0.8/0.99 quantiles, the `noise > high_intensity/10` fallback chain, `max(5·noise, high_intensity/5)`), `boxfilterline!` (with its asymmetric edge normalization), `getboxsizes`/`checkboxsizes!`, the two smoothing passes (`nbox=1, boxsizes=[5]`, threshold 0.4; then `nbox=2, boxsizes=[3,3]`, threshold 0.6), and `fill_holes` = `!imfill(!mask, (1, N/20))` with 6-connectivity.

**Exit criteria:**
- The four quantiles, `high_intensity`, `noise`, and `threshold` match the oracle scalars to ≤ 1e-6 relative.
- Each of the four mask stages (initial threshold, first smoothing, hole fill, final smoothing) compares exactly (`--compare 0`) to the oracle.
- `<base>_mask` passes `--compare 0` against `results0/mask.nii`, `results1/mask.nii`, and `results/mask.nii`.
- `-k qualitymask`, `-k qualitymask 0.25`, `-k <file>`, and `-k nomask` all behave as ROMEO does (the `qualitymask` path reuses M3 + the `threshold` argument of `robustmask`). The scaled external-mask fixture proves truth is based on raw stored values.

### M5 — Spatial unwrapping, single echo — **the headline milestone**

Port `PQueue` (LIFO buckets), `getseedqueue`, `findseed!`, `seedcorrection!` (both branches), `getseedfunction`, `grow_region_unwrap!` including `getvoxelsfromedge`/`getnewedge`/`unwrapedge!`/`unwrapvoxel`, and the 3D `unwrap!` with `correctglobal`. `maxseeds=1` only.

**Exit criteria:**
- The seed voxel index and `new_seed_thresh` match the oracle exactly.
- `visited` (region labels) compares exactly to the oracle.
- `niimath phase0.nii.gz -romeo mag0.nii.gz -t 16.8 out` → `out.nii.gz` vs `results0/unwrapped.nii` at `--compare 1e-4` **passes**, and the reported max |diff| is ≤ 1e-5 (i.e. no 2π errors anywhere, only float noise). Same for `phase1`/`results1`.
- The harness reports maximum direct difference and maximum residual after subtracting the nearest multiple of 2π. Direct difference must stay ≤1e-5; a wrap-count difference is always a failure.
- The `unwrap(phase; mag, correctglobal=true)` ≈ `phase_uw` property test from `ROMEO.jl/test/features.jl` is reimplemented in C and passes for `l ∈ 7:5:20`, offsets `2π·[0,2,-1,10,-50]`, and each of the three magnitude variants.
- Record median runtime and peak RSS for the validation case after a warm filesystem cache. Performance is reported, not a gate, until a strict-C baseline exists; fidelity work must not be distorted to meet an unmeasured 200 ms target.

### M6 — Multi-echo

Port the 4D `unwrap!`: template selection, `p2ref`, per-echo `phase2`/`TEs` pairs, the multi-echo `seedcorrection!` offset search over `off1 ∈ -2:2, off2 ∈ -1:1` with the `(|off1|+|off2|)/100` tie-break penalty, temporal unwrapping via `unwrapvoxel`, `unwrap_individual!`, and `correct_multi_echo_wraps!`. `temporal_uncertain_unwrapping` (with `unwrapped_quality`, `getseededges`, `initqueue`, and the re-entrant `grow_region_unwrap!`) lands here too even though the CLI defaults it to 0.

**Exit criteria:**
- `niimath phase.nii -romeo mag.nii -t '[16.8,38.56]' out` vs `results/unwrapped.nii` at `--compare 1e-4`, max |diff| ≤ 1e-5, on both volumes independently.
- `-i` (individual) matches a Julia `--individual-unwrapping` oracle run at the same tolerance.
- `-temporal-uncertain-unwrapping 0.5` matches a Julia oracle run at the same tolerance.
- `-template 2` matches a Julia oracle run.

### M7 — Remaining CLI behaviour

`--threshold` (zero out |phase| > n·2π), `-u` mask-unwrapped (including the "`nomask` becomes `robustmask`" promotion), `-g` correct-global for both 3D and multi-echo, `-e` echo selection, `-no-phase-rescale`, `-fix-ge-phase`, `-Q` all-quality-maps (including the "skip if all 1.0" rule), external 4D weights, `-w bestpath` (Abdul-Rahman weights: ordered neighbor offsets, `getD`, `getbestpathweight`, the `/10` scaling), and the experimental `-max-seeds`/`-merge-regions`/`-correct-regions`/`-wrap-addition` path.

**Exit criteria:** one pinned Julia oracle run per flag, each matching at the appropriate tier in §6, plus ports of the relevant upstream `features.jl`, `specialcases.jl`, `voxelquality.jl`, and `mri.jl` cases. For `-max-seeds > 1` with region merging, first determine whether the pinned Julia `Dict` iteration is stable for the fixture. If it is not a portable contract, document that limitation and gate only this experimental path on region count, merge adjacency, and wrap-consistent output rather than pretending byte equality is guaranteed.

### M8 — B0 (optional output)

Port `calculateB0_unwrapped` and `get_B0_snr` with all six `-B0-phase-weighting` modes, the `exp(-TE/20)` synthetic-magnitude fallback, and the non-finite→0 cleanup.

**Exit criterion:** `-B` and `_B0_snr` outputs vs a Julia `--compute-B0 --phase-offset-correction off` oracle at `--compare 1e-4`, for each weighting mode, both with and without a magnitude. The explicit `off` is essential because ordinary Julia `-B` enables the out-of-scope MCPC correction on multi-echo data.

### M9 — Hardening, performance, integration

- OpenMP over the weight computation only (`Threads.@threads for dim in 1:3` in Julia is over 3 dims; parallelising over voxels is also race-free). The MST loop stays serial — it is inherently sequential and it is where determinism lives.
- **Thread-count parity is a hard requirement:** identical uncompressed voxel bytes at `-p 1`, `-p 2`, `-p 4`, `-p 8` (plus `--compare 0`).
- Repository regression tests pass. A serial `make OMP=0 ubsan` build is clean on synthetic and validation cases. Then run a normal binary with `MallocScribble=1 MallocGuardEdges=1`, followed by `MallocNanoZone=0 leaks --atExit -- <command>`. Do not use Apple-Clang ASan by default on this machine.
- Every size product is checked. Buffers installed into `nim->data` use `nii_malloc`/`nii_calloc`; internal scratch uses checked plain `malloc`/`calloc` and returns a clean error. Defer side-output writes until computations succeed so an allocation failure cannot leave a misleading mask beside a missing main result.
- No `qsort`-with-comparator anywhere in `romeo.c` (WASM rule); quantiles and medians use quickselect.
- Add `-romeo` cases to `.github/scripts/release_smoke.py` — synthesized fixtures only, stdlib-only, minimum-Python-clean (use the local `_prod()` helper, no `math.prod`); `release-checks.yml`'s `vermin` guard must stay green. Include a tiny committed synthetic golden output or exact scalar/hash oracle so CI tests the result, not merely successful execution.
- WASM: `bun run build` succeeds with `-romeo` present; add a case to the js test suite.
- Docs: `README.md` op list, `niimath.c` help text, an `AGENTS.md` section covering the module, the strict-FP requirement, and the fidelity hazards in §3.

## 6. Validation harness

`test/romeo_compare.sh <niimath-binary>` drives the whole thing:

```bash
cd /Users/chris/src/ROMEO.jl/romeo
julia ../romeo.jl phase0.nii.gz -m mag0.nii.gz -t '[16.8]'       -o results0
julia ../romeo.jl phase1.nii.gz -m mag1.nii.gz -t '[38.56]'      -o results1
julia ../romeo.jl phase.nii     -m mag.nii     -t '[16.8,38.56]' -o results
# ...then, per case:
niimath phase0.nii.gz -romeo mag0.nii.gz -t 16.8 nm0
niimath nm0.nii.gz      --compare 1e-4 results0/unwrapped.nii
niimath nm0_mask.nii.gz --compare 0    results0/mask.nii
```

Three tolerance tiers, used deliberately:

- **raw byte comparison plus `--compare 0`** — integer-valued intermediates: UInt8 weights, masks, region labels. `--compare 0` is value-oriented rather than a file-byte comparator, so raw dumps/hashes are the authoritative byte-exact check for weights.
- **`--compare 1e-6`** — the quality map and other pure float32 reductions.
- **`--compare 1e-4`** — the unwrapped phase. Wide enough to absorb float noise, far narrower than the 6.28 that a single wrap error would produce. Always accompanied by the reported max |diff| so a "pass" at 9e-5 is visible rather than silent.

The current `--compare` implementation explicitly counts non-finite mismatches and fails them, including NaN-vs-finite and opposite infinities. Keep an independent finite/NaN/Inf count in the ROMEO harness anyway, because it makes the report diagnostic and guards future compare regressions.

The supplied validation data is small (76×76×46, 2 echoes, ~1 MB per Float32 volume) and lives outside this repo. Do not vendor it; the harness takes the ROMEO.jl checkout path as an argument and skips with a clear message if absent. Commit only synthetic fixtures and small derived golden artifacts with recorded provenance.

## 7. Out of scope

Explicitly not ported, and rejected with a clear message rather than silently ignored:

- **MCPC-3D-S phase-offset correction** (`--phase-offset-correction`, `--write-phase-offsets`, `--phase-offset-smoothing-sigma-mm`) and multi-channel coil combination — these live in `MriResearchTools.mcpc3ds` and are a substantial separate port with their own complex-smoothing machinery. ROMEO's ordinary multi-echo `-B` activates monopolar correction; niimath must state that it computes B0 without MCPC and require an explicit acknowledgement/option for multi-echo B0 rather than silently implying full CLI equivalence.
- 5D (multi-channel) input.
- Memory mapping (`--no-mmap`) — an I/O strategy, not an algorithm.
- The `citations_romeo.txt` / `settings_romeo.txt` sidecar files. Put the ROMEO citation (Dymerska et al. 2020, doi:10.1002/mrm.28563) in help/README and print it in verbose mode; do not add unsolicited stderr output to every successful chain.

## 8. Decisions recorded by this plan

1. **Enabled by default, still optional:** native/CMake/full WASM builds define `HAVE_ROMEO` by default; `ROMEO=0` / `ENABLE_ROMEO=OFF` disables it. `tiny`/`nano` omit it.
2. **Filename-based side outputs:** use `<out>_mask` and the other postfixes in §4.2. ROMEO's fixed names inside an output directory do not fit niimath's chain model.
3. **Chain op only for the first implementation:** expose `-romeo`, not a second `--romeo` subcommand. Later niimath operations may follow it; phase rescaling may not follow an earlier mutating operation unless disabled, as specified in §4.3.
4. **MIT compliance is a deliverable, not an open wording question:** add an attribution header to `romeo.c` naming the pinned ROMEO.jl and MriResearchTools sources, preserve both upstream copyright/permission notices in a shipped `src/romeo.LICENSE` (or repository-wide third-party notice), and include that file in source/binary/npm release inventories where applicable. Algorithm citations remain separately documented in help/README.
