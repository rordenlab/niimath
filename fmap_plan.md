# fmap_plan.md — B0 fieldmap preparation and EPI unwarping in niimath

> **HISTORICAL — superseded by `fmap_bench/test/fmap_reference_manifest.md`.**
>
> This is the plan the work was executed against, kept because the manifest cites it by name and because two of its predictions were overturned *by measurement*, which is worth being able to read. It is **not** a contract and must not be updated: a second document that can drift from the manifest is exactly the failure AGENTS.md already records for `skullstrip_plan.md`.
>
> Two decisions here were deliberately reversed during implementation:
> * **Decision 4** said the apply stage would reuse `medic_unwarp`'s `md_pull`. It does not — `md_pull` is a 3D Lanczos-5 pull over a millimetre map writing to a second buffer, while `-fugue` is a 1D linear pull over a voxel shift that rewrites in place. See the note at the foot of `src/fmap.c`.
> * **M4** said `ROMEO=0` must imply `FMAP=0`. It does not — only `-fmapprep` needs ROMEO, so a `ROMEO=0` build keeps a working `-fugue`, and CMake emits a `message(STATUS)` rather than erroring.
>
> Its M3 prediction that the regularisation chain would be the schedule risk was also wrong in an interesting way: there is no regularisation chain at all.

**Goal.** Give niimath a BSD-2 replacement for the two FSL tools that carry the fieldmap half of a FUGUE-style susceptibility-distortion pipeline: `fsl_prepare_fieldmap` (Siemens phase-difference → rad/s fieldmap) and `fugue` (fieldmap + dwell + phase-encoding direction → unwarped EPI). `bet` is explicitly out of scope. Two new chain ops, `-fmapprep` and `-fugue`, validated end to end against the FSL binaries on `/Users/chris/src/fmap_bench/generic`, with the evidence published in the `fmap_bench` repository the way `moco_bench` and `medic_bench` already do for `-moco`, `-stc` and `--medic`.

**Why this is worth doing.** On the benchmark dataset FSL spends 12.7 s and 650 MB peak, of which `fugue` alone is 11.3 s and all 650 MB. The apply stage is a separable 1D resample along one axis over a 76×76×45×254 series (66 M voxels, 264 MB as float32) — memory-bandwidth-bound work that niimath should do several times faster in roughly half the RAM. `fsl_prepare_fieldmap` (819 ms) is small, but folding it in removes the last non-niimath step apart from brain extraction.

---

## Strategic decisions (locked — do not relitigate mid-run)

1. **Strict clean-room, with a hard guard.** FSL's source for these tools *is* present on this machine at `/Users/chris/fsl/src/fsl-fugue/` (`fugue.cc`, `prelude.cc`, `unwarpfns.cc/.h`, `fsl_prepare_fieldmap.tcl`). It is licensed under `/Users/chris/fsl/LICENCE.FSL` — Oxford, non-commercial, "no part of the Software may be reproduced, modified, transmitted or transferred" — which is incompatible with niimath's BSD-2. **Nothing under `/Users/chris/fsl/src/` may be read, grepped, opened, summarised or quoted at any point in this plan, by any agent.** This is the same posture AGENTS.md records for the MCW AFNI files behind `-moco` and `-stc`, and unlike those it has *not* been lifted by any relicense.
2. **The FSL binaries are black-box oracles.** Running them and measuring their inputs and outputs is the intended method and is what every FSL user does. The line: *allowed* — normal invocation, whatever the tool prints to stdout/stderr in ordinary use, the files it writes, and anything derivable by feeding it synthetic inputs; *forbidden* — reading shipped source, and equally `set -x`, `strace`/`dtruss`, disassembly, or any other instrumentation whose purpose is to recover the internal command sequence. A runtime trace of `fsl_prepare_fieldmap` would be a transliteration of the `.tcl` we have agreed not to open, so it buys the same information by another door. If a convention cannot be pinned down without it, record that in the manifest as an open question and settle it by fitting synthetic inputs instead.
3. **BET stays FSL's, for now.** `fsl_prepare_fieldmap` needs a brain-extracted magnitude, and both pipelines will use the *same* `bet` output. This is deliberate experimental control: holding the mask identical means every measured difference is attributable to our two ops and not to a different brain mask. Replacing `bet` is a follow-on, out of scope here; the ops must therefore accept a caller-supplied mask/brain-extracted magnitude as an argument rather than deriving one internally.
4. **The apply stage reuses the existing resampler.** niimath already has `-unwarp <map> <axis>` (`medic_unwarp`, medic.c), which resamples a 3D/4D float32 image through a scalar displacement map in **millimetres** along one axis using a Lanczos-3 pull. `-fugue` converts rad/s + dwell + unwarpdir into exactly such a map and calls that kernel, so the tree keeps one interpolation implementation. If a measured FUGUE convention turns out to be incompatible with `md_pull` (see M1), extend `md_pull` under a flag rather than forking it — and say so in the manifest.
5. **The gate is the final unwarped EPI.** ROMEO and PRELUDE are different unwrapping algorithms and will legitimately disagree in wrapped, noisy and near-mask-edge voxels, so gating the intermediate rad/s fieldmap tightly would fail for the wrong reason. Primary acceptance is on the end product, via `niimath --compare`, with the fieldmap itself checked more loosely as a secondary signal. Evidence lives in `/Users/chris/src/fmap_bench`.

### Naming (low-cost, override freely)

`-fmapprep` and `-fugue`. Naming an op after the tool it emulates is established here — `-allineate` and `-qwarp` are AFNI tool names, `-deface` and `-moco` are descriptive. Both are ordinary chain ops (single `-`), not terminal `--` subcommands, so they compose with the rest of the calculator.

```
niimath phasediff.nii.gz -fmapprep <mag_brain> <deltaTE_ms> fmap_rads.nii.gz
niimath bold.nii.gz      -fugue <fmap_rads> <dwell_s> <unwarpdir> bold_unwarped.nii.gz
```

`-fmapprep`'s argument order mirrors `fsl_prepare_fieldmap SIEMENS <phase> <mag_brain> <out> <deltaTE_ms>` minus the pieces niimath gets from the chain. Both must use the `ac + 1 <= argc` peek documented under the op-loop gotcha in AGENTS.md if either gains an optional trailing sub-option, so a trailing flag can never become the output filename.

---

## Immediate housekeeping (do first, before any code)

- **`fugue.pdf` must never be committed.** It is currently untracked in the niimath repo root and `.gitignore` does not cover it. It is Jenkinson, *Magn Reson Med* 49:193–197 (2003), © Wiley-Liss — a copyrighted article. Move it out of the working tree (`/Users/chris/src/fmap_bench/papers/` or anywhere outside a git repo) **or** add it to `.gitignore`; moving it out is safer, because `.gitignore` does not protect against `git add -f` or a future `git add -A` in a fresh clone. Verify with `git status --porcelain` that it no longer appears.
- **Note what that PDF actually is.** It is the **PRELUDE** paper — the region-merging N-dimensional phase unwrapper — not a description of FUGUE. It is the published basis for the *unwrapping* stage, which we are replacing with ROMEO rather than reimplementing, so it is background reading, not a specification. FUGUE's own method traces to Jezzard & Balaban, *Magn Reson Med* 34:65–73 (1995), plus the public FSL documentation at <https://fsl.fmrib.ox.ac.uk/fsl/docs/registration/fugue.html>.
- **Add the guard mechanically.** Before any milestone that spawns subagents, confirm the prohibition in decision 1 is restated verbatim in every subagent prompt that could plausibly go looking for an implementation.

---

## Background: what the two stages actually do

Stated from the public documentation and standard EPI physics only.

**Prepare.** A Siemens `gre_field_mapping` acquisition yields a magnitude image and a *phase difference* between two echoes separated by ΔTE. The phase difference is wrapped into a 2π range and stored as scaled integers. Recovering a field offset means: rescale the stored phase to radians, unwrap it (PRELUDE in FSL, ROMEO for us) using the magnitude as an anatomical guide and a brain mask to bound the problem, then divide by ΔTE in seconds to get rad/s. The result is regularised — the documentation names median filtering, despiking and smoothing as the available tools — and its median inside the mask is removed, because a constant field offset is unobservable and would otherwise translate the entire EPI. Some form of edge handling extrapolates the fieldmap outward so that unwarping near the mask boundary does not read zeros.

**Apply.** In EPI the phase-encoding direction is traversed slowly, so an off-resonance frequency `Δf` accumulates phase across the readout and displaces signal along that one axis. The displacement in voxels is `Δf [Hz] × dwell [s] × N_pe`, equivalently `Δf × TotalReadoutTime` — with the fieldmap in rad/s that is `fmap / (2π) × dwell × N_pe`. The benchmark's own driver already computes its expected maximum shift this way (`bids_fmap.py:530`), which is a useful independent check on whatever M1 measures. Unwarping resamples the distorted image back onto the true grid; forward-warping (`-w`) does the opposite and is used to push an undistorted target into EPI space. Intensity (Jacobian) correction, which compensates for the signal pile-up that accompanies compression, is a separate documented option (`--icorr`) and is **not** on by default — M1 must confirm this rather than assume it.

Two consequences worth holding on to. First, both unwrappers are only determined up to a global 2π multiple per connected region, and the demedian step is what makes that harmless — which is precisely why an end-to-end gate can pass while the raw unwrapped phase differs. Second, a constant error in the fieldmap becomes a constant shift in the EPI, so a *systematic* offset is exactly the failure a bounded max-diff gate catches and a bare correlation floor does not.

---

## Benchmark data and reference numbers

Dataset: `/Users/chris/src/fmap_bench/generic`, subject `sub-fm`.

| Image | Shape | Stored type | Key metadata |
| --- | --- | --- | --- |
| `fmap/..._phasediff.nii.gz` | 76×76×45 | int16 | EchoTime1 4.00 ms, EchoTime2 6.46 ms → **ΔTE 2.46 ms** |
| `fmap/..._magnitude1.nii.gz` | 76×76×45 | — | input to `bet` |
| `func/..._part-mag_sbref.nii.gz` | 76×76×45 | uint16 | PE `j`, dwell 0.00054001 s, TRT 0.0405008 s |
| `func/..._part-mag_bold.nii.gz` | 76×76×45×**254** | uint16 | same; 66 M voxels, 264 MB float32 |

The `sbref` is a single volume on the same grid as the `bold` and is the right first target: it makes the apply stage a 3D problem, so a convention error shows up in one image instead of 254.

FSL baseline from `fmap_bench/README.md` (`fugue` appears twice — the 4D `bold` and the 3D `sbref`):

| Tool | Time ms | Peak RAM MB |
| --- | --- | --- |
| bet | 474 | 18 |
| fsl_prepare_fieldmap | 819 | 20 |
| fugue (bold, 254 vol) | 11343 | 650 |
| fugue (sbref, 1 vol) | 110 | 19 |
| **TOTAL** | **12745** | **650** |

Note that FSL's 650 MB is ~2.5 copies of the float32 series. `medic_unwarp` as written allocates a full second copy of `nim->data` plus the map, so a naive `-fugue` lands near 530 MB — better than FSL but not by much. Beating it decisively is a per-frame streaming change, which is why performance is its own milestone (M6) rather than an afterthought.

Also note: `/Users/chris/src/fmap_bench/bids_fmap.py` and `bids_niimath_fmap.py` are currently byte-identical copies. The first stays the FSL reference pipeline and should not drift; the second becomes the niimath pipeline in M5.

---

## Milestones

Each milestone has an explicit exit criterion. Do not start the next one until the current one's criterion is demonstrated, and record the demonstration in the manifest. Where a criterion says "measured", it means a number captured from a run, not an expectation.

### M0 — Harness, oracle capture, and the manifest skeleton

Stand up `/Users/chris/src/fmap_bench` as a real repository (it currently has no git history at all — see the standing complaint in AGENTS.md that `strip_bench` and `skullstrip_bench` are unpublishable, and do not repeat it here). Deliverables:

- `git init` plus a first commit; the ~300 MB of derivatives stay untracked or under LFS-free ignore rules, but the scripts, manifest and measured tables are committed.
- `test/fmap_reference_manifest.md`, following the shape of `moco_bench/test/moco_reference_manifest.md`: one section per measured convention, each with the experiment that established it. Open with the provenance statement (decision 1 above) and the exact FSL version string from `cat $FSLDIR/etc/fslversion`.
- A capture script that regenerates the FSL reference tree from scratch and records wall time and peak RSS per tool, so the baseline table above is reproducible rather than quoted.
- A comparison helper wrapping `niimath --compare`, which is the project's established oracle. **Heed the `--compare` gotcha in AGENTS.md**: its verdict is magnitude-based, so a non-finite difference has no magnitude and an all-NaN-versus-finite pair reads as *equal* and exits 0. Count non-finite mismatches separately, and assert the count is zero.
- Confirm the `-p 8` thread banner appears and that timing runs pass `-gz 0` on both sides, per the benchmark traps AGENTS.md records for `moco_bench`.

**Exit:** the FSL reference tree regenerates from one command, the baseline timing/RAM table is reproduced within noise, and `--compare` of a file against itself and against a deliberately NaN-poisoned copy both give the expected verdicts.

### M1 — Measure the apply-stage conventions

Hold the fieldmap constant by feeding FSL's *own* `fmap_rads.nii.gz` into both sides. That isolates the apply stage completely — any difference is interpolation, sign or scaling, never unwrapping.

Use synthetic fieldmaps so each convention is separable rather than tangled: a constant field (isolates the sign and the shift constant — a constant rad/s must produce a pure whole-image translation whose size you can read off directly), a linear ramp along each axis in turn (isolates the axis mapping and tests compression/expansion), a single-voxel delta (exposes the interpolation kernel's footprint and hence its identity), and a field large enough to push samples out of the FOV (exposes the edge/fill rule). Run each through `fugue` on the `sbref`.

Conventions to pin down and record, each with its experiment:

- The exact shift constant. Is it `fmap/(2π) × dwell × N_pe`, or `× (N_pe − 1)`, or something referenced to the image rather than the acquisition matrix? A constant-field test answers this in one number. Cross-check against `bids_fmap.py:530`.
- Sign per `--unwarpdir` value, for all six of `x y z x- y- z-`, and how that maps onto NIfTI axis order and the sform/qform determinant. The benchmark data has PE `j` for the EPI and `j-` for the fieldmap, so the two are not interchangeable and a sign error is very easy to make here.
- Pull versus push: is `-u` output(v) = input(v + shift) or input(v − shift)? `medic_unwarp` is a pull; confirm the direction matches or invert at the map-construction step.
- Interpolation kernel and whether it is 1D along PE only or 3D. `md_pull` is a 3D Lanczos-3 separable pull. If FUGUE is 1D-along-PE, a 3D kernel will differ measurably in the other two axes even where the shift is zero — this is the single most likely source of a systematic mismatch, so measure it early with the delta-field test.
- Out-of-FOV behaviour: zero fill, edge clamp, or wrap.
- Whether intensity/Jacobian correction is applied by default (documentation says `--icorr` is opt-in; verify with a ramp field, where a Jacobian term produces an intensity gradient a pure resample does not).
- What `fugue` does with a 4D input — per-volume identical treatment, presumably, but confirm, and confirm the fieldmap is broadcast rather than expected per-frame.
- The output datatype and any scaling applied on write.

**Exit:** a table in the manifest giving each convention, its measured value, and the synthetic experiment that established it — sufficient that someone could implement the op from the table alone without touching FSL.

### M2 — Implement `-fugue` (apply)

New `src/fmap.c` / `fmap.h`, `HAVE_FMAP`, DT32 only, ordinary FP, ordinary chain op. Deliberately **not** in `kHugeSafeOps` (see the huge-image section of AGENTS.md; `-unwarp` is already excluded and `-fugue` inherits the same INT_MAX limit through it). Diagnostics through `printfx`/stderr only, never bare `printf` — niimath writes the image to stdout when the output name is `-`, and the `-skullstrip` bullet in AGENTS.md records exactly what a stray `printf` does to a piped NIfTI.

Structure: validate that the fieldmap shares the input's grid (reuse the `max_displacement_mm()`-based same-grid test `md_same_grid` already applies, at the same 0.001 mm tolerance); build a millimetre displacement map from rad/s, dwell, unwarpdir and the measured constant; hand it to the `medic_unwarp` resampler. Fail closed on a singular world transform, a mismatched grid, a non-finite fieldmap voxel, an unrecognised `unwarpdir`, and a missing argument (`NII_NEED_ARGS`, per the twelve-op sweep recorded in AGENTS.md).

**Exit:** with FSL's own fieldmap as input, `-fugue` on the `sbref` matches `fugue -u` to a correlation ≥ 0.999 inside the brain mask and a bounded max relative difference, with **zero** non-finite mismatches; the same on the 4D `bold`; and byte-identical output across `-p 1 / -p 2 / -p 8`. Thread-count reproducibility is not optional — every op in this tree except `-qwarp` holds it, and `-fugue` has no reduction that would make it hard.

### M3 — Measure the prepare-stage conventions

Same discipline, holding `bet`'s mask constant. What to pin down:

- The phase rescaling actually expected on input. `bids_fmap.py` already knows three encodings (`0..4095`, `-4096..4094`, `-π..π`) and rescales onto the range `fsl_prepare_fieldmap SIEMENS` wants; confirm what the tool does with each by feeding it all three forms of the same field and comparing outputs.
- Where ΔTE enters and in what units (the CLI takes milliseconds; the output is rad/s).
- What the mask does: is it a hard restriction on the unwrap, on the output, or both? Does the output mask equal the input mask, or is it eroded, or dilated?
- The regularisation chain and its parameters. Documentation names median filtering, despiking and smoothing as available; the shipped defaults must be *fitted* from oracle behaviour on synthetic fields with known spikes and known noise, not guessed and not recovered from a trace (decision 2). Budget real time for this — it is the least well specified part of the job and the most likely to need several rounds.
- The demedian step: median of what, over which voxels (mask interior only?), subtracted where.
- Edge extrapolation outside the mask: measure by comparing a fieldmap prepared with a small mask against one prepared with a large mask over the same field.

**Exit:** manifest table as in M1. Where a convention cannot be established to better than a stated tolerance without instrumentation, say so explicitly and record the residual — an honest open question beats a confident guess, and AGENTS.md's MEDIC entry sets that precedent (`--rank 0`, the one residual that could not be reproduced).

### M4 — Implement `-fmapprep`

Uses `romeo_unwrap_frame()` (romeo.h) in memory — do not shell out and do not duplicate ROMEO, exactly as `--medic` does. Single echo, `neco = 1`, `TEs = [ΔTE]`, `mask_in` = the caller-supplied brain mask so ROMEO's own `robustmask` is bypassed and the mask stays identical to FSL's. Rescale stored phase to radians first; support at least the three encodings `bids_fmap.py` knows, and **detect rather than assume** — a wrong assumption here scales the entire fieldmap.

`-fmapprep` requires `HAVE_ROMEO`, so `ROMEO=0` must imply `FMAP=0` the way it already implies `MEDIC=0` (Makefile silently, CMake with a hard error). Do not let a `ROMEO=0` build advertise the op.

**Exit (two-tier, per decision 5).** *Primary:* the full niimath chain — `-fmapprep` then `-fugue` — matches FSL's final unwarped `sbref` and `bold` to a correlation ≥ 0.99 inside the brain mask, with a max difference bounded at a threshold justified in the manifest and zero non-finite mismatches. *Secondary:* the rad/s fieldmap itself agrees with FSL's inside the mask to a looser, explicitly stated tolerance, with the disagreement characterised — if it is concentrated at mask edges and in low-magnitude voxels that is the expected ROMEO-versus-PRELUDE signature; if it is a global scale or offset, that is a bug in the ΔTE or demedian step and the primary gate must not be used to paper over it. Report both numbers; gate on the primary.

### M5 — The niimath pipeline

Rewrite `/Users/chris/src/fmap_bench/bids_niimath_fmap.py` to call `bet` (FSL, held constant) plus `niimath -fmapprep` and `niimath -fugue`, keeping the same BIDS discovery, sidecar handling, unsupported-suffix reporting and derivative layout as `bids_fmap.py` so the two are directly comparable. Keep `bids_fmap.py` frozen as the reference.

**Exit:** both pipelines run to completion on `generic/`; a comparison report over every published derivative meets the M4 gate; the README carries a side-by-side time and peak-RAM table for FSL versus niimath, generated by the M0 capture script rather than typed in.

### M6 — Performance

Only now, with correctness pinned and reproducible. The apply stage is the target: 11.3 s and 650 MB for 66 M voxels is bandwidth-bound work with headroom.

- **RAM.** The current `medic_unwarp` allocates a full `nvox` output buffer and `memcpy`s back. For 4D input, process frame by frame into a single-frame scratch buffer and write back in place — peak drops from ~2 copies to ~1 copy plus one frame (264 MB → ~270 MB total including the map, against FSL's 650 MB). This changes `medic_unwarp`, which `--medic` also uses, so `--medic`'s own outputs must be verified unchanged.
- **Time.** Frame-parallel OpenMP is already there (`frame_parallel = nt >= omp_get_max_threads()`); 254 frames will take it. Check whether the Lanczos-3 3D kernel is doing more work than the measured FUGUE convention requires — if M1 finds the shift is 1D along PE, a separable 1D pass is dramatically cheaper than a 3D one and is not an approximation.
- Any reduction introduced here must be thread-count invariant. `max` is exact; `+` on doubles is not — use the fixed-chunk pattern from `coreg_fast.c` (`CF_CR_NCHUNK`) if a sum is unavoidable.

**Exit:** measured wall time and peak RSS for both ops on `bold` and `sbref`, published in the README table; output byte-identical to M4's across `-p 1/2/8` and unchanged from the pre-optimisation build; `--medic`'s regression coverage still green.

### M7 — Integration

The feature-list drift problem is documented in AGENTS.md as known issue 2 and it is a real recurring cost. `-fmapprep`/`-fugue` are structurally identical to MOCO and STC — ordinary FP, no strict-FP object, no separate flag stamp — so follow those two exactly:

- `src/Makefile` (`FMAPFLAGS` across `all` / `native-noomp` / `static` / `debug` / `verbose` / `ubsan` / `sanitize`), plus `FMAP_WASM` for the Emscripten target, kept separate rather than folded into another target's flags so a `ROMEO=0`/`MEDIC=0` build cannot silently drop it.
- `src/CMakeLists.txt` and `SuperBuild/SuperBuild.cmake` (`ENABLE_FMAP`, forwarded).
- `src/notarize.sh` — and **never** add GPL or FSL-derived sources there; that script builds the BSD-2-branded shipped macOS binary.
- Deliberately absent from `WASI_SRCS`, `wasm-emcc-core`, `tiny` and `nano`, matching MEDIC/MOCO/STC.
- `#ifdef`-paired help lines in `niimath.c`, so a `FMAP=0` build still tells the user how to re-enable it.
- `release_smoke.py` coverage: dispatch, rejection paths (4D fieldmap, mismatched grid, bad `unwarpdir`, missing argument, `-dt double`), closed-form numeric checks on a synthetic linear field where the correct answer is known analytically, and a `-p 1` versus `-p 8` byte-equality check. Keep it stdlib-only and minimum-Python clean — it runs inside every wheel Python from `requires-python` up, and `release-checks.yml` guards this with vermin.
- An AGENTS.md entry under Source Architecture recording the licence boundary (decision 1), the measured conventions' home (`fmap_bench`, quoted not restated), the ROMEO dependency, and the `medic_unwarp` sharing.

**Exit:** `make test` green; a `ROMEO=0` build, an `OMP=0` build and a `MOCO=0 STC=0` build all configure, compile and link; CMake and Makefile binaries agree on the benchmark to the M4 gate.

---

## Risks, ranked

1. **The prepare-stage regularisation chain (M3) is the schedule risk.** It is the least publicly specified part and the part we have most deliberately denied ourselves a shortcut to. Mitigation: it sits behind the end-to-end gate, and a fieldmap that is under- or over-regularised relative to FSL's still produces a usable correction — so if it stalls, ship `-fugue` (M2) alone and let `-fmapprep`'s tolerance be looser and stated. `-fugue` is where all the time and memory is anyway.
2. **Interpolation mismatch (M1).** If FUGUE turns out to do something `md_pull` cannot express, decision 4 says extend rather than fork. Watch for this becoming a slow rewrite of the resampler under a different name.
3. **ROMEO versus PRELUDE divergence exceeding the primary gate.** Plausible in low-SNR data. If it happens, characterise *where* before touching thresholds — mask-edge disagreement is expected and can be excluded by eroding the comparison mask; a global offset is a bug. Do not loosen the gate to make a bug pass.
4. **Single-dataset validation.** One subject, one manufacturer, one PE direction. The measured sign conventions for the other five `unwarpdir` values will be established synthetically (M1), which is the right way, but real data with `i`/`k` phase encoding would be worth acquiring before this is described as general. State this limitation in the manifest rather than letting it be inferred.
5. **Accidental provenance contamination.** The forbidden source sits in a plausible search path on the same machine. Restate the prohibition in every subagent prompt; treat any transcript showing a read under `/Users/chris/fsl/src/` as grounds to discard the work product downstream of it.

## Non-goals

Brain extraction (`bet`); `prelude` as such (ROMEO substitutes for it); pepolar/`topup`-style fieldmaps and the `_epi` BIDS suffix; two-separate-phase-image (`_phase1`/`_phase2`) fieldmaps; registering a fieldmap into EPI space when the grids differ; forward warping (`fugue -w`) and `--phaseconj`; k-space unwarping (`--nokspace` is the image-space method, and the k-space path is out of scope); GE/Philips fieldmap conventions.
