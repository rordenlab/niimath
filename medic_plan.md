# Plan: clean-room MEDIC emulation using niimath and ROMEO

## Status (2026-07-25)

**M0-M12 implemented and merged** on the `romeo` branch: `src/medic.c`/`medic.h` behind `HAVE_MEDIC` (requires `HAVE_ROMEO`), `--medic` and `-unwarp`, the stdlib-only `medic.py` wrapper, build wiring across Makefile/CMake/SuperBuild/notarize.sh, and an analytic check in `release_smoke.py`. The black-box measurement record is `test/medic_reference_manifest.md`; the MCPC-3D-S patent analysis is `prior_art.md`.

**M3 (`-unwarp`) passes its gate outright**: fed the reference's own displacement map it reproduces the reference's corrected magnitude at nrmse 3.5e-5 / 4.7e-5, corr 0.999999999, zero non-finite.

**`--medic` does not yet match the reference end-to-end, and no equivalence is claimed** (bit identity is a stated non-goal, §11). Given a shared mask the native field map matches to p99 = 0.0027 Hz — so the regression, MCPC-3D-S, unwrapping, rescaling and echo handling are all correct — but ~0.24 % of voxels land on a different 2*pi branch and the inversion smears that along the phase-encoding line. With the built-in `robustmask` default the divergence is larger (p95 ~46 Hz) because the two tools' masks differ. Per the §7.2 directive we are NOT chasing the reference's mask; `--mask` is the answer.

**M5's displacement gate is currently missed for `j`.** End to end on the sbref demo with the reference mask supplied, displacement p95 is 0.059 mm against a 0.05 mm gate (`j-` is at 0.029 mm and passes). The isolated inversion/scaling convention checks do pass the same threshold (manifest §3.4, §3.5b), but those feed the reference's own native field through one formula and are not pipeline results — manifest §5.4 attributes each figure. Stated here so the gate table below is not read as fully green.

Open, deliberately not guessed: a broadband residual survives the reference's rank-10 truncation on real 170-frame data while ours is strictly rank 10 (§7.5). `--rank 0` disables the filter, and a non-finite field-map series is now a hard error — checked on the field series itself immediately after the regression, so it applies whether or not the low-rank filter runs — rather than a silently all-zero output set.

**One audit round has since landed** (`audit_response.md`): the phase-encoding polarity is now honoured by `--medic` (manifest §3.5b — the highest-value finding, since the supplied three-echo data is `j-`), `--weights` governs both unwrapping stages, the temporal correction implements the paper's Eq. 6 cumulative fit for three or more echoes from an order-independent snapshot, and the memory model, grid matching, unit handling and output atomicity were corrected as recorded below. Residual items are ranked at the end of `audit_response.md`; the top one is test coverage for those fixes.

### Decided: memory model is in-RAM, and that is not a defect

Streaming (the old §5.2) is a **non-goal**. A 4D `.nii.gz` cannot be seeked, so essentially every tool — including this one and the reference — reads and writes whole volumes in RAM anyway; a streaming layer would buy nothing for the dominant gzip case while adding a large validated surface. The requirement is instead to be **fast and honest about the RAM cost**: document the working-set formula, print it at startup, and state the wasm ceiling.

- Work arrays: `phase + mag + fields + fu + disp` = `n3 * T * (2*echoes + 3) * 4` bytes (1.18 GiB on the 170-frame two-echo run), printed to stderr at startup. That banner is a **budget, not a peak**, and the load ordering is what keeps the two close: inputs are validated from their headers alone, so no payload is resident when the arrays are allocated, and the repack then loads and frees **one echo pair at a time** — a transient overshoot of one echo pair regardless of echo count. (An earlier revision read all `2*echoes` payloads during validation and freed them only after the allocation, making the real peak `(4*echoes + 3)` series, 1.85 GiB there, against a `(2*echoes + 3)` banner. Do not reinstate that ordering.) Phase is unwrapped **in place**; the separate unwrapped-phase series was deleted during the audit, taking the gzipped peak from 2.53 GB to **2.35 GB**. Measured peak RSS on the 170-frame run: **2.04 GB** single-threaded and **2.24 GB** at 8 threads writing uncompressed, **2.35 GB** gzipped, against the reference's **3.40 GB** for the same work (manifest §5.3, authoritative for every MEDIC timing and memory figure) — all measured **before** the per-echo-loading change, so they now bound the current binary rather than describing it.
- **wasm32 has a 4 GiB address space** and every wasm target sets `-DFORCE_INT32_MAX`. `--medic` ships in the default Emscripten build but long multi-echo runs will not fit — roughly 10 GiB for 5 echoes × 600 frames at this resolution: estimate natively, apply `-unwarp` in the browser. MEDIC is omitted from `tiny`/`nano` and from the WASI reactor.

### Build requirement for any timing claim

Build against **zlib-ng** (or zlib-cloudflare), not system zlib — see M12.

niimath's ROMEO port is complete and is the foundation for this work. It already provides strict-FP, multi-echo ROMEO unwrapping, robust masks, quality maps, and B0 estimation with all six ROMEO weighting modes. The current parity suite reports 602/602 checks. MEDIC must reuse these kernels; it must not duplicate ROMEO or invoke the niimath CLI once per frame.

## 1. Goal and scope

Implement the useful core of Warpkit's MEDIC workflow without reading or copying its institutionally licensed implementation:

1. Estimate one B0 field map for every frame of a multi-echo EPI run.
2. Convert each field map to a phase-encoding displacement map.
3. Invert the displacement map onto the undistorted grid.
4. Apply the frame-matched displacement map to each magnitude series.
5. Provide a small, dependency-free Python front end that discovers BIDS files and metadata.

The scientific specification is Van et al., *Imaging Neuroscience* 4 (2026), doi:10.1162/IMAG.a.1262 (`warpkit.pdf`). The phase unwrapping foundation is the completed MIT-licensed ROMEO/MriResearchTools port in `src/romeo.c`.

The first release targets:

- coil-combined magnitude and phase images;
- at least two echoes;
- 3D single-frame and 4D BOLD inputs;
- one scalar displacement per voxel along a BIDS phase-encoding axis;
- magnitude correction only;
- float32 output;
- one-file NIfTI input (`.nii`, `.nii.gz`, and `.nii.zst` where supported).

It does not target every Warpkit utility. The goal is MEDIC-compatible behavior, not a clone of Warpkit's package structure.

Proposed user-facing commands:

```bash
python3 medic.py /path/to/sub-fm/ses-1/func --out-dir derivatives/medic

niimath --medic \
  --magnitude echo-1_mag.nii.gz echo-2_mag.nii.gz \
  --phase echo-1_phase.nii.gz echo-2_phase.nii.gz \
  --te-ms 16.8,38.56 \
  --total-readout-time 0.02025 \
  --phase-encoding-direction j \
  --out-prefix out/sub-fm

niimath echo-1_mag.nii.gz \
  -unwarp out/sub-fm_displacementmaps.nii j \
  echo-1_mag_undistorted.nii.gz
```

`--medic` is a terminal subcommand, like `--dtifit` and `--qc`. `-unwarp` is an ordinary chain operation.

## 2. Clean-room boundary

Warpkit's implementation, tests, build products, source archives, and debug symbols are out of bounds. Do not inspect them even if locally accessible.

Allowed sources:

- `warpkit.pdf`, under CC-BY 4.0;
- the MIT-licensed ROMEO.jl and MriResearchTools.jl sources already used for `romeo.c`;
- user-facing `wk-* --help` output;
- the user-provided demo scripts and README;
- supplied input data;
- output images, headers, logs, timings, and exit behavior produced by invoking the installed Warpkit commands as a black box;
- independently designed synthetic experiments.

Do not derive behavior from Warpkit unit tests. If an implementation convention is not specified by the paper or an allowed upstream project, measure it through the public executable and record the experiment.

Create `test/medic_reference_manifest.md` before implementation. For every reference run, record:

- Warpkit version and executable path;
- exact command line and environment affecting results;
- input and output SHA-256 values;
- output dimensions, datatype, transforms, units, scaling, intent, and finite-value range;
- wall time and peak RSS;
- whether the result is specified by the paper or only observed from the black box.

Reference NIfTIs are large and should remain gitignored. Commit only manifests, compact numeric summaries, and redistributable synthetic fixtures.

## 3. Paper-defined algorithm

For each frame, with echo times `t_e`, wrapped phases `phi_e`, and magnitudes `m_e`:

1. **MCPC-3D-S phase-offset correction (§2.1.2, Eq. 4).**
   Estimate the zero-echo phase offset from the first two echoes. ROMEO spatially unwraps the wrapped phase difference. Subtract the estimated offset from every echo before multi-echo unwrapping.

2. **Multi-echo ROMEO unwrapping (§2.1.2).**
   Unwrap all echoes together so later echoes are constrained by approximately linear phase accumulation.

3. **Temporal phase correction (§2.1.3, Eqs. 5–6).**
   For each frame, group frames whose magnitude images correlate at least 0.98. Correct the first echo to the nearest `2*pi` branch of the group mean. Correct later echoes to the branch predicted by the previously corrected echoes.

4. **Weighted field-map regression (§2.1.4, Eq. 7).**
   Fit phase against echo time through the origin using squared magnitude weights:

   ```text
   omega = sum(m_e^2 * t_e * phi_e) / sum(m_e^2 * t_e^2)   [rad/time]
   field_hz = omega / (2*pi)
   ```

5. **Low-rank approximation (§2.1.5, Eqs. 8–9).**
   Reshape the field-map series to `Nvox × T` and compute its rank-10 truncated SVD. The paper does not state that the matrix is centered; do not center it unless a black-box experiment proves that convention.

6. **Displacement and inversion (§2.1.6).**
   Convert the native-space field map to an EPI displacement using total readout time and phase-encoding voxel size. Invert the field with behavior equivalent to ITK's `InvertDisplacementFieldImageFilter`, yielding a pull transform and a field map on the undistorted grid.

The paper does not fully specify masking, phase scaling, temporal grouping details, border filtering, displacement signs, inversion parameters, interpolation, extrapolation, or output headers. These are M0 measurements, not assumptions.

## 4. Reuse of the completed ROMEO feature

### 4.1 What is already available

`romeo.c` already supplies:

- raw phase rescaling;
- `robustmask`, `qualitymask`, file masks, and no mask;
- all ROMEO weight presets and quality maps;
- single-echo and multi-echo unwrapping;
- template, individual, global, and temporal-uncertain modes;
- B0 output in Hz;
- `phase_var` weighting.

For echo times in milliseconds, ROMEO's `phase_var` formula is algebraically the MEDIC weighted regression:

```text
B0 = (1000 / 2*pi) *
     sum(m_e^2 * t_e * phi_e) /
     sum(m_e^2 * t_e^2)
```

This is formula reuse, not proof of complete MEDIC equivalence. MEDIC performs phase-offset and temporal corrections before this regression, and exact float-width conventions may differ. Component-level validation remains required.

### 4.2 Required refactor

The current public entry point, `romeo_run()`, is file-oriented and owns loading, side-output creation, and replacement of `nim->data`. MEDIC needs an in-memory, one-frame interface.

Refactor without changing existing `-romeo` results:

```text
romeo_prepare_frame(...)
romeo_unwrap_frame(...)
romeo_compute_b0_frame(...)
romeo_frame_cleanup(...)
```

The names are illustrative. The interface should:

- accept caller-owned float32 echo-major phase and magnitude buffers;
- accept validated dimensions, echo times, and ROMEO options;
- return caller-owned or explicitly owned unwrapped phase, mask, and optional B0;
- perform no file I/O and create no side outputs;
- keep all strict-FP, numeric-width, masking, and ownership rules in one implementation;
- leave `romeo_run()` as a thin adapter over the same kernels.

Do not expose every internal helper. One context type and two or three operations are enough.

MCPC-3D-S belongs beside the strict-FP ROMEO preprocessing code because its corrected phase controls ROMEO's integer edge weights. Port only the required monopolar path from the pinned MIT MriResearchTools source. Preserve its license notice and validate every intermediate against Julia. Do not add bipolar or multi-channel support in v1.

## 5. Architecture

### 5.1 Modules

- `src/medic.c`, `src/medic.h`
  - MEDIC orchestration;
  - temporal correction;
  - weighted field-map driver;
  - low-rank filtering;
  - field-to-displacement conversion;
  - scalar displacement inversion;
  - scalar-axis resampling.

- `src/romeo.c`, `src/romeo.h`
  - existing ROMEO kernels;
  - narrow in-memory frame API;
  - MCPC-3D-S phase-offset correction.

- `src/nifti_io.c`, `src/nifti_io.h` — **not needed, never added.** The opaque sequential float32 volume reader/writer planned here belonged to the superseded streaming design (§5.2); `medic.c` uses the ordinary `nifti_image_read`/`nifti_save` path through its own `md_read_f32()`/`md_write()` helpers.

- `medic.py`
  - BIDS discovery and JSON parsing;
  - subprocess orchestration only.

Keep validation at the public `--medic` boundary. Internal kernels may assume the context has already validated dimensions, counts, echo times, transforms, and allocation sizes. Avoid repeating defensive checks at every layer.

### 5.2 Streaming and bounded memory — **SUPERSEDED, not implemented**

> Retained for history only. Streaming is a decided **non-goal** — see the Status section above: a 4D `.nii.gz` cannot be seeked, so every tool including the reference holds whole volumes in RAM anyway. What shipped instead is the documented, printed working-set budget and the fail-atomic output rule below. Do not implement the API described here without reopening that decision.

Loading every 4D echo at once defeats the design. The 170-frame, two-echo demo already requires hundreds of megabytes for the four inputs; a five-echo, 600-frame run can require tens of gigabytes.

Add a narrow streaming API that:

- reads the NIfTI header once;
- keeps the payload open;
- reads the next 3D volume as scaled, host-endian float32;
- writes a float32 header and appends complete 3D volumes;
- detects short reads/writes;
- supports sequential gzip access without reopening the stream;
- uses the existing zstd temporary-file behavior;
- closes cleanly on every error.

Do not expose `NIIFILE` or compression internals. Use an opaque stream handle.

Processing layout:

1. Open all magnitude and phase inputs.
2. Read a small block of corresponding frames in lockstep.
3. Run MCPC and ROMEO independently for each frame.
4. Write unwrapped phases, masks, and native field maps to uncompressed scratch files.
5. Run temporal correction and low-rank processing from scratch.
6. Write final outputs sequentially.
7. Remove scratch files only after all outputs close successfully.

Outputs should be fail-atomic: write sibling temporary files and rename them only when the complete command succeeds. Scratch files go in a user-selectable directory and are removed on handled failure. Do not hold all corrected phases or field maps in RAM.

The block size is a memory budget, not necessarily the thread count. Default it from volume size, echo count, and a conservative fixed working-set target; allow `--block-frames` for benchmarking.

### 5.3 Parallelism

Parallelize independent frames during MCPC/ROMEO. Do not create nested OpenMP teams. Either:

- run the frame loop in parallel and keep ROMEO's internal regions serial, or
- run frames serially and use ROMEO's internal parallelism.

Benchmark both on the 170-frame reference before choosing. Thread-count parity is required within documented float tolerances; bit identity is desirable but not a requirement for SVD reductions.

**Resolved: the first option shipped** — the per-frame MCPC/ROMEO loop is the `#pragma omp parallel for`, ROMEO's own regions stay serial inside it, and no nested teams are created. Measured scaling 1→8 threads is 3.47× (manifest §5.3).

Do not compile all of `medic.c` strict-FP merely because ROMEO requires it. Keep branch-sensitive phase preprocessing in the existing strict-FP ROMEO translation unit. Compile SVD and resampling under the normal project policy unless measurement identifies a real correctness issue.

### 5.4 Low-rank implementation

Do not plan to reuse `tensor.c` as a general eigensolver. In the default build its active readable implementation is specialized to 3×3 even though a function accepts `n`. That still holds: `md_lowrank()` carries its own `md_jacobi_eigh()`.

**What shipped, replacing the on-disk design below.** The whole field-map series is already resident (§5.2 is a non-goal), so the `T×T` Gram matrix is accumulated in one in-memory pass, diagonalised by a deterministic Jacobi eigensolve, and the series is projected onto the leading `r` eigenvectors in place. Memory is `O(T^2)` **independent of voxel count** — no tiling, no scratch files, no Lanczos. The truncation is **uncentered** (measured, manifest §3.8). `md_lowrank()` scratch is preflighted before `F` is mutated, so a mid-way OOM cannot leave a half-filtered series.

> **Superseded, retained for history.** The original design assumed the series lived on disk: store the native field maps uncompressed and frame-major; read spatial tiles across all frames with bounded random reads; accumulate the symmetric `T×T` Gram in double; solve the largest `r` eigenpairs with a deterministic block-Lanczos or subspace iteration; accumulate `Nvox × r` coefficient volumes; reconstruct one output frame at a time — `O(T^2 + Nvox*r + block*T)`. None of that was built, because the streaming premise it rested on was dropped.

Before committing to a solver, compare it against a trusted offline SVD on synthetic matrices with clustered singular values, rank deficiency, constant frames, and `T < 10`. If convergence or runtime is poor, use a separate permissively licensed small symmetric eigensolver rather than enlarging `tensor.c` with another hidden mode.

Rank is `min(10, T, positive numerical rank)`; `--rank 0` disables the filter. A non-finite field series is a hard error raised on the series itself immediately after the regression, **not** inside `md_lowrank()` — that function returns early whenever `rank >= T`, which would make the guard frame-count dependent.

### 5.5 Python boundary

`medic.py` uses only:

```text
argparse, json, os, pathlib, re, shutil, subprocess, sys
```

It does not import numpy, nibabel, pybids, scipy, or Warpkit. niimath performs all NIfTI I/O.

V1 wrapper behavior:

- accept a BIDS `func/` directory or a dataset root;
- discover `*_echo-<N>_part-{mag,phase}_bold.nii[.gz|.zst]`;
- group files after removing only `echo` and `part` entities;
- require identical numeric echo sets for magnitude and phase;
- ignore `sbref` and non-`bold` suffixes;
- read the exact phase sidecar for each echo;
- convert BIDS `EchoTime` seconds to milliseconds;
- require consistent `TotalReadoutTime` and `PhaseEncodingDirection` across echoes;
- if `TotalReadoutTime` is absent, derive it only when both `EffectiveEchoSpacing` and `ReconMatrixPE` are present;
- call `niimath --medic` once per run;
- call `-unwarp` once per magnitude echo;
- support `--dry-run`, `--niimath`, `--n-cpus`, and `--overwrite`.

**As shipped** `medic.py` takes `input`, `--out-dir` (required), `--niimath`, `--n-cpus`, `--noise-frames`, `--rank`, `--dry-run`, `--overwrite` and `--no-apply` (estimate only, skip `-unwarp`). There is no `--scratch-dir`: it belonged to the superseded streaming design and nothing writes scratch files. There is no `--jobs` either, as planned below.

Do not claim full BIDS inheritance support in v1. Exact sidecars are present in both supplied datasets. Add inheritance only if a real target dataset needs it.

Run different BIDS runs sequentially by default. A `--jobs` option plus OpenMP would create easy oversubscription and is not needed initially.

## 6. CLI and output contract

### 6.1 `--medic`

Required:

```text
--magnitude <one file per echo>
--phase <one file per echo>
--te-ms <comma-separated echo times>
--total-readout-time <seconds>
--phase-encoding-direction <i|j|k|i-|j-|k->
--out-prefix <path>
```

Useful controls:

```text
--n-cpus <N>
--noise-frames <N>
--rank <N>                 default 10; 0 disables
--temporal-correction <0|1>
--phase-offset <mcpc|none>
--mask <file>              external mask, used verbatim by every stage; in-mask is `>= 1` (see section 7.2)
--weights <sel>            ROMEO weight preset: romeo|romeo2|romeo3|romeo4|romeo6 (governs BOTH unwrapping stages)
--save-intermediates
--gz <0|1>
```

`--block-frames` and `--scratch-dir` were part of the superseded §5.2 streaming design and were never implemented; see the Status section.

Avoid exposing speculative border-filter controls until M0 establishes their semantics and importance.

Outputs:

```text
<prefix>_fieldmaps_native.nii[.gz]   Hz, distorted grid
<prefix>_displacementmaps.nii[.gz]   mm, pull map on output grid
<prefix>_fieldmaps.nii[.gz]          Hz, undistorted grid
```

`--save-intermediates` additionally writes per-echo unwrapped phase and masks for component validation. It is not the default because those files are large.

A diagnostic `--medic-from-unwrapped` subcommand may be added during development if it materially simplifies component comparison. Do not commit it as public API until its maintenance value is demonstrated.

### 6.2 `-unwarp`

```text
niimath <input> -unwarp <displacement-map> <axis> <output>
```

Contract:

- scalar displacement in millimeters;
- one voxel axis only;
- output grid equals the input grid in v1;
- a 3D map broadcasts over a 4D input;
- a 4D map must have either one frame or the same frame count as the input;
- map and input must have matching spatial dimensions and world frame;
- out-of-FOV fill and interpolation are set by M0 black-box measurements;
- no Jacobian modulation unless M0 proves otherwise;
- reject oversized inputs; this is not a huge-image-safe operation.

The signed phase-encoding direction may be needed only while generating the displacement map. M0 must determine whether the stored map already contains the sign. Do not make `-unwarp` negate it a second time. **Answered (manifest §3.5/§3.5b):** the stored map already carries the sign, so `-unwarp` accepts a trailing `-`/`+` and ignores it — while `--medic --phase-encoding-direction` **honours** it, because the polarity drives the inversion and the displacement sign. The two are deliberately different and both are measured.

## 7. Unknown conventions and decisive experiments

Run these before implementing the affected component:

1. **Phase scaling.**
   Feed integer phase with known stored range and scaling. Determine whether the reference uses header scaling, observed extrema, a fixed scanner range, or radians declared in JSON.

2. **Mask generation and use.**
   Capture reference masks with public debug/intermediate output. Determine *what the mask controls* — unwrapping, regression, SVD, border processing, output zeroing, or some subset — and whether one mask is shared across stages.

   **Do not attempt to reproduce the reference's mask exactly.** Both implementations use crude, heuristic masks and neither is authoritative, so bit-matching one to the other buys nothing scientific and is an open-ended reverse-engineering task. The requirement is instead:

   - `--mask <file>` accepts an **external mask**, used verbatim by every stage that needs one. **In-mask means `>= 1`**, the measured reference contract (manifest §3.7) — not "non-zero". A fractional probability map must be thresholded first; NaN is treated as outside; a mask with no voxel `>= 1` is an error, not an empty run.
   - That option is what makes exact cross-validation possible: supply the *same* mask to both implementations and any remaining difference is a real algorithmic difference, not a masking difference.
   - It also lets users supply a **better** brain mask than either tool's built-in heuristic (e.g. mindgrab), which is the more useful capability in practice.
   - The built-in default remains ROMEO's `robustmask`; no attempt is made to match the reference's own construction.

   Consequently, parity gates below are evaluated **with a shared supplied mask**. Differences attributable solely to mask choice are reported, not chased.

3. **MCPC smoothing.**
   Compare raw Eq. 4 and the pinned MriResearchTools MCPC-3D-S result against reference unwrapped phases.

4. **Temporal grouping.**
   Create a short series with controlled magnitude correlations and injected `2*pi` shifts. Determine which echo's magnitude defines correlation, how ties and empty groups behave, whether correction is per voxel or a global frame offset, and the rounding direction.

5. **Low-rank details.**
   Use small synthetic field-map matrices where centered and uncentered SVD differ. Probe rank zero, `T < 10`, masks, non-finite values, and frames excluded as noise.

6. **Border filter.**
   Ablate the public border-filter option if available. Localize changed voxels and quantify its effect on displacement inside the brain. Implement it only if it materially affects the target outputs.

7. **Hz-to-mm sign and units.**
   Use constant ±field maps, anisotropic voxels, all six BIDS phase-encoding directions, and permuted/flipped voxel orientations.

8. **Displacement inversion.**
   Use constant shifts, affine-varying shifts, smooth nonlinear shifts, and a deliberately folding field. Determine iteration tolerance, boundary behavior, sign, output domain, and failure behavior.

9. **Interpolation and fill.**
   Resample an impulse, a half-voxel Gaussian, a ramp, and a constant image. Distinguish linear, B-spline, and windowed-sinc interpolation; distinguish zero, clamp, mirror, and nearest boundary behavior.

10. **Jacobian modulation.**
    Resample a constant image through a nonuniform but invertible field. Any intensity variation indicates modulation.

11. **Headers and naming.**
    Record all three reference output headers and determine which stage is represented by `_fieldmaps_native` versus `_fieldmaps`.

12. **Noise frames.**
    Confirm whether trailing noise frames are omitted from all outputs or excluded only from estimation.

Each experiment should be a small script or documented command with an analytic expected result. Do not retain a conclusion that cannot be reproduced.

## 8. Milestones and gates

### M0 — Freeze the observable contract

- Inventory existing `demo/out`, `demo/out1`, and `demo/out170` outputs.
- Regenerate only when the command/version/input hashes differ.
- Run the experiments in §7.
- Capture the two-echo 170-frame run and the supplied three-echo dataset.
- Record Warpkit timing and peak RSS.

**Gate:** every implementation-sensitive convention is either measured and documented or explicitly deferred behind a non-default option. No code based on guessed sign, scaling, interpolation, or inversion behavior.

### M1 — In-memory ROMEO API

- Extract the narrow frame interface described in §4.2.
- Keep `romeo_run()` behavior and all existing outputs unchanged.
- Add ownership and failure-path tests for the new API.

**Gate:** existing ROMEO parity remains 602/602; real CLI outputs remain equal at their current thresholds; repeated in-memory calls have clean UBSan, malloc diagnostics, and `leaks`.

### M2 — Sequential NIfTI volume I/O — **DROPPED with §5.2**

Not implemented and not needed: the streaming premise was retired (see Status), so `medic.c` reads whole images through `md_read_f32()` (header preflight, then payload, converting datatype **or** scaling) and writes through `md_write()`. The original text — opaque scaled-float32 reader/writer handles, `.nii`/`.nii.gz`/byte-swapped/short-input tests, and a gate of memory independent of frame count — is retained here only as history.

### M3 — `-unwarp`

- Implement scalar-axis pull resampling using M0's measured sign, interpolation, and fill.
- Test 3D, 4D frame pairing, 3D broadcast, anisotropic voxels, and orientations.

**Gate:** using Warpkit's displacement map as input, reproduce its corrected magnitude output with normalized RMSE and absolute-error percentiles below the agreed threshold. Correlation alone is insufficient.

### M4 — Single-frame MEDIC estimate

- Implement MCPC-3D-S from the pinned MIT source.
- Run in-memory ROMEO.
- Reuse the existing `phase_var` B0 kernel.
- Write native field map and validation intermediates.

**Gate:** MCPC intermediates match Julia at the established ROMEO numeric policy; reference unwrapped phases differ only by tolerated rounding and integer `2*pi` branch choices; native field maps meet absolute Hz error thresholds inside the common mask.

### M5 — Displacement conversion and inversion

- Implement the measured Hz-to-mm convention.
- Implement scalar-axis inversion with behavior equivalent to the reference/ITK result.
- Resample the native field map onto the undistorted grid.

**Gate:** constant and analytic fields pass exact/property tests; reference displacement error is below 0.05 mm at the 95th percentile inside the valid mask; non-finite mismatch count is zero.

**Status: MET as a convention check, NOT met end to end for `j`.** The inversion and Hz→mm formulas, fed the reference's own native field, reproduce the reference's displacement map at p95 0.006 mm (manifest §3.4, `j`) and 0.045/0.024 mm (§3.5b, `j`/`j-`). The full `--medic` pipeline with the reference mask supplied is at p95 **0.059 mm for `j`** — over the gate — and 0.029 mm for `j-`; with the shipping `robustmask` default it is 2.40 mm, which is the mask difference of §7.2 and is not being chased. Manifest §5.4 attributes each figure. Non-finite mismatch count is zero throughout.

### M6 — Multi-frame estimator (**not** streaming; §5.2 dropped)

- Process every frame of every echo in one resident working set, repacked frame-major/echo-minor.
- No scratch storage: unwrapped phases, masks and native field maps stay in RAM, and `--save-intermediates` writes them only on request.
- Establish the chosen OpenMP structure: the frame loop is parallel and ROMEO's internal regions stay serial.

**Gate:** results match independent single-frame M4/M5 runs; peak RSS stays within the documented budget on the 170-frame run; `-p 1`, `-p 2`, and `-p 8` agree within the documented tolerance.

### M7 — Temporal phase consistency

- Implement the paper-defined correction after M0 resolves its remaining conventions.
- Read the phase from the resident working set (the on-disk/tiled access this bullet used to require went with §5.2), taking group means from an immutable snapshot of the first-echo series so the result cannot depend on traversal order.
- Optimize repeated or nearly identical correlation groups only after profiling.

**Gate:** injected `2*pi` errors are removed and untouched frames remain unchanged; the real runs agree with the reference on which frames/echoes are corrected and on the resulting field maps.

### M8 — Low-rank and border cleanup

- Implement and validate the rank-limited solver in §5.4.
- Add border processing only if M0 shows material impact.

**Gate:** synthetic matrices meet residual and subspace-angle tolerances against a trusted offline SVD; reference denoised field maps meet absolute Hz and displacement-mm thresholds; runtime and memory are recorded (there is no scratch usage — nothing is written to disk).

### M9 — End-to-end `--medic`

- Connect frame estimation, temporal correction, low-rank filtering, inversion, and output writing.
- Make output creation fail-atomic — as shipped, sibling-temp-plus-`rename` via `md_write_temp()`.
- Add clear stage-specific errors and concise progress reporting.

**Gate:** the single-frame demo, 170-frame two-echo run, and three-echo run complete end to end. Compare native field, displacement, undistorted field, and corrected magnitudes separately. Report finite mismatch counts, median/95th/max absolute errors, normalized RMSE, and spatial correlation.

### M10 — Minimal BIDS wrapper

- Implement §5.5.
- Test file grouping, numeric echo ordering, missing pairs, inconsistent metadata, dry-run, overwrite, and paths containing spaces.

**Gate:** the wrapper discovers and processes both supplied BIDS layouts without nibabel/numpy and produces the same niimath commands as their manual equivalents.

### M11 — Build, CI, documentation, and diagnostics

- Add `HAVE_MEDIC`, `MEDIC=0`, and `ENABLE_MEDIC`, requiring `HAVE_ROMEO`.
- Cover Makefile, CMake/SuperBuild, notarization, native release, debug, UBSan, no-OpenMP, and supported wasm inventories.
- Decide explicitly whether MEDIC belongs in tiny/nano/WASI; do not add it mechanically.
- Add citations and clean-room provenance to help and documentation.
- Add only compact synthetic tests to CI.

**Gate:** normal regression tests pass; serial UBSan passes; `MallocScribble=1 MallocGuardEdges=1` passes the small end-to-end fixture; Darwin `leaks` is clean; build-off configurations do not reference MEDIC symbols.

### M12 — Performance

**Build against zlib-ng (or zlib-cloudflare), not system zlib.** Plain `make` links `-lz`; use `make -C src ZLIBNG_ROOT=<path-to-zlib-compat-build>`, or the CMake release path, which defaults to `ZLIB_IMPLEMENTATION=zlib-ng`. On the 170-frame run this is not a rounding error: gzip of the three float32 output series dominates the serial tail, and switching to zlib-ng took the estimate stage from 16.2 s to 10.6 s wall. Any timing comparison made with a system-zlib build understates niimath by ~35 % and should be rejected.


- Profile decompression, ROMEO, temporal correlations, SVD, inversion, and output compression separately.
- Compare wall time and peak RSS with the recorded Warpkit baseline.
- Optimize measured hotspots only.

**Gate:** memory accounted for on the target five-echo/600-frame geometry by calculation and on the 170-frame dataset by measurement. No performance optimization may weaken the component error gates. **Note what this gate can and cannot say now that streaming is a non-goal:** the five-echo/600-frame calculation (~10 GiB of work arrays) is a statement of cost, not a bound — there is no mechanism that keeps a large run inside a fixed budget, which is exactly why `--medic` is documented as native-scale-only.

## 9. Validation policy

Use four layers:

1. **Analytic properties**
   - known phase slope;
   - known displacement;
   - warp/inverse round trip;
   - orientation equivariance;
   - constant-image resampling;
   - injected `2*pi` temporal errors;
   - exact low-rank matrices.

2. **MIT upstream parity**
   - existing ROMEO suite;
   - new MCPC-3D-S intermediates against pinned MriResearchTools.

3. **Black-box component comparison**
   - unwrapped phase and masks;
   - native field maps;
   - displacement maps;
   - undistorted field maps;
   - corrected magnitudes.

4. **End-to-end real data**
   - supplied single-frame, 170-frame/two-echo, and three-echo datasets.

`niimath --compare` remains useful, but do not use its exit status alone:

- count finite/non-finite mismatches separately;
- compare only a documented common-validity mask;
- report absolute errors in physical units;
- report normalized RMSE and correlation for corrected magnitude;
- inspect output headers independently;
- keep tolerances tied to displacement error, preferably well below 0.1 voxel.

Float32 is the target storage type. SIMD/FMA differences are acceptable only when the physical-unit gates pass. Integer `2*pi` branch errors are not ordinary rounding error and must be reported separately.

## 10. Failure and ownership rules

- Validate all input headers before reading payloads.
- Require equal spatial dimensions, frame count, and agreeing world transforms across all echoes and parts.
- Require finite, positive, strictly increasing echo times.
- Require at least two echoes and at least one retained frame.
- Validate allocation products with existing checked-size helpers.
- Assign only plain `malloc`/`calloc` buffers to `nim->data`.
- Every stream/context owns an explicit close/free operation.
- A failed stage writes no final output and does not leave apparently valid partial NIfTIs.
- Treat stdin as unsupported for `--medic`; multiple synchronized inputs make its semantics ambiguous.

Keep these checks at orchestration boundaries. Kernels should not duplicate them.

## 11. Non-goals

- Reading or translating Warpkit implementation or tests.
- Coil combination or 5D coil-channel input.
- Bipolar MCPC correction.
- Unwrapping phase after resampling.
- Correcting wrapped phase images.
- Motion correction, slice timing, TOPUP, or the paper's downstream fMRI preprocessing.
- General 3-vector ITK/FSL/ANTs/AFNI displacement-field conversion.
- Jacobian images unless required for measured MEDIC behavior.
- Full BIDS validator or metadata inheritance in v1.
- Exact bit identity with Warpkit.
- New external runtime dependencies.

## 12. First implementation sequence

1. Complete M0 and commit the reference manifest.
2. Refactor the in-memory ROMEO API with no behavioral change.
3. ~~Add streaming NIfTI volume I/O~~ — dropped with §5.2/M2; whole-image reads are used instead.
4. Implement and validate `-unwarp` using reference displacement maps.
5. Implement single-frame MCPC + ROMEO + weighted B0.
6. Add displacement conversion and inversion.
7. Scale to the multi-frame resident estimator.
8. Add temporal correction, then low-rank filtering.
9. Join the stages under `--medic`.
10. Add the minimal Python wrapper last, when the C command contract is stable.
