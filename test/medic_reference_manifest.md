# MEDIC reference manifest (M0)

Everything below was obtained by treating Warpkit as a **black box**: running the installed `wk-medic` / `wk-apply-warp` executables and measuring their outputs. No Warpkit implementation, test, build product, or debug symbol was read. Provenance for every row is marked **P** (specified by the paper, `warpkit.pdf`, doi:10.1162/IMAG.a.1262), **U** (specified by the MIT upstream ROMEO.jl / MriResearchTools.jl already ported in `src/romeo.c`), or **B** (observed only from the black box).

The reference NIfTIs are large and gitignored (`/test/medic_ref/`). This file plus the scripts in `test/medic_experiments/` are the committed record; every conclusion here is reproducible by re-running the named script.

## 1. Environment

| item | value |
| --- | --- |
| Warpkit executables | `~/src/warpkit/.venv/bin/wk-medic`, `wk-apply-warp` |
| version string | `wk-medic 1.4.1` |
| host | Darwin 26.5.2, Apple Silicon |
| reference oracle python | `~/src/warpkit/.venv/bin/python` (numpy; **never** `import warpkit`) |
| niimath binary under test | `src/niimath`, `v1.0.20260724 OpenMP Clang17.0.0 BSD (64-bit MacOS)` |

Experiment scripts (`test/medic_experiments/`, analysis-only, not shipped, not in CI):

| script | covers |
| --- | --- |
| `nii.py` | minimal NIfTI-1 read/write helper |
| `exp09_interp_fill.py` | §7.9 interpolation, §7.9 fill, §7.10 Jacobian, sign |
| `exp07_axis_orientation.py` | §7.7 letter → voxel axis vs world axis |
| `exp07b_axis_sweep.py` | §7.7 full 10-grid × 3-letter physical-displacement sweep |
| `exp01_08_known_field.py` | §7.1 phase scaling, §7.8 inversion, analytic field |
| `exp04_05_12_multiframe.py` | §7.4 temporal, §7.5 low-rank, §7.12 noise frames |

## 2. Reference run inventory

Regenerate only if a hash below changes. Command lines are `~/src/warpkit/demo/run.sh` and `run170.sh` verbatim.

### Inputs

| sha256 (first 16) | bytes | file |
| --- | --- | --- |
| `ba23f0164a13c567` | 449997 | `demo/data/echo-1_part-mag.nii.gz` |
| `62e949ba4d5d0a58` | 449206 | `demo/data/echo-1_part-phase.nii.gz` |
| `20bfa56042703ab4` | 426550 | `demo/data/echo-2_part-mag.nii.gz` |
| `ec9b7d86366bbdc5` | 450552 | `demo/data/echo-2_part-phase.nii.gz` |
| `d32479cfacc8bf8c` | 75271515 | `inputs/.../echo-1_part-mag_bold.nii.gz` |
| `15c034104766dea2` | 76275776 | `inputs/.../echo-1_part-phase_bold.nii.gz` |
| `c8d2edfc482b382a` | 70453558 | `inputs/.../echo-2_part-mag_bold.nii.gz` |
| `b77b3fd11f3ef7f9` | 76530748 | `inputs/.../echo-2_part-phase_bold.nii.gz` |

Geometry: 76×76×46, pixdim `(-1, 2.80263, 2.80263, 2.8)`, `xyzt_units` 10 (mm/s), qform=sform=1, oblique (det(srow) = −21.993, LAS-ish). sbref = 1 frame, BOLD = 170 frames. TEs 16.8 / 38.56 ms, `TotalReadoutTime` 0.02025 s, `PhaseEncodingDirection` `j`. Phase stored uint16 with `scl_slope=2, scl_inter=-4096` → scaled range −4096…4094.

### Outputs

| sha256 (first 16) | bytes | file |
| --- | --- | --- |
| `265dd4002f71a83f` | 531744 | `demo/out/sub-fm_fieldmaps_native.nii` |
| `a57db07f394f3c56` | 531744 | `demo/out/sub-fm_fieldmaps.nii` |
| `1373e354b6256e67` | 531744 | `demo/out/sub-fm_displacementmaps.nii` |
| `33e1dca3be4b6c56` | 461008 | `demo/out/echo-1_part-mag_undistorted.nii.gz` |
| `41d6a6824b7c04be` | 437810 | `demo/out/echo-2_part-mag_undistorted.nii.gz` |
| `ea27901938436c84` | 90336992 | `demo/out170/sub-fm_fieldmaps_native.nii` |
| `c6f9d3ea33ff9ef3` | 90336992 | `demo/out170/sub-fm_fieldmaps.nii` |
| `8d34e966d4bfcff3` | 90336992 | `demo/out170/sub-fm_displacementmaps.nii` |
| `d426f6d896fbb47f` | 76991561 | `demo/out170/echo-1_part-mag_bold_undistorted.nii.gz` |
| `4d02098130c34c4a` | 71998767 | `demo/out170/echo-2_part-mag_bold_undistorted.nii.gz` |

`out1/` is a partial rerun of `out/` (echo-2 only) and carries byte-identical map files — confirming `wk-medic` is deterministic across runs on identical input.

Timing (sbref, 1 frame, `--debug`, 4 CPUs): 0.71 s wall, 111 % CPU. The 170-frame run is documented by the demo README as ~1 min at `NCPUS=10`; see M12 for our own measurement.

### Output headers — B

| output | datatype | scl_slope / scl_inter | qform / sform |
| --- | --- | --- | --- |
| `_fieldmaps_native` | uint16 (512) | 0.00534369 / −162.100 | 1 / 1 (copied from input) |
| `_fieldmaps` | uint16 (512) | 0.00499006 / −150.629 | **0** / 2 |
| `_displacementmaps` | uint16 (512) | −0.000283202 / 8.54869 | **0** / 2 |

The `srow` is byte-identical across all three and equal to the input's; only `qform_code`/`sform_code` differ, so the "undistorted grid" **is the input grid** — no regridding occurs.

**Output datatype follows the PHASE input datatype** (B, `exp01_08_known_field.py`): float32 phase in → float32 maps out (`scl 1/0`); uint16 phase in → uint16 maps out with computed scaling. The magnitude datatype does not affect it.

> **niimath decision:** we write **float32** for all three outputs (plan §6.1). uint16 + `scl` is lossy — the quantum is 5.3e-3 Hz / 2.8e-4 mm on the demo — and every gate below is stated in physical units, not stored codes.

## 3. Settled conventions

### 3.1 Phase scaling — §7.1 — B (agrees with U)

`wk-medic` logs `Estimated min phase: -4096.0 / Estimated max phase: 4094.0` and rescales the **observed min/max of the scaled data** linearly onto [−π, π]. Writing the identical phase three ways — float32 radians, uint16 with `scl 2/−4096`, uint16 with `scl 1/0` — produced field maps agreeing to uint16 quantization (`exp01_08_known_field.py`). The header `scl_slope`/`scl_inter` is therefore **not** consulted beyond producing the scaled values whose extrema are taken.

This is exactly ROMEO's `readphase`, already implemented and Julia-validated in `src/romeo.c`. **No new code needed.**

### 3.2 Weighted field-map regression — §3.4 — P, confirmed B

```text
omega    = sum(m_e^2 * t_e * phi_e) / sum(m_e^2 * t_e^2)
field_Hz = omega / (2*pi)
```

Reproducing this from `--debug`'s own `phase0.nii`/`phase1.nii` and the input magnitudes matches `_fieldmaps_native` to **max 0.0052 Hz** — one uint16 quantum (0.00534). Exact. `t_e` is in **seconds**; the regression runs over **all** voxels (no mask term), and is zero outside the mask only because the unwrapped phase is zero there.

This is algebraically ROMEO's `phase_var` B0 mode, already in `src/romeo.c`.

### 3.3 Displacement: Hz → mm — §7.7 — B

Exact identity on the demo outputs (`exp08_inversion.py`), residual 3.8e-7 mm at p95 versus a 5.0e-3 quantum:

```text
displacement_mm = -field_undistorted_Hz * TotalReadoutTime * pixdim[PE axis]
```

Note it is the **undistorted-grid** field, not `_fieldmaps_native`, and the length scale is the NIfTI `pixdim` of the phase-encoding voxel axis (2.80263 mm here), not `|srow[:,j]|` — they coincide on this data.

### 3.4 Displacement inversion — §7.8 — B

The three outputs are mutually consistent with a single scalar fixed point along the PE voxel axis (units: voxels):

```text
f_undistorted(y) = f_native( y + f_undistorted(y) * TotalReadoutTime )
```

i.e. the undistorted-grid field is the native field sampled at the **distorted** location. Solving this by direct iteration from a zero start with **linear** interpolation along the PE axis reproduces `_fieldmaps` on the demo data to p50 = 1e-4 Hz, p95 = 0.041 Hz. Cubic/nearest sampling are strictly worse (p95 2.05 / 4.07), so the field resampling is **linear**, unlike `wk-apply-warp` (§3.6).

Two honest caveats:

- **We do not reproduce Warpkit's inverter exactly, by design.** It is iteration-limited (ITK-style), not converged: on a synthetic linear ramp with analytic answer `1/(1 − b·TRT) = 1.14490`, Warpkit returns a ratio of `1.14911` (+0.37 %), which is *above* the converged fixed point and so cannot be reached by running our iteration longer or shorter. On real data 1.8 % of voxels differ by >1 Hz, concentrated where `|f_native|` is large (p50 53 Hz there vs 0 Hz overall) — i.e. near folds where the inverse is genuinely multi-valued.
- **It passes the plan's M5 gate with 8× margin anyway.** Converged fixed point vs Warpkit's `_displacementmaps`, inside a crude magnitude-based brain mask (106 268 voxels): p50 = 0.0002 mm, **p95 = 0.006 mm** (gate: < 0.05 mm), p99 = 1.4 mm, max 9.0 mm — the tail being exactly the folded voxels. Displacement range is −10.0…8.5 mm.

**niimath decision:** implement the converged fixed point (iterate to tolerance with a cap). It is the mathematically correct inverse, deterministic, and ~15 lines. Gate on displacement p95 inside the common mask, per plan §9.

### 3.5 `--phase-encoding-axis` semantics — §7.7 — B

This one is not guessable and was worth the sweep. Measured over 10 grids × 3 letters (`exp07b_axis_sweep.py`), plus the real oblique data, with a **constant** map to remove every other ambiguity.

The letter names a **voxel index axis** (`i`/`x`→0, `j`/`y`→1, `k`/`z`→2). But the displacement applied is a **physical vector along the canonical world axis that voxel axis is most aligned with** — *not* along the image's own column direction. Writing `A` for the 3×3 `srow` (voxel→world RAS), `u = A[:,m]/|A[:,m]|`, `w = argmax|u|`, `sigma = sign(u[w])`, and `kappa = (-1,-1,+1)` (the RAS→LPS sign of world axis `w`):

```text
delta_RAS = d * sigma * kappa[w] * e_w          # physical displacement, mm
s         = A^-1 @ delta_RAS                    # voxel-space offset
out(v)    = in(v + s)                           # pull
```

All 30 sweep rows match exactly, and on the real oblique demo data this reproduces `echo-{1,2}_part-mag_undistorted.nii.gz` at **nrmse 3.5e-5**. Modelling the displacement along the image's own `j` column instead gives nrmse 4.3e-2 — three orders of magnitude worse — because that column sits 14.4° off world +y here, so the correct model carries genuine `i` and `k` voxel components.

Two consequences:

- **The `-`/`+` suffix is ignored.** `j` and `j-` produce byte-identical output. The sign already lives in the stored map (§3.3 carries the minus). Plan §6.2's open question is answered: **`-unwarp` must not negate again.**
- This is arguably a Warpkit quirk for oblique acquisitions — the physically correct EPI shift is along the voxel PE column, and using the canonical world axis costs a factor `cos(14.4°) = 0.968` plus spurious off-axis components on this data. We implement the measured convention because `--medic` and `-unwarp` must be self-consistent with the reference; recorded here so the choice is deliberate and reversible.

### 3.6 Interpolation, fill, Jacobian — §7.9, §7.10 — B

From the impulse response at a half-voxel shift (`exp09_interp_fill.py`), taps at |x| = 0.5…4.5 measured as
`0.62620, −0.18216, 0.08106, −0.033457, 0.0077309`:

- **Kernel: Lanczos-windowed sinc, radius 5, UNNORMALIZED.** Fitting `sinc(x)·window(x/R)` over R ∈ {4,5,6} and windows {lanczos, hamming, cosine, welch, blackman}, lanczos/R=5 matches to **1.3e-8** (float32 noise); the runners-up are off by 9e-3 to 3.8e-2. The weights sum to 0.998746, not 1 — the filter is not normalized, which is directly visible as a ≤0.13 % dip when resampling a constant image at fractional offsets.
- **Out-of-FOV fill: zero.** A ramp pulled ±5 voxels reads exactly 0 beyond the edge.
- **No Jacobian modulation.** A constant image through a field with 0.1 voxel/voxel gradient comes back flat to within the kernel's own normalization dip; correlation between the observed ratio and the field Jacobian is 0.003.
- The output is **not clipped to the input range** (undistorted magnitude reaches −1654 on non-negative input) — consistent with an unclamped sinc-family kernel.

Because the displacement has components on all three voxel axes for oblique data (§3.5), the kernel is applied **separably in 3D**, not only along the PE axis.

### 3.7 Masking — §7.2 — B, **not reproduced**

`--debug` writes `masks.nii` with three levels: 0 (165 738 voxels), 1 (35 081), 2 (64 877). Magnitude increases monotonically with level (p50 = 291 / 1530 / 11 667), but the levels are not a pure intensity threshold — the ranges overlap, so a spatial step (component labelling or hole filling) is involved.

What the mask **does**: the per-echo unwrapped phase is nonzero exactly on `mask >= 1` (99 925 of 99 958 voxels; the 33 exceptions are voxels whose phase is genuinely ≈0, and `nz \ (mask>=1)` is empty). So the mask gates **where phase is unwrapped**, and hence where the field map is nonzero. It does not appear as a weight in the regression (§3.2 is exact without one).

Hypotheses tested and **rejected**:

| hypothesis | result |
| --- | --- |
| niimath `robustmask` (4D input) == `mask>=1` | dice 0.868 |
| niimath `robustmask` == `mask==2` | dice 0.885 |
| `robustmask(echo1) + robustmask(echo2)` == masks.nii | 88.1 % exact, counts `[183886, 9564, 72246]` vs `[165738, 35081, 64877]` |
| niimath `qualitymask` == `mask==2` | dice 0.415 |
| `robustmask & qualitymask` == `mask==2` | dice 0.889 |

niimath's `robustmask` on the 4D input is byte-identical to `robustmask(echo-1 magnitude)`, confirming ROMEO takes its mask from the first echo.

**niimath decision (deferred, deliberate):** use niimath's existing Julia-validated `robustmask` and compare only on the **common validity mask**, exactly as plan §9 prescribes. Reverse-engineering a 3-level mask that the paper does not specify is not worth the risk of encoding a guess; revisit only if an M9 gate fails because of it.

### 3.8 Low-rank filter — §7.5 — B

Feeding a field series built with 15 nonzero singular values decaying 2× each (T = 20) and reading the singular values back (`exp04_05_12_multiframe.py`):

```text
input  numerical rank : 13   (sv: 2.6e4, 1.1e2, 5.3e1, 2.8e1, 1.5e1, 6.9, 3.6, 1.4, 0.70, 0.41, 0.13, ...)
output numerical rank : 10   (sv 11..20 collapse to ~5e-4, five orders below sv 10)
```

- **Rank-10 truncation IS applied, and it lands in `_fieldmaps_native`** — so the written native field map is already filtered; `_fieldmaps` and `_displacementmaps` derive from the filtered series.
- **Uncentered.** The temporal mean survives but is not preserved exactly (‖mean frame‖ 5917.87 → 5904.88, −0.22 %), which is the signature of truncating the raw matrix. Centering would have preserved the mean exactly on reconstruction. This confirms plan §3.5's default: **do not center.**
- For T = 1 (the sbref demo) rank-10 truncation of a one-column matrix is the identity, consistent with §3.2 matching exactly there.

### 3.9 Temporal phase correction — §7.4 — B, partially characterised

A whole-2π injection into already-wrapped phase is the identity, so it cannot probe anything (my first attempt, corrected). The working probe adds `1/TE_1 = 59.524 Hz` to frame 3's field only — exactly one 2π cycle at echo 1, 2.295 cycles at echo 2:

```text
frame 0,1,2,4..7 : err vs truth p50 0.163 Hz, max 0.368   (untouched)
frame 3          : err vs truth p50 45.95 Hz, max 46.29   median(frame3 - frame0) = +13.71
```

So a genuine 59.5 Hz single-frame excursion is **suppressed to 13.7 Hz** — the correction pulls the frame back toward the group's 2π branch, and does so as a **spatially uniform per-frame offset** (p50 ≈ max across voxels). That is the behaviour the paper's §2.1.3 describes, and it confirms the correction is aggressive: it will flatten real single-frame field changes that happen to sit near a 2π multiple at the first echo.

The grouping threshold (magnitude correlation ≥ 0.98), which echo's magnitude defines it, and tie/empty-group handling remain **paper-specified but not black-box-measured**. Implement per the paper (plan §3.3) and validate end-to-end on the 170-frame run at M7/M9.

### 3.10 Noise frames — §7.12 — B

`-f N` **drops the N trailing frames from the outputs entirely** (T = 8, `-f 2` → 6 output frames), rather than merely excluding them from estimation.

### 3.11 Border filter — §7.6 — not applicable

`wk-medic --help` exposes no border-filter option in 1.4.1 (only `--wrap-limit`, "turns off some heuristics for phase unwrapping"). There is nothing to ablate, so no border processing will be implemented. Plan §8/M8's "add border processing only if M0 shows material impact" resolves to **no**.

### 3.12 Pipeline order — B, inferred from §3.2/§3.4/§3.8

```text
per frame:  rescale phase -> MCPC-3D-S offset -> multi-echo ROMEO unwrap (mask>=1)
            -> weighted regression                       [= raw native field]
across frames: temporal 2*pi correction -> rank-10 truncation
            -> _fieldmaps_native
            -> scalar fixed-point inversion  -> _fieldmaps
            -> * -TRT * pixdim_PE            -> _displacementmaps
```

## 4. What still needs porting

`--debug` also writes `phase_offset0.nii` (range ±π, the MCPC-3D-S zero-echo offset) and `phase{0,1}.nii` (per-echo unwrapped phase). Comparing niimath's current `-romeo` against `phase{0,1}.nii` shows the expected large disagreement — median 4.40 rad at echo 1, with 49 441 of 64 877 in-mask voxels off by a whole 2π — because niimath does **not** yet remove the phase offset before unwrapping. Once offsets are removed the unwrapped phases are near-perfectly linear in TE: `median(phi_2/phi_1) = 2.295230` versus `TE_2/TE_1 = 2.295238`.

**MCPC-3D-S is the one genuinely new numeric kernel** and is the M4 deliverable. Port the monopolar path only, from the pinned MIT MriResearchTools source, beside the strict-FP ROMEO code (plan §4.2).

## 5. M0 gate

Every implementation-sensitive convention is now measured and recorded, or explicitly deferred behind a documented decision:

| convention | status |
| --- | --- |
| phase scaling | measured, already implemented (§3.1) |
| weighted regression | measured exact (§3.2) |
| Hz→mm sign and units | measured exact (§3.3) |
| inversion | measured; converged fixed point adopted, residual quantified, gate passes 8× (§3.4) |
| PE-axis semantics, displacement sign | measured over 30 configurations (§3.5) |
| interpolation kernel, fill, Jacobian | measured exact (§3.6) |
| masking | measured *behaviour*; construction not reproduced — deliberate, gate on common mask (§3.7) |
| low-rank rank and centering | measured (§3.8) |
| temporal correction | behaviour measured; parameters from the paper (§3.9) |
| noise frames | measured (§3.10) |
| border filter | not exposed; not implemented (§3.11) |
| output headers/datatype | measured; niimath deliberately writes float32 (§2) |

No code below rests on a guessed sign, scaling, interpolation, or inversion behaviour. The two open items (§3.7 mask construction, §3.9 grouping parameters) are named, bounded, and carry a validation plan.
