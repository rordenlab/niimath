# MEDIC reference manifest (M0)

Everything below was obtained by treating Warpkit as a **black box**: running the installed `wk-medic` / `wk-apply-warp` executables and measuring their outputs. No Warpkit implementation, test, build product, or debug symbol was read. Provenance for every row is marked **P** (specified by the paper, `warpkit.pdf`, doi:10.1162/IMAG.a.1262), **U** (specified by the MIT upstream ROMEO.jl / MriResearchTools.jl already ported in `src/romeo.c`), or **B** (observed only from the black box).

The reference NIfTIs are large and gitignored (`/test/medic_ref/`). This file plus the scripts in `test/medic_experiments/` are the committed record; every conclusion here is reproducible by re-running the named script.

**No equivalence is claimed with the reference tool, and none should be inferred from any table below.** `-unwarp` reproduces the reference's corrected magnitudes to nrmse 3.5e-5 (§5.1) — that stage does match. `--medic` does **not** match end to end: given a shared mask its native field map agrees to p99 0.0027 Hz (§5), but the mask construction is deliberately not reproduced (§3.7, §4.3), the inverter is deliberately different (§3.4), a broadband low-rank residual is unexplained (§5.2), and the current end-to-end displacement p95 misses the plan's M5 gate for `j` (§5, §5.4). Bit identity with the reference is an explicit non-goal (`medic_plan.md` §11).

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
| `exp08_inversion.py` | §3.3 Hz→mm identity, §3.4 fixed-point inversion (linear vs cubic vs nearest) |
| `exp04_05_12_multiframe.py` | §7.4 temporal, §7.5 low-rank, §7.12 noise frames |
| `bench170.py` | M12 wall time / CPU time / peak RAM head-to-head (gzipped output, plus the apply stage) |
| `bench_threads.py` | M12 like-for-like head-to-head at 1 and 8 threads, both tools writing uncompressed `.nii` |

**Two measurements below have no committed script** and are therefore not reproducible from this directory: the §3.5b polarity sweep (`j` vs `j-`, including its 0.045/0.024 mm figures — this is why its support is unrecorded) and the end-to-end agreement tables of §5, which were run by hand against `test/medic_ref/`. Both are flagged where they appear.

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

i.e. the undistorted-grid field is the native field sampled at the **distorted** location. Solving this by direct iteration from a zero start with **linear** interpolation along the PE axis reproduces `_fieldmaps` on the **sbref demo (1 frame, `j`)** to p50 = 1e-4 Hz, p95 = 0.041 Hz. Cubic/nearest sampling are strictly worse (p95 2.05 / 4.07), so the field resampling is **linear**, unlike `wk-apply-warp` (§3.6).

Two honest caveats:

- **We do not reproduce Warpkit's inverter exactly, by design.** It is iteration-limited (ITK-style), not converged: on a synthetic linear ramp with analytic answer `1/(1 − b·TRT) = 1.14490`, Warpkit returns a ratio of `1.14911` (+0.37 %), which is *above* the converged fixed point and so cannot be reached by running our iteration longer or shorter. On real data 1.8 % of voxels differ by >1 Hz, concentrated where `|f_native|` is large (p50 53 Hz there vs 0 Hz overall) — i.e. near folds where the inverse is genuinely multi-valued.
- **As a CONVENTION check the inversion clears the plan's 0.05 mm M5 threshold with 8× margin.** This row is *not* a `--medic` end-to-end result — see §5.4, which attributes every displacement figure in this file. Measured: the converged fixed point plus the §3.3 scaling, applied to the **reference's own `_fieldmaps_native`**, compared against the reference's `_displacementmaps`; **sbref, 1 frame, polarity `j`, no polarity term**, inside a crude magnitude-based brain mask built for this experiment (106 268 voxels): p50 = 0.0002 mm, **p95 = 0.006 mm** (threshold < 0.05 mm), p99 = 1.4 mm, max 9.0 mm — the tail being exactly the folded voxels. Displacement range is −10.0…8.5 mm.

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

### 3.5b Phase-encoding POLARITY is load-bearing for `--medic` — B

§3.5 establishes that `wk-apply-warp` ignores the `-` suffix. **`wk-medic` does not.** Running it on the sbref demo with `--phase-encoding-direction j` and then `j-`, everything else identical:

| output | max abs difference | corr(j, j-) |
| --- | --- | --- |
| `_fieldmaps_native` | **0.0000** | +1.000000 |
| `_fieldmaps` | 128.79 Hz | +0.917584 |
| `_displacementmaps` | 16.72 mm | **−0.917584** (median ratio −0.973) |

So the native field map does not depend on polarity — correctly, it is just the weighted regression — but the **inversion direction and the displacement sign do**. The model that reproduces both polarities is §3.4 and §3.3 with a polarity term `s = ±1`:

```text
f_undistorted(y) = f_native( y + s * f_undistorted(y) * TRT )
displacement_mm  = -s * f_undistorted * TRT * pixdim[PE axis]
```

Verified as a **convention check, on the same footing as §3.4 and again NOT end-to-end**: this model, fed the reference's own `_fieldmaps_native`, reproduces the reference's `_displacementmaps` on the **sbref demo (1 frame)** at displacement p95 **0.045 mm** for `j` and **0.024 mm** for `j-`, both inside the 0.05 mm threshold. **Attribution caveat, flagged rather than reconciled:** the support (mask) over which those two percentiles were taken was not recorded at the time, so this row and §3.4's 0.006 mm — nominally the same model on the same data at `j` — differ by 7× for a reason this file cannot document. Both are quoted as measured; neither has been adjusted to agree with the other. See §5.4.

This matters in practice: the supplied three-echo dataset is acquired `j-`. Discarding the sign would apply the correction backwards and roughly double the distortion instead of removing it. Found by external review after the first implementation dropped the suffix — every M0 experiment had used `j` only, so the black-box coverage had a genuine hole.

`-unwarp` still ignores the suffix, and must: by then the sign is already baked into the stored map.

**Support matters, and both earlier figures were right.** §3.4 quotes p95 = 0.006 mm and §3.5b quotes 0.045 mm for what looks like the same measurement. Re-run with the support stated explicitly (200 iterations, reference's own native field, `test/medic_experiments/exp08_inversion.py` conventions):

| polarity | support | p50 | p95 |
| --- | --- | --- | --- |
| `j` | reference mask >= 1 (99 958 vox) | 0.00030 mm | **0.0392 mm** |
| `j` | crude magnitude p60 (106 268 vox) | 0.00020 mm | **0.0059 mm** |
| `j-` | reference mask >= 1 | 0.00024 mm | 0.0213 mm |
| `j-` | crude magnitude p60 | 0.00016 mm | 0.0041 mm |

The 7x spread is entirely the **choice of support**: the reference mask extends further into low-signal regions where the inversion is least well determined. All four are inside the 0.05 mm gate. Quote the support with the number, always.

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

**The measured in-mask test is `>= 1`, not "nonzero", and `--mask` implements exactly that** (`medic.c`, the `--mask` load): a supplied mask voxel is in-mask iff its scaled value is `>= 1.0`. Consequences, all deliberate: a fractional probability map is **not** a mask here (0.9 is *out*; threshold it first, e.g. `niimath p.nii -thr 0.5 -bin m.nii`); NaN fails the comparison and is therefore excluded, which is the safe direction; and a mask with no voxel `>= 1` is a hard error naming the remedy rather than a silently all-zero output set. The reference's `masks.nii` is a 3-level integer image, so `>= 1` is what "in mask" means for it too — passing it verbatim to `--mask` selects levels 1 and 2, which is the region §4.3 shows makes the MCPC offset match exactly.

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

## 4. MCPC-3D-S, ROMEO weights, and the mask -- measured during M4

`--debug` also writes `phase_offset0.nii` (range +-pi, the MCPC-3D-S zero-echo offset) and `phase{0,1}.nii` (per-echo unwrapped phase). These turn out to pin three more conventions exactly.

### 4.1 The MCPC-3D-S formula -- B, exact

```text
hip       = m1*m2 * exp(i*(phi2 - phi1))          (Hermitian inner product of echoes 1,2)
d_uw      = ROMEO_unwrap(angle(hip), mag=|hip|)
offset    = wrap( phi1 - TE1/(TE2-TE1) * d_uw )
```

with **no spatial smoothing**. Against `phase_offset0.nii`, after removing whole-2*pi branch differences in `d_uw`, the residual is p50 2.4e-5 / p95 4.6e-5 rad against a stored quantum of 9.6e-5 rad -- i.e. exact. A smoothed offset could not match at p50 = 0.

Patent position for this stage: see [prior_art.md](../prior_art.md). Summary: US10605885B2 claim 1 requires multi-channel coil data, `TE1:TE2 = n:(n+1)`, and an integer `n`-fold subtraction with no unwrapping; this implementation has none of those (coil-combined input, 16.8:38.56 = 1:2.295, a 0.7721-fold subtraction, and it *does* unwrap). MCPC-3D-S itself is prior art to the 2016 priority date.

### 4.2 ROMEO weight preset -- B

The reference uses **romeo4** weights, at BOTH unwrapping stages -- not ROMEO's own `romeo` default (which resolves to romeo3 when a magnitude is present). Measured as the fraction of in-mask voxels landing on the same 2*pi branch as `phase{0,1}.nii`, with the mask held fixed:

| `-w` | same-branch fraction |
| --- | --- |
| romeo2 | 0.8791 |
| romeo3 / romeo (default) | 0.9301 |
| romeo6 | 0.9817 |
| **romeo4** | **0.9976** |

`--medic` therefore defaults to romeo4; `--weights` overrides it.

### 4.3 One shared mask, and what it costs us -- B

The reference uses a **single mask for both stages** (the MCPC phase-difference unwrap and the multi-echo unwrap), at level >= 1. Supplying that exact mask via `--mask` makes the MCPC offset match **perfectly**:

| HIP-unwrap configuration | offset frac exact (in mask) | p95 |
| --- | --- | --- |
| robustmask(\|hip\|), romeo3 | 0.8236 | 1.4322 rad |
| nomask, romeo4 | 0.8889 | 1.4322 rad |
| dilate(robustmask(\|hip\|), 5), romeo4 | 0.9337 | 1.4322 rad |
| dilate(robustmask(mag1), 4), romeo4 | 0.9318 | 1.4322 rad |
| **reference mask (level >= 1), romeo4** | **1.0000** | **0.0000** |

The 1.4322 rad quantum is exactly `|wrap((TE1/dTE)*2*pi)|` -- these are whole-2*pi branch differences in the phase-difference unwrap, not formula error.

So the mask is the **only** thing separating this implementation from the reference at the MCPC stage, and §3.7's deferral is what now binds. Additional hypotheses tested and rejected since §3.7 (all against `masks.nii`):

| hypothesis | best result |
| --- | --- |
| robustmask(\|hip\| = m1*m2) vs level 2 | dice 0.9754, exact 0.9880 -- **close, but not it** |
| robustmask of sum / mean / geometric-mean / max / RMS magnitude vs level >= 1 | dice <= 0.868 |
| union of per-echo robustmasks vs level >= 1 | dice 0.868 |
| robustmask internal stages (threshold / smooth / fill / final) on \|hip\| vs level >= 1 | dice <= 0.810 |
| binary dilation of level 2 (6- and 26-connected, k = 1..5) vs level >= 1 | dice <= 0.971 |

Level 2 is very nearly `robustmask(|hip|)` (dice 0.975) and level >= 1 is a *looser* 99 958-voxel region that no dilation of level 2 reproduces exactly. `--medic` therefore ships ROMEO's `robustmask` of the first echo's magnitude as the default and exposes `--mask` so exact parity is reachable, and demonstrable, on demand.

### 4.4 `--wrap-limit` -- B, no effect

`wk-medic --wrap-limit` ("turns off some heuristics for phase unwrapping") produces **byte-identical** field maps to the default on the sbref demo (max difference 0.000 Hz). It is not the source of any remaining divergence.

## 5. Measured agreement of this implementation

Everything in this section is **`niimath --medic` end to end** versus the reference's outputs — a different measurement from the convention checks of §3.4/§3.5b, which feed the reference's own native field through one formula. §5.4 attributes every displacement figure in this file to its dataset, stage, mask, polarity and code revision; read it before quoting any number from here.

**sbref demo (1 frame), polarity `j`, errors inside the reference's own mask. Measured BEFORE the mask-gating change** (the unwrapped phase is now zeroed outside the mask before anything reads it), so the two `_displacementmaps` p95 entries are superseded — see the note under the table:

| configuration | output | p50 | p95 | p99 | corr |
| --- | --- | --- | --- | --- | --- |
| shipping default (built-in `robustmask` + romeo4) | `_fieldmaps_native` | 0.0016 Hz | 45.96 Hz | 91.91 Hz | 0.726 |
| | `_displacementmaps` | 0.0004 mm | 2.40 mm | 4.46 mm | 0.771 |
| **`--mask` = reference mask (`masks.nii >= 1`)** | `_fieldmaps_native` | **0.0013 Hz** | **0.0025 Hz** | **0.0027 Hz** | **0.988** |
| | `_fieldmaps` | 0.0056 Hz | 3.47 Hz | 40.30 Hz | 0.958 |
| | `_displacementmaps` | 0.0003 mm | 0.197 mm | 2.29 mm | 0.958 |

**Current values for the `--mask` = reference-mask row after mask gating** (`audit_response.md`, re-audit item 1): `_displacementmaps` p95 **0.059 mm** for `j` (from 0.197 mm) and **0.029 mm** for `j-` (from 0.096 mm) — a 3.3× improvement, and the reason peak RSS rose ~50 MB (§5.3). Only these two percentiles were re-measured; the p50/p99/corr columns and every `_fieldmaps*` row above have **not** been re-run since the change and are therefore pre-gating figures. They are left as measured rather than adjusted.

> **M5 gate status, stated plainly.** The plan's M5 gate (§8) is "reference displacement error below 0.05 mm at the 95th percentile inside the valid mask". End to end, with the reference's own mask supplied, `--medic` is at **p95 0.059 mm for `j` — the gate is NOT met**, by 18 %. `j-` is at 0.029 mm and **does** meet it. The gate *is* met by the isolated inversion/scaling convention checks (§3.4 0.006 mm, §3.5b 0.045/0.024 mm), which is a statement about the formulas, not about the pipeline. With the shipping `robustmask` default the end-to-end figure is 2.40 mm and the gate is missed by ~48×; that difference is the mask, per §3.7/§4.3, and is not being chased.

Read the table carefully: **given the reference's mask, the native field map is exact to 0.0027 Hz at p99** -- the regression, MCPC, unwrapping, rescaling and echo handling are all right. What remains is (a) the mask, and (b) 0.24 % of voxels on a different 2*pi branch, which the inversion then smears along the phase-encoding line (hence `_fieldmaps` p99 40 Hz from `_fieldmaps_native` p99 0.0027 Hz).

`-unwarp` is unaffected by any of this and passes its own gate outright (§5.1).

**170-frame run, shipping default (built-in `robustmask`), polarity `j`, sampled every 17th frame, also PRE-mask-gating:** `_fieldmaps_native` p50 0.37 Hz, p95 26.6 Hz, corr 0.893; `_displacementmaps` p50 0.022 mm, p95 1.37 mm, corr 0.900. The 2.74 GB peak RSS quoted with that run is from an older revision still (before the `uw` buffer was deleted); **§5.3 carries the current memory and timing figures and is authoritative for both.**

### 5.1 M3 gate -- PASSED

Feeding the reference's own `_displacementmaps` to `niimath -unwarp` and comparing against its `echo-{1,2}_part-mag_undistorted.nii.gz`:

| echo | nrmse | p50 | p95 | max | corr | non-finite |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 3.456e-05 | 0.168 | 0.312 | 0.330 | 0.999999999 | 0 |
| 2 | 4.733e-05 | 0.153 | 0.295 | 0.312 | 0.999999999 | 0 |

on data ranging to 41 482, i.e. 8e-6 relative. The residual is the reference map's own uint16 quantization. Ranges match to the last digit (`[-1654.4, 41482.4]` both).

### 5.2 Open item: the low-rank residual

§3.8 established rank-10 truncation from a synthetic series, where the reference's output came back at numerical rank exactly 10. On the **real** 170-frame run it does not. Full-volume singular values of the reference's `_fieldmaps_native`:

```text
95793.2  2830.7  2321.5  2113.6  1526.5  1093.4  954.8  851.1  806.1  770.6 | 383.3  347.2  279.4  243.3 ...
```

There is a clean factor-2 drop at index 10 (770.6 -> 383.3), so a rank-10 truncation *was* applied -- but a broadband residual survives it, three orders of magnitude above the uint16 quantization floor (0.55). This implementation's output is strictly rank 10 (sv[10..] = 2.4e-3).

Hypotheses not yet discriminated: truncation applied per temporal-correlation group rather than globally; a residual add-back (Eq. 9 read as a correction rather than a replacement); or truncation applied before a later full-rank stage. **Deliberately not guessed.** `--rank 0` disables the filter for anyone who wants the raw regression.

### 5.3 Performance vs the reference (M12)

Same workload as `demo/run170.sh` (76x76x46 x 170 frames, 2 echoes, magnitude + phase), shipping default mask, Apple Silicon (10P+4E, 48 GB). Reproduce with `test/medic_experiments/bench_threads.py` (the like-for-like table, 1 and 8 threads, both tools uncompressed) and `test/medic_experiments/bench170.py` (the gzipped-output and apply-stage tables, 8 threads). **These four tables are the authoritative timing and memory figures for MEDIC; `AGENTS.md`, `README.md` and `medic_plan.md` quote them and must not contradict them.**

**Estimate stage, like for like** — both writing UNCOMPRESSED `.nii`, so this compares compute rather than gzip (`test/medic_experiments/bench_threads.py`):

| threads | tool | wall | CPU | parallelism | peak RAM |
| --- | --- | --- | --- | --- | --- |
| 1 | `wk-medic` | 44.37 s | 53.29 s | 1.2x | 3.40 GB |
| 1 | **`niimath --medic`** | **14.52 s** | **13.96 s** | 1.0x | **1.70 GB** |
| 8 | `wk-medic` | 15.23 s | 65.51 s | 4.3x | 3.40 GB |
| 8 | **`niimath --medic`** | **3.87 s** | **15.16 s** | 3.9x | **1.89 GB** |

**3.06x faster single-threaded, 3.94x faster at 8 threads**, using ~4x less CPU and ~1.8x less RAM — while writing float32 (172 MB/series) against the reference's uint16 (86 MB/series). Thread scaling 1→8 is 3.75x for niimath and 2.91x for the reference. niimath now scales BETTER as well as running faster (3.75x versus 2.91x). Note also that the reference spends 53 s of CPU to do single-threaded what niimath does in 14.5 s: its parallel gain is largely recovering its own overhead.

These are the current revision. Two changes moved them since the previous round, in opposite directions and both deliberately:

- **Peak RAM fell 2.04 → 1.70 GB (1 thread) and 2.24 → 1.89 GB (8 threads)** once geometry validation moved to headers alone and each echo pair is loaded, repacked and freed in turn. Previously all `2E` input payloads were still resident when the work arrays were allocated, so the true peak was ~`(4E+3)` series while the banner reported `(2E+3)`.
- **Wall time fell despite that round ADDING correctness work** (4.15 → 3.42 s at 8 threads as measured then; the table above is the current figure) (mask gating, per-voxel temporal validity counting, fold detection), because the group-mean accumulation is hoisted when every frame falls in one temporal group — exact, the same values summed in the same order, removing ~7.7e9 float adds on this dataset.

The mask retention that costs a little memory bought a 3.3x improvement in end-to-end displacement agreement on the sbref demo with the reference mask supplied (p95 0.197 → 0.059 mm for `j`, 0.096 → 0.029 mm for `j-`; §5, §5.4).

**Treat the RATIO as the robust quantity.** Runs taken while other jobs competed for cores showed both tools ~1.5x slower with the ratio preserved (3.1x single-threaded, 4.1x at 8 threads). Re-measure on an idle machine before quoting absolute seconds.

**Estimate stage, gzipped output** (`wk-medic` vs `niimath --medic`):

| tool | wall | CPU | parallelism | peak RAM |
| --- | --- | --- | --- | --- |
| `wk-medic` | 15.63 s | 65.07 s | 4.2x | 3.40 GB |
| `niimath --medic` (gz out) | **10.64 s** | **19.95 s** | 1.9x | **2.35 GB** |
| `niimath --medic` (`--gz 0`) | **4.29 s** | 13.39 s | 3.1x | 2.19 GB |

(Those three rows are **historical** — separate runs, taken before the per-echo loading, the removal of `md_write()`'s output copy and the OpenMP level change. Every figure in them is superseded by the like-for-like table above, which is the authoritative one; the `--gz 0` peak of 2.19 GB there reads **1.89 GB** now. The gzipped peak has not been re-measured since those changes, so no current gzipped figure is quoted anywhere — use the uncompressed numbers.)

**Apply stage** (one echo, 170 frames):

| tool | wall | CPU | parallelism | peak RAM |
| --- | --- | --- | --- | --- |
| `wk-apply-warp` | 50.72 s | 485.62 s | 9.6x | 1.14 GB |
| `niimath -unwarp` | **4.68 s** | **17.47 s** | 3.7x | **0.59 GB** |

So: 1.5x faster and 3.3x less CPU on the estimate, 10.8x faster and 28x less CPU on the apply, at 1.3-1.9x less RAM — while writing float32 (172 MB/series) against the reference's uint16 (86 MB/series). The reference parallelises harder (4.2x / 9.6x) but spends far more total CPU to get there.

Two build gotchas that invalidate this table if ignored:

- **Verify OpenMP is actually linked.** `./src/niimath <img> -p 8 -s 1 out.nii` must print `Using 8 threads`. A stale object from a `make OMP=0 / ROMEO=0 / MEDIC=0` build silently yields a serial binary; a first run of this benchmark measured 1.0x parallelism for exactly that reason and had to be discarded.
- **Build against zlib-ng, not system zlib** (`make -C src ZLIBNG_ROOT=...`, or the CMake release path, which defaults to it). Output gzip is the dominant serial tail: system zlib gives 16.21 s where zlib-ng gives 10.64 s. `--gz 0` isolates it at 4.29 s, and is the like-for-like comparison since the reference writes uncompressed `.nii`.

### Is everything held in RAM?

**niimath: yes, by design, and that is a decided position rather than an omission.** A 4D `.nii.gz` cannot be seeked, so every tool — including the reference — reads and writes whole volumes in RAM anyway; a streaming layer would buy nothing for the dominant gzip case. The requirement is instead to be fast and honest about the cost.

The **work arrays** are `phase + mag + fields + fu + disp` = `n3 * T * (2*echoes + 3) * 4` bytes — 1.18 GiB on this run, and that is the figure `--medic` prints at startup. It grows linearly with frames x echoes, so a 5-echo/600-frame run at this resolution needs roughly 10 GiB of work arrays — fine natively, and **impossible in wasm32**, whose 4 GiB address space (and `-DFORCE_INT32_MAX` on every wasm target) is why `--medic` is documented as native-scale-only while `-unwarp` is browser-friendly.

**The banner is the work-array budget, not the peak — and the ordering that keeps the two close is load-bearing.** In the current `medic.c` every input is validated from its **header alone**, so no payload is resident when the five work arrays are allocated, and the repack loop then loads, rescales and frees **one echo pair at a time**. The transient input overshoot is one echo pair (2 series) regardless of echo count; `md_write()` now lends each resident output directly to the synchronous writer instead of copying another complete series. An earlier revision read all `2*echoes` payloads during validation and freed them only during the repack, i.e. after the work allocation, so its true peak was `(4*echoes + 3)` series against a `(2*echoes + 3)` banner.

An earlier revision carried a separate unwrapped-phase series; it was removed during the audit (phase is unwrapped in place), taking the gzipped peak from 2.53 GB to 2.35 GB.

**The reference: also yes, and more.** Peak footprint is 3.40 GB against inputs that are 0.67 GiB as float32 and 1.35 GiB as float64, so it is holding every series resident (float64, on the evidence of the ratio) plus unwrapped phase, intermediates and Python overhead. Neither tool streams; niimath simply keeps a smaller resident set by working in float32 throughout.

### 5.4 Where every displacement figure in this file comes from

Six different displacement percentiles appear above and they are **not** measurements of the same thing. Quote a figure only with its row.

| p95 | § | what was compared | dataset | mask / support | polarity | revision |
| --- | --- | --- | --- | --- | --- | --- |
| **0.006 mm** | §3.4 | *Convention check.* Converged fixed-point inversion + §3.3 Hz→mm applied to the **reference's own** `_fieldmaps_native`, vs the reference's `_displacementmaps`. No `--medic` pipeline involved. | sbref, 1 frame | crude magnitude-based brain mask, 106 268 voxels | `j` (no polarity term) | M0, pre-implementation |
| **0.045 mm** (`j`), **0.024 mm** (`j-`) | §3.5b | *Convention check*, same as above with the polarity term `s = ±1` added and run once per polarity. | sbref, 1 frame | **not recorded** — see the caveat in §3.5b | `j` and `j-` | audit round 1 |
| **0.197 mm** | §5 | *End to end.* `niimath --medic` `_displacementmaps` vs the reference's. | sbref, 1 frame | `--mask` = reference `masks.nii >= 1` | `j` | before mask gating |
| **0.059 mm** (`j`), **0.029 mm** (`j-`) | §5, §5.3 | *End to end*, same configuration, after the unwrapped phase is zeroed outside the mask. **This is the current end-to-end number, and 0.059 mm misses the 0.05 mm M5 gate.** | sbref, 1 frame | `--mask` = reference mask | `j` and `j-` | current |
| **2.40 mm** | §5 | *End to end*, shipping default mask. | sbref, 1 frame | built-in `robustmask` of echo 1 | `j` | before mask gating |
| **1.37 mm** | §5 | *End to end*, shipping default mask, every 17th frame. | 170 frames | built-in `robustmask` per frame | `j` | before mask gating |

Two things this table deliberately does **not** do. It does not reconcile §3.4's 0.006 mm with §3.5b's 0.045 mm — same model, same data, same polarity, 7× apart, and the missing support makes the difference undiagnosable from the record; re-running §3.5b with a stated mask is the only honest fix. And it does not back-fill the pre-gating rows with post-gating values: only the two `--mask` `_displacementmaps` p95 entries were re-measured after that change.

## 5.5 FP policy for medic.c — MEASURED, no strict FP needed

Three audits flagged that `md_rescale_phase()` and `md_mcpc3ds()` run under `-ffast-math` while feeding ROMEO's bin-quantised 8-bit edge weights — the exact mechanism that forces `romeo.c` strict. Settled the way ROMEO's own policy was settled, with `test/medic_experiments/fp_policy_medic.py`.

Two binaries differing **only** in medic.c's FP flags (65 fused multiply-adds versus 0; the fast object is bit-identical to the shipped LTO build, so this measures what ships).

| measurement | result |
| --- | --- |
| configurations | 32 (sbref × 8 option sets, 170-frame BOLD × 4, 12 adversarial synthetics) |
| voxels compared | 508 051 136 |
| **whole-2π branch differences in unwrapped phase** | **0** |
| **mask differences** | **0** |
| largest unwrapped-phase difference | 3.8e-6 rad (1–4 float32 ULP) |
| largest `_fieldmaps_native` difference | 3.05e-5 Hz |

Adversarial cases included ±π wrap boundaries, a non-representable 2π/3 rescale slope, a gradient sweep driving weights through all 255 `rescale()` bin edges, 1-ULP-tied and 1e±30 magnitudes, and raw-Siemens-scaled versions so `md_rescale_phase()` actually runs.

**Mechanism isolated:** every `offset=none` configuration is **bit-identical**, including raw-scale synthetics where the rescale slope and intercept are both non-terminating in binary and the FMA is live. `md_rescale_phase()` did not move a single float32 value anywhere in the corpus. All divergence comes from `md_mcpc3ds()`, at 1–4 ULP on ~1e-3 of voxels.

**Saturating probe** (the answer to "how hard did you try"): perturbing the phase handed to ROMEO by ±1 float32 ULP at **100 %** of voxels — one binary run twice, on smooth, residue-laden and fully inconsistent fields — still gave 0 branch flips and 0 mask changes. Measured in `romeo_plan.md`'s own units via `-romeo-dump`:

| case | weight bytes differing | dropped to bin 0 (edge deleted) | 2π branch diffs |
| --- | --- | --- | --- |
| smooth | 1 / 221 184 | **0** | 0 |
| fully inconsistent | 4 / 221 184 | **0** | 0 |
| *(`romeo.c` built `-ffast-math`, for scale)* | *360 / 797 088* | *yes — edges deleted* | *66 voxels* |

The medic path is 25–90× quieter at the byte level and **never deletes an edge**, which is the specific mechanism that breaks `romeo.c`.

**One separate finding, not an FP defect.** On the 170-frame run the post-inversion `_fieldmaps`/`_displacementmaps` differ by up to 29.2 Hz / 1.66 mm — at **21 voxels of 45 168 320**, 20 of them in folded PE columns, while the pre-inversion native field differs by ≤3.05e-5 Hz. The saturating probe reproduces the same amplification (130 Hz / 7.9 mm) from a pure ±1 ULP perturbation with a single binary. This is `md_invert()` fixed-point conditioning where the inverse is genuinely multi-valued — the run already warns about it — so strict FP would pick a *different* arbitrary branch, not a well-defined one. The lever is `MD_INVERT_ITERS` and fold handling, not compiler flags.

## 6. What still needs porting

> **M0 snapshot, since delivered.** MCPC-3D-S landed in M4 (`md_mcpc3ds()`, ordinary-FP in `medic.c` rather than beside the strict-FP ROMEO unit) and reproduces the reference's own `phase_offset` exactly under a shared mask — see §4. Retained because the measurements below are the record that motivated it.

`--debug` also writes `phase_offset0.nii` (range ±π, the MCPC-3D-S zero-echo offset) and `phase{0,1}.nii` (per-echo unwrapped phase). Comparing niimath's current `-romeo` against `phase{0,1}.nii` shows the expected large disagreement — median 4.40 rad at echo 1, with 49 441 of 64 877 in-mask voxels off by a whole 2π — because niimath does **not** yet remove the phase offset before unwrapping. Once offsets are removed the unwrapped phases are near-perfectly linear in TE: `median(phi_2/phi_1) = 2.295230` versus `TE_2/TE_1 = 2.295238`.

**MCPC-3D-S is the one genuinely new numeric kernel** and is the M4 deliverable. Port the monopolar path only, from the pinned MIT MriResearchTools source, beside the strict-FP ROMEO code (plan §4.2).

## 7. M0 gate

Every implementation-sensitive convention is now measured and recorded, or explicitly deferred behind a documented decision:

| convention | status |
| --- | --- |
| phase scaling | measured, already implemented (§3.1) |
| weighted regression | measured exact (§3.2) |
| Hz→mm sign and units | measured exact (§3.3) |
| inversion | measured; converged fixed point adopted, residual quantified, convention check clears the 0.05 mm threshold 8× (§3.4 — not the end-to-end figure, see §5.4) |
| PE-axis semantics, displacement sign | measured over 30 configurations (§3.5) |
| interpolation kernel, fill, Jacobian | measured exact (§3.6) |
| masking | measured *behaviour*; construction not reproduced — deliberate, gate on common mask (§3.7) |
| low-rank rank and centering | measured (§3.8) |
| temporal correction | behaviour measured; parameters from the paper (§3.9) |
| noise frames | measured (§3.10) |
| border filter | not exposed; not implemented (§3.11) |
| output headers/datatype | measured; niimath deliberately writes float32 (§2) |

No code below rests on a guessed sign, scaling, interpolation, or inversion behaviour. The two open items (§3.7 mask construction, §3.9 grouping parameters) are named, bounded, and carry a validation plan.
