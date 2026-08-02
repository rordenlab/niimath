# `-skullstrip` — restart document

**This replaces the original execution plan.** The native implementation is complete and the deformation is at parity with AFNI. Historical experiments live in `skullstrip_bench/test/skullstrip_reference_manifest.md`; that repository is currently local-only, so the evidence is not independently available from a clean niimath checkout.

Read this file, then the manifest, then the `skullstrip.c/.h` bullet in `AGENTS.md`. Do not start from memory.

---

## Where things stand

An AFNI-style deformable-surface skull stripper: no template, no mask, no network. **OFF by default** (`SKULLSTRIP=1 make`, `cmake -DENABLE_SKULLSTRIP=ON`), **64-bit native only**, absent from every released binary.

The faithful kernel's original five-image parity set, measured with the harness built at the shipped binary's flags:

| image | Dice vs AFNI `-no_use_edge` | Dice vs external reference masks | AFNI vs external reference masks |
| --- | --- | --- | --- |
| T1w | 0.9788 | 0.9522 | 0.9434 |
| T1w1mm | 0.9825 | 0.9652 | 0.9596 |
| T1w2mm | 0.9688 | 0.9507 | 0.9530 |
| T1w_ARC2017 | 0.9639 | 0.8241 | 0.8055 |
| T1w_MICCAI2017 | 0.9543 | 0.8625 | 0.8304 |
| **mean** | **0.9697** | **0.9109** | 0.8984 |

The current command compiles the deformation loop twice from `skullstrip_kernel.h`. `-faithful` reproduces the pre-optimisation kernel; the default removes code-generation costs and parallelises independent nodes. On the newer nine-image set, faithful/default mean Dice against AFNI is 0.9407/0.9409. Best-of-three M4 Pro timings are 37.4 s for the original baseline, 17.4 s faithful, and 9.6 s default. Default output is byte-identical at `-p 1/3/8` for the measured corpus.

`make test`, CTest, serial UBSan, malloc diagnostics, and `leaks` are clean. The local evidence repository is at `62e4d9b`.

---

## Measurement traps — read before quoting any number

Each of these produced a wrong published figure at least once.

1. **`3dSkullStrip -mask_vol` is NOT a binary mask.** It writes a graded code in which **1 means OUTSIDE the surface**. AFNI's own rule is `isin >= 3`. Thresholding at `> 0` inflates the reference by 211 mL on T1w and wrecks every Dice against it — same mask scores 0.9227 at `>0` and 0.9812 at `>=3`.
2. **Build the harness with the SHIPPED binary's flags.** A `-O2` driver and the `-O3 -ffast-math` build disagree by 36,596 voxels on T1w (Dice 0.991). Wrong flags measure a different program.
3. **Never derive a mask as `output > min`.** `-skullstrip` blanks to the image minimum and 18,969 in-brain voxels on T1w sit exactly there. Compare masks, not stripped volumes.
4. **Do not require byte identity across builds or between fast and faithful.** The stage loop's integer convergence threshold amplifies sub-ULP code-generation changes. CMake vs Makefile differ by 26,707 voxels, and merely adding a warning block moved three of five masks. Use the faithful specialization for exact within-build regression work; assess intentional fast-kernel changes with the nine-image aggregate and per-image Dice gates in `AGENTS.md`.
5. **A synthetic fixture needs a plausible physical SIZE.** Normalisation measures the head in millimetres over a fixed 167x212x175 grid, so a 48³ phantom at 1 mm is correctly rejected as degenerate. Use 4 mm voxels to make the same array a ~190 mm head.

## Verification traps — how *I* got these wrong, not just the code

6. **`cmp` on two `.nii.gz` files is not a content comparison.** gzip embeds an mtime, so two byte-identical images written a second apart compare as different. Always `-gz 0` when byte-comparing, or compare the decompressed arrays. This produced a false "the verbose flag changes the segmentation" alarm.
7. **Reuse of scratch filenames across a loop produces phantom diffs.** `rm -f` the outputs and check the exit code of every run before comparing. A stale file from the previous iteration reads as a difference in the current one. Combined with trap 6 this nearly caused a correct change to be reverted.
8. **When a comparison surprises you, bisect before concluding.** The fastest disproof here was to make the suspect code path unconditional and re-test; that took one rebuild and settled it.

## Code traps

- **A diagnostic must not gate WORK, only OUTPUT — but the current `ss_verbose() ? &xs : NULL` gate on the exhaustive self-intersection scan IS safe and was verified.** All five images plus the synthetic fixture produce byte-identical output with and without `SKULLSTRIP_VERBOSE`, and `release_smoke.py` asserts it on every run. The scan is genuinely pure and its result reaches only a log line. Keep that assertion: it is what distinguishes this from the bug below.
- **Never put a function call inside `SSV(...)`.** The macro guards its body on `ss_verbose()`, so arguments are not evaluated in a normal run. Two pipeline stages were once written that way and silently did not run unless `SKULLSTRIP_VERBOSE` was set — a diagnostic that changed the segmentation, and the whole test suite stayed green through it. `release_smoke.py` now byte-compares verbose vs non-verbose output; keep that check.
- **`ss_afni_self_intersect` deliberately reproduces an AFNI quirk** (it only sees edges pointing into the +,+,+ octant, so it detects roughly one fold in eight). It decides the RETRY; `ss_mesh_self_intersections` is the correct test and is for reporting only. Do not "fix" it.
- **An OOM sentinel must not collide with a valid result**, and one sentinel must not be read two ways by two callers. Both mistakes shipped here.
- **`ss_mesh_normals` must not allocate** — it is `void` and cannot report failure. Keep the divisor from `nbr_off`.
- **Keep `skullstrip.c` self-contained.** It is linked alone by the mesh selftest, the CMake selftest and the benchmark drivers; do not add a `core.c` dependency for three lines.
- **No `//` comments inside the line-continued cell-range macros** — they swallow the continuations.
- The full list, with the measurement behind each, is the `skullstrip.c/.h` bullet in `AGENTS.md`.

---

## Uncommitted state — READ BEFORE COMMITTING

The working tree contains **two unrelated bodies of work**. Do not commit them together.

1. **`-skullstrip`** — everything in this document. Reviewed, validated, ready.
2. **A separate GPL/bandpass change**, which predates this effort and was never part of it. It has since gone all the way: `-bandpass` is **RETIRED** and `bw.c`/`bw.h` are deleted from the project *and* from the `src/GPL` submodule — not moved, not `GPL=1`-only, gone. Exstrom's LGPL-3 grant was the only LGPL-3 component in niimath and the sole reason a `GPL=1` binary resolved to GPL-3; with it gone the copyleft payload is SPM's `spm_coreg` alone and such a binary is **GPL-2-or-later**. `-bptf`/`-bptfm`, which are different (BSD) code, are unaffected. This is a user-visible behaviour change with its own release implications — the npm `bandpass()` method disappears when the generated API is next rebuilt from the help text.

`js/bun.lock` and `js/package.json` are also modified and unrelated to both.

## Release gates — four separate ones, deliberately

Earlier versions mixed implementation, distribution, and WASM readiness. They are independent gates.

| gate | status | what closes it |
| --- | --- | --- |
| **Native beta** — off by default, usable by anyone who opts in | NOT met for public distribution | implementation and tests are ready; publish or vendor the cited provenance evidence |
| **Native release** — on by default in shipped binaries | NOT met | first close provenance, then make the explicit product decision on enabling this optional example |
| **Provenance** — independently auditable | NOT met | `thd_coords.c` was resolved by clean room, but `skullstrip_bench` has no remote and niimath does not vendor the evidence |
| **WASM** — the plan's original stated product | NOT met, NOT started | `SKULLSTRIP_WASM` is an unreferenced placeholder; every wasm/tiny/nano build rejects the request explicitly |

Both `src/Makefile` and `src/CMakeLists.txt` already state the current rationale ("OFF by default pending the native-release gate ... NOT because of accuracy -- the two weakest images now sit inside the agreement band"). The older "two validation images miss the agreement band" wording is gone from both; do not re-add it.

## Next steps, in priority order

### 1. ~~Stored-datatype heuristic~~ — RESOLVED

The pre-promotion NIfTI datatype is now plumbed from `coreFLT.c` (`ihdr.datatype`) through `skullstrip_run()` and `ss_normalize()` into `ss_to_rai_short()`, which branches on the STORED type as AFNI does instead of guessing from the values. The value inspection survives only as a `DT_NONE` fallback for standalone harnesses.

Demonstrated: the same voxel data stored as int16 and as float32 now takes **different** normalisation paths (42,812 voxels, Dice 0.9895) — which is AFNI's behaviour. Under the old heuristic they were identical, because it tested a property of the values and the values are the same. All five validation images are stored integer, so the corpus could never have caught this; that is why it needed a paired fixture rather than another benchmark run.

### 2. Publish or vendor `skullstrip_bench` — OPEN

Committed locally (`62e4d9b`), with no remote. The `.gitignore` excludes image volumes. Publish the repository and link it, or vendor its redistributable manifest and clean-room kit into niimath. Until then, a clean checkout cannot inspect the derivation or reproduce the licence conclusion.

### 3. Two real but non-blocking accuracy gaps

- `T1w2mm` fails the harness's `contains_head_com` check — **but AFNI's own mask fails it too** when binarised correctly. The gate is finding a genuine inferior-extent limitation both strippers share, not our defect.
- ARC2017 and MICCAI2017 sit at Dice 0.82-0.86 against external reference masks where the other three are 0.95. AFNI is also low there (0.81, 0.83), so parity is fine — but if you want absolute quality rather than parity, that is where it is.

### 4. WASM — the stated product, still not built

The original plan called WebAssembly "the reason the project exists". It is **not implemented**. `SKULLSTRIP_WASM` is an unreferenced placeholder, and every wasm/tiny/nano build rejects an explicit request rather than silently dropping it. Treat this as a separate product milestone; it does not block the optional native example.

### 5. Deferred deliberately

- Unifying the two spatial-grid builders (~120 lines → ~45). Real duplication, but one of them decides the retry and therefore the whole result. Preserve the arithmetic verbatim and use the five-image byte comparison as the check.
- Milestone 4, the edge-assisted AFNI default: measured worth only 0.004-0.06 Dice, and it is the GPL-3.0-adjacent path. Probably never.
- The dedicated skullstrip workflow runs enabled Make and CMake smoke paths on Linux, CTest, and an enabled MSVC compile/smoke job.

---

## Licence position — implementation resolved, evidence not yet published

Normalisation, intensity prep, deformation and touchup are **adapted from public-domain AFNI** (`thd_brainormalize.c`, `thd_automask.c`, `SUMA_BrainWrap.c`). The basis is affirmative rather than the absence of a notice: AFNI's `LICENSE.txt` declares the tree a US Government Work and states that "contributions without explicit licensing will be assumed to be entered into the public domain", and `README.copyright` dates the rule to after 15 Jan 2001 — the adapted files' first commits are 2001-2004 by NIH authors.

**The one carve-out is `SUMA_3dedge3`**, which wraps Malandain's GPL-3.0 code. Our `-no_use_edge` contract never reaches it. Do not read it.

**`thd_coords.c` was the one real problem and it is now resolved.** An MCW file from 1998, outside both rules, reached indirectly by following a call out of the vetted `SUMA_BrainWrap.c`. (It was GPL-2 when this was clean-roomed; MCW relicensed its 1994-2000 AFNI code to CC BY 4.0 on 2026-05-12, which removes the copyleft bar but adds attribution and change-notice duties -- so the clean-room result is still the preferable outcome. See the AFNI licensing section in AGENTS.md.) The index conversion derived from it has been **replaced by a clean-room reimplementation** — see `skullstrip_bench/clean_room/`. An implementer who had read neither codebase reproduced all 19,139 rows of a measured behaviour table, and the resulting `ss_world_to_index()` gives byte-identical masks. Caveat recorded in the manifest: that implementer was an LLM, so "has not read the original" is weaker than it would be for a person. The kit is designed to be handed to a human unchanged if you want the stronger guarantee.

**Lesson worth keeping: vet callees, not just files.** The provenance table covered the files we meant to adapt and missed one their code called.

---

## Picking the work back up

```sh
cd /Users/chris/src/niimath/src && SKULLSTRIP=1 make test      # ~30 s, must be green
```

The measurement harness is in `skullstrip_bench/scripts/`: `ss_strip.c` (dev driver), `run.sh`, `score.py`, `maskradii.py`, `cmp_node.py`, `traj.py`. Build the driver with `-O3 -ffast-math -fno-finite-math-only` — see trap 2.

For per-node parity against AFNI: build with `-DSS_DEV_HOOKS=1` to re-enable `SS_NODE_DBG` (the developer hooks are compile-gated so no environment variable can change a shipped result), and drive AFNI with `-node_dbg N -debug 3`. **`-debug 3` is required** — the per-iteration write is gated on `LocalHead`, which only `> 2` sets. Column layout and the frame mapping are in the manifest.
