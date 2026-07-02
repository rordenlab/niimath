# Audit response — release hardening after Windows/MSVC merge

This supersedes the previous audit response. I reviewed the current working-tree changes plus the recent Windows/MSVC merge, with emphasis on allocator ownership, fail-closed behavior, WASM/browser packaging, and duplicate build paths.

## Changes made

1. Fixed a JS worker output edge case. `workerImpl.ts` now reads `outName` and then `outName + ".gz"` when the requested output was not already gzip-named, matching niimath/FSL gzip defaults and the direct WASM test helper. It also reports the actual output name and unlinks both possible output spellings so failed or gzip-renamed runs do not leave stale MEMFS files.

2. Restored GPL TypeScript packaging. The earlier `tsconfig` exclusion kept BSD-only type checks green but stopped `dist/index-gpl.d.ts` from being emitted for a GPL build even though `package.json` exports it. I removed that exclusion and added tiny source declarations for the generated Emscripten modules (`js/src/niimath.d.ts`, `js/src/niimath-gpl.d.ts`) so `tsc` works even when ignored generated JS artifacts are absent.

3. Made BSD-only JS builds non-dangling. `esbuild.config.ts` now writes a clear GPL-unavailable stub for `dist/index-gpl.js`/`dist/niimath-gpl.js` when the GPL WASM artifacts are absent, removes any stale `dist/niimath-gpl.wasm`, and removes stale `js/corresponding-source/`. The test helper now treats GPL as built only when both `dist/niimath-gpl.js` and `dist/niimath-gpl.wasm` exist, so a stub does not accidentally run GPL tests.

4. Fixed the browser package build. A real Playwright browser smoke showed `dist/index.js` loading and then requesting `/dist/core` (404): `tsc --project tsconfig.json` was overwriting esbuild's bundled browser output with extensionless native ESM. `esbuild.config.ts` now cleans `dist/` at the start and runs `tsc --emitDeclarationOnly`, so esbuild owns JS output while TypeScript still emits declarations.

5. Fixed Emscripten 6 resizable-memory browser output. The generated module used `wasmMemory.toResizableBuffer()`, then Emscripten's own UTF-8 helper passed HEAP subarrays backed by that resizable buffer to `TextDecoder.decode()`, which current Chromium rejects. `scripts/pre-build.ts` now rewrites `getMemoryBuffer()` to return `wasmMemory.buffer`; Emscripten's existing `updateMemoryViews()` still refreshes views after memory growth, but browser `TextDecoder` sees normal ArrayBuffers.

6. Finished fail-closed allocation hardening in `allineate.c`. `al_scalar_ransetup()` now returns an error on allocation failure and checks the coarse-start buffers (`kpar`, `all_wpar`, `all_isrand`, `tvals`, `tpars`, `tisrand`) before use. The refinement and fine candidate `cand_wpar` arrays are also checked before OpenMP workers can dereference them.

7. Fixed an OOM cleanup leak introduced by early `goto al_cleanup` paths during the downsampled coarse pass. If setup/startup allocation fails before the normal full-resolution restore block runs, cleanup now frees the downsampled source buffer, frees any smoothed downsample buffer, restores source geometry fields, and restores `stup.ajim_orig` so the source automask backup is freed by the standard cleanup path.

8. Removed the stale `uf_isort` helper and updated the `unifize.c` comment to match the quickselect-based trimmed mean implementation. This avoids future confusion about using insertion sort on the large per-voxel neighborhood.

9. Hardened the PyPI/AppVeyor release paths after the reported PyPI `v1.0.20260315 Clang21.0.0` short-read case. I reproduced the user's command against current source and an installed local wheel using `/Users/chris/src/tmp/brainchop/strokeSubject.nii.gz`; current builds do not emit `++ WARNING: read ...` and produce a valid 256^3 uint8 NIfTI (`16782752` bytes including the 5536-byte header/extension area). The old warning expected `67108864` bytes, i.e. a 256^3 float payload, so the new release smoke explicitly catches any future header/data-size mismatch.

10. Added `.github/scripts/release_smoke.py`, a stdlib-only packaged-binary smoke test that synthesizes NIfTI fixtures, checks the exact `-conform -gz 0 ... -odt char` class, verifies gzip read/write, verifies zstd read/write when expected, checks key feature dispatch in the packaged help, and confirms BSD wheels reject `-spm_coreg` with the GPL-module message.

11. Changed PyPI wheel CI from a help-only smoke (`niimath`) to the release smoke and removed the forced `ENABLE_ZSTD=OFF`. The GitHub Actions wheel matrix now installs zstd/libomp where needed, passes `ENABLE_ZSTD=ON`, enables macOS OpenMP wheels with `OPENMP_XCODE=ON`, and uses vcpkg's static zstd on Windows. I also checked cibuildwheel's current platform behavior: Linux wheels can include manylinux and musllinux containers, so the zstd install step now handles RPM (`yum`/`dnf`/`microdnf`), Debian (`apt-get`), and Alpine (`apk`) package managers instead of assuming a manylinux-only RPM image.

12. Fixed zstd portability for CMake/scikit-build packaging. `src/CMakeLists.txt` now prefers static `libzstd.a` when available, recognizes Windows static names (`zstd_static`/`libzstd_static`), and honors an explicit `ZSTD_ROOT`; `SuperBuild/SuperBuild.cmake` forwards `ZSTD_ROOT` and `CMAKE_PREFIX_PATH` into the nested `src` configure so CI-provided zstd installations are visible to the executable build. A first local OpenMP wheel linked to `/opt/homebrew/opt/zstd/lib/libzstd.1.dylib`; after this fix, `otool -L` on the wheel's `niimath` binary reports only `/usr/lib/libSystem.B.dylib`.

13. Wired the same release smoke into AppVeyor GitHub-release ZIP builds. Linux/test jobs install `libzstd-dev`; Windows builds a static zstd from upstream source and passes it through CMake; the macOS universal artifact continues using its existing universal static zstd path and now smokes the final lipo'd binary before upload.

## Residual notes

The broader project still has known unchecked allocations outside this touched path, especially in older mesh code. I did not sweep those here because they are already tracked as a project-wide issue and are outside the release-facing Windows/MSVC/WASM/registration changes under review.

The optional GPL SPM module remains intentionally absent from the PyPI `niimath` wheel because that package is declared BSD-licensed and exposes a single `niimath` console script. Shipping GPL code there would make the wheel a GPL-2 combined work and require changing PyPI metadata/license posture. The npm package already handles this with a separate `@niivue/niimath/gpl` export and corresponding-source bundle; native GPL builds remain covered by `GPL=1 make`, `cmake -DENABLE_GPL=ON`, and `gpl-build.yml`.

## Verification

- `make` in `src` passed.
- `make GPL=1` in `src` passed.
- `cmake -S src -B /private/tmp/niimath-cmake-src-bsd-audit` and `cmake --build /private/tmp/niimath-cmake-src-bsd-audit` passed.
- `cmake -S src -B /private/tmp/niimath-cmake-src-gpl-audit -DENABLE_GPL=ON` and `cmake --build /private/tmp/niimath-cmake-src-gpl-audit` passed.
- CMake BSD `-spm_coreg` smoke failed with the expected GPL-required message; CMake GPL `-spm_coreg` smoke passed.
- `bun x tsc --project tsconfig.json --noEmit` passed with generated Emscripten JS present and also with `js/src/niimath.js`/`js/src/niimath-gpl.js` temporarily absent.
- `bun run makeWasm`, `bun run scripts/pre-build.ts -i src/niimath.js -o src/niimath.js`, and `bun run makeWasmGpl` passed.
- GPL-enabled `bun run esbuild.config.ts` passed and emitted `dist/index-gpl.d.ts`.
- GPL-enabled `bun test ./tests` passed: 8 pass, 0 fail.
- Simulated BSD-only `bun run esbuild.config.ts` passed after temporarily hiding GPL source artifacts; it removed stale GPL WASM/source output, emitted GPL stubs, and `bun test ./tests` passed with 5 pass, 4 skip.
- Mocked worker harness passed: requested `out.nii`, found `out.nii.gz`, returned the gzip bytes, and left 0 staged MEMFS files.
- Playwright/Chromium browser smoke passed through the actual `dist/` package and module workers: BSD default gzip fallback, BSD explicit `.nii.gz`, BSD `-allineate`, BSD rejection of `-spm_coreg`, BSD `resliceNN` + `mulImage` file operands, and GPL `-spm_coreg`.
- `leaks --atExit -- ./niimath ... -spm_coreg ...` reported 0 leaks.
- Current native GPL binary: `./src/niimath /Users/chris/src/tmp/brainchop/strokeSubject.nii.gz -conform -gz 0 /private/tmp/brainchop-conform-current.nii.gz -odt char` passed with no short-read warning; resulting `/private/tmp/brainchop-conform-current.nii` is datatype `2`, `nbyper=1`, `nvox=16777216`, size `16782752`.
- New release smoke passed against `./src/niimath --expect-zstd`.
- Top-level SuperBuild path passed: `cmake -S . -B /private/tmp/niimath-superbuild-release-audit -DENABLE_ZSTD=ON -DOPENMP_XCODE=OFF`, `cmake --build /private/tmp/niimath-superbuild-release-audit`, and `python3 .github/scripts/release_smoke.py /private/tmp/niimath-superbuild-release-audit/bin/niimath --expect-bsd --expect-zstd`.
- Local scikit-build wheel path passed: `python3 -m pip wheel . -w /private/tmp/niimath-wheel-audit --no-deps`, install into `/private/tmp/niimath-wheel-venv`, and release smoke with `--expect-bsd --expect-zstd`.
- Release-config macOS wheel path passed: `python3 -m pip wheel . -w /private/tmp/niimath-wheel-audit-staticzstd --no-deps -Ccmake.define.OPENMP_XCODE=ON -Ccmake.define.ENABLE_ZSTD=ON`, install into `/private/tmp/niimath-wheel-staticzstd-venv`, release smoke with `--expect-bsd --expect-zstd`, and `otool -L` confirmed no Homebrew zstd runtime dependency.
- Installed release-config wheel against the user's real file passed: `/private/tmp/niimath-wheel-staticzstd-venv/bin/niimath /Users/chris/src/tmp/brainchop/strokeSubject.nii.gz -conform -gz 0 /private/tmp/brainchop-wheel-conform.nii.gz -odt char` emitted no short-read warning; output header is datatype `2`, `nbyper=1`, `nvox=16777216`.
- Cibuildwheel documentation check confirmed platform-specific variables such as `CIBW_CONFIG_SETTINGS_WINDOWS` and `CIBW_BEFORE_ALL_LINUX`, and confirmed Linux builds use manylinux/musllinux containers with different package managers; release CI was adjusted accordingly.
- Release workflow YAML parses with Ruby/Psych after the multiline Linux zstd install command.
- Direct `src` CMake reconfigured after the static-name expansion: `cmake -S src -B /private/tmp/niimath-cmake-src-zstd-names-audit -DENABLE_ZSTD=ON` selected `/opt/homebrew/opt/zstd/lib/libzstd.a`, built successfully, and passed the release smoke.
- JS package tests still pass after the release review: `bun test ./tests` reports 8 pass, 0 fail.
- `git diff --check` passed.
