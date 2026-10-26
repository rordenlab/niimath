# Releasing niimath

This is the procedure for cutting a release. It is kept out of `AGENTS.md` because it is needed a few times a year, not every session. Each guard below exists because its failure was paid for once.

## The three release outputs

Each output has its own mechanism, and the three are independent. The npm version is NOT the CLI `kMTHdate`.

| Output | Mechanism | Trigger |
|---|---|---|
| Python wheels to PyPI | `release.yml` (cibuildwheel, `deploy-pypi`) | a `v*.*.*` tag. A push to `master` runs the same jobs but uploads nothing |
| Standalone CLI binaries as GitHub Release assets (`niimath_{win,lnx,macos}.zip`) | `release-binaries.yml`, with the built-in `GITHUB_TOKEN` | a `v*.*.*` tag, or `workflow_dispatch` with the tag as input |
| npm package `@niivue/niimath` (BSD-only) on npmjs | `js.yml`, npm trusted publishing (OIDC) | a push to `master` whose `js/package.json` version is not yet on npm |

AppVeyor only builds and tests. It publishes nothing.

## Procedure

### Before you tag

1. Check that `release-checks.yml` passed on the pull request. It runs `vermin --target=<min> --violations` on `.github/scripts/release_smoke.py`, reading `<min>` from `pyproject.toml`. It runs on every `pull_request`, so it blocks before a merge to `master` triggers a release, and it runs again on the release push or tag.
2. If you changed `release_smoke.py`, keep it stdlib-only and clean for the package-minimum Python (currently `>=3.7`). See the guards below.
3. To reproduce the wheel test locally, run `pipx run cibuildwheel --only <id>` or `python3.7 .github/scripts/release_smoke.py <binary>`.

### Python wheels and CLI binaries

1. Push the `v*.*.*` tag. `release.yml` builds the wheels and uploads them to PyPI. `release-binaries.yml` attaches the CLI binaries to the GitHub release.
2. PyPI rejects a duplicate version. Do NOT re-push an existing tag and expect a new upload. Cut a new version instead.
3. To attach binaries to an existing release again WITHOUT re-tagging, and without re-triggering the one-shot PyPI upload, run `release-binaries.yml` via `workflow_dispatch` with the tag as input.

### npm package

`js.yml` publishes on a push to `master` when `js/package.json` holds a version npm does not have. After the WASM build and every `bun` test pass, a separate publish job receives the tested tarball. Only that job has OIDC permission; publishing runs no package scripts. Registry lookup errors stop the job. Auth is npm trusted publishing: no token, no OTP. Configure the trusted publisher on npmjs.com: owner `rordenlab`, repository `niimath`, workflow `js.yml`, with direct `npm publish` allowed. Provenance requires the `repository` field in `js/package.json`.

1. Bump `version` in `js/package.json` in a pull request.
2. Merge it. The push to `master` builds, tests and publishes.
3. A version already on npm is skipped, so other pushes publish nothing. Do not lower the version: npm would tag the lower version `latest`.

To publish by hand (a fallback), run `cd js && bun run build && npm publish --access public` from a clean `master`. `neurolabusc` and `thanayik` have `read-write` access through the `@niivue` org, and 2FA is required (E403 otherwise).

## Guards and why they exist

- **`release_smoke.py` MUST run on the package-minimum Python.** `release.yml` builds wheels with cibuildwheel across five OS/arch runners and runs `release_smoke.py` as the per-wheel `CIBW_TEST_COMMAND` (`pyproject.toml`), inside EVERY target Python from `requires-python` up. This is a different failure surface from the rest of CI, which runs the smoke test only on the runner's modern `python3`. A stdlib symbol newer than the minimum (for example `math.prod` 3.8+, `str.removeprefix` 3.9+, `math.isqrt` 3.8+) passes ordinary CI but fails only on the OLDEST wheel at release. Use the local `_prod()` helper instead of `math.prod`, and avoid other post-minimum stdlib features. `from __future__ import annotations` keeps type hints deferred, so PEP 585/604 annotations are fine.
- **`release-checks.yml` is the fast pre-deploy guard.** Do NOT pass `--eval-annotations` to vermin. The `__future__` import already makes annotations non-evaluated.
- **`release.yml` also runs on every pull request and every Monday, in a reduced shape: ubuntu-latest only, `cp37-*` and `cp38-*` only.** The wheel build used to run only after a merge to master, so it failed at release time. Twice the cause was build-dependency drift with no PR involved: an upstream release changed what pip could resolve under cibuildwheel's per-Python constraint pins. The oldest Pythons carry the most aged pins, and the pins are per Python, not per OS, so one Linux leg reproduces the failure the whole matrix would hit. Keep the reduced shape in the same workflow as the full one, so the two cannot drift. `fail-fast` is off so a release reports every failing leg. The PyPI job runs on pushes only.
- **`pyproject.toml` build requires stay unconditional and extra-free.** `scikit-build-core[pyproject]` named an extra that scikit-build-core 1.x no longer defines, and pip 25 then failed to backtrack around an unrelated pin. The floor stays at `scikit-build-core>=0.5`: raising it to 0.11 broke the Linux cp37 wheels, because 0.11 requires Python 3.8 while cp37 resolves to 0.10.7.
- **Only `release_smoke.py` is minimum-Python-constrained.** `make_fixtures.py` and `reg_quality.py` run only on the runner's modern `python3`, never inside a wheel.
- **Never gate a release on a rotating PAT.** Prefer the auto-provisioned `GITHUB_TOKEN`, as `release-binaries.yml` does.
- **The macOS binary is a universal build with static zstd and no OpenMP.** The `release-binaries.yml` GitHub Actions job builds separate x86_64 and arm64 slices and joins them with `lipo`. Homebrew ships only the runner's native arch of libzstd, so the job builds a **universal static `libzstd.a` from source** (`-arch x86_64 -arch arm64`) and points both slices at it via `ZSTD_ROOT`/`CMAKE_PREFIX_PATH`. The binary is self-contained, with no runtime `libzstd.dylib`. `src/CMakeLists.txt` resolves the configured library to a full path. The universal build passes `-DUSE_OPENMP=OFF`, because a single-arch libomp cannot link into a universal binary.
