# Releasing niimath

This is the procedure for cutting a release. It is kept out of `AGENTS.md` because it is needed a few times a year, not every session. Each guard below exists because its failure was paid for once.

## The three release outputs

Each output has its own mechanism, and the three are independent. The npm version is NOT the CLI `kMTHdate`.

| Output | Mechanism | Trigger |
|---|---|---|
| Python wheels to PyPI | `release.yml` (cibuildwheel, `deploy-pypi`) | push to `master`, or a `v*.*.*` tag |
| Standalone CLI binaries as GitHub Release assets (`niimath_{win,lnx,macos}.zip`) | `release-binaries.yml`, with the built-in `GITHUB_TOKEN` | a `v*.*.*` tag, or `workflow_dispatch` with the tag as input |
| npm package `@niivue/niimath` (BSD-only) on npmjs | manual `npm publish` from `js/` | by hand; there is no CI publish |

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

Publish by hand from a clean `master`. `neurolabusc` and `thanayik` have `read-write` access through the `@niivue` org. **2FA is required to publish** (E403 otherwise).

1. Bump `version` in `js/package.json`.
2. Run `cd js && bun run build`. `dist/` is gitignored, so a fresh build is REQUIRED before you publish.
3. Run `npm pack --dry-run --json` and check that the tarball is BSD-only, with no GPL or test artifacts. `tests/pack.test.ts` asserts this.
4. Run `npm publish --access public`. With a configured **passkey**, bare `npm publish` succeeds interactively. Otherwise pass `--otp=<code>`, or put an npm Automation or "bypass 2FA" token in `~/.npmrc`.
5. Keep the committed `js/package.json` version equal to the published version.

## Guards and why they exist

- **`release_smoke.py` MUST run on the package-minimum Python.** `release.yml` builds wheels with cibuildwheel across five OS/arch runners and runs `release_smoke.py` as the per-wheel `CIBW_TEST_COMMAND` (`pyproject.toml`), inside EVERY target Python from `requires-python` up. This is a different failure surface from the rest of CI, which runs the smoke test only on the runner's modern `python3`. A stdlib symbol newer than the minimum (for example `math.prod` 3.8+, `str.removeprefix` 3.9+, `math.isqrt` 3.8+) passes ordinary CI but fails only on the OLDEST wheel at release. Use the local `_prod()` helper instead of `math.prod`, and avoid other post-minimum stdlib features. `from __future__ import annotations` keeps type hints deferred, so PEP 585/604 annotations are fine.
- **`release-checks.yml` is the fast pre-deploy guard.** Do NOT pass `--eval-annotations` to vermin. The `__future__` import already makes annotations non-evaluated.
- **Only `release_smoke.py` is minimum-Python-constrained.** `make_fixtures.py` and `reg_quality.py` run only on the runner's modern `python3`, never inside a wheel.
- **Never gate a release on a rotating PAT.** Prefer the auto-provisioned `GITHUB_TOKEN`, as `release-binaries.yml` does.
- **The macOS binary is a universal build with static zstd and no OpenMP.** The AppVeyor macOS job builds separate x86_64 and arm64 slices and joins them with `lipo`. Homebrew ships only the runner's native arch of libzstd, so the job builds a **universal static `libzstd.a` from source** (`-arch x86_64 -arch arm64`) and points both slices at it via `PKG_CONFIG_PATH`. The binary is self-contained, with no runtime `libzstd.dylib`. `src/CMakeLists.txt` resolves the pkg-config result to a full library path. The universal build passes `-DOPENMP_XCODE=OFF`, because a single-arch libomp cannot link into a universal binary, so `OPENMP_XCODE` must remain a consumed, override-able option.
