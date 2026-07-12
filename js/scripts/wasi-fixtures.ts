/* wasi-fixtures.ts — portable clean-checkout gate for the WASI test suite.
 *
 * Replaces the earlier POSIX-shell one-liner (which failed under the project's Windows shell and
 * only checked the large fixture). Runs via `bun run` on every platform. It: (1) verifies EVERY
 * expected fixture, regenerating the full set through the deterministic Python generators if any is
 * missing (trying `python3` then `python`); (2) records per-fixture SHA-256 and warns on drift from
 * a prior run (informational — regeneration is deterministic per numpy build, so hashes are not a
 * hard gate across numpy versions); (3) reports the native-binary and WASI-reactor prerequisites the
 * suite needs. Kept out of runtime code — this is test tooling only.
 */
import { existsSync, readFileSync, writeFileSync } from "node:fs";
import { spawnSync } from "node:child_process";
import { createHash } from "node:crypto";
import { fileURLToPath } from "node:url";
import { resolve } from "node:path";

const here = fileURLToPath(new URL(".", import.meta.url)); // js/scripts/
const WASI = resolve(here, "../src/wasi");
const FX = resolve(WASI, "fixtures");
const NII = resolve(here, `../../src/niimath${process.platform === "win32" ? ".exe" : ""}`);
const REACTOR = resolve(here, "../src/niimath-wasi.wasm");
const MANIFEST = resolve(FX, ".sha256.json");
const requireBuilds = process.argv.includes("--require-builds");

// gen_fixture.py writes the first; gen_test_fixtures.py writes the rest.
const EXPECTED = [
  "t2w_256x256x192_int16.nii",
  "bold_4x4x4x20.nii", "dwi_6x6x6x7.nii", "dwi.bvec", "dwi.bval", "dwi_mask.nii",
  "t1_small.nii", "seg_small.nii", "mov_small.nii", "tmpl_small.nii", "mask_small.nii",
  "t2w_small_nifti2.nii",
];

function runGen(gen: string): boolean {
  for (const py of ["python3", "python"]) {
    const r = spawnSync(py, [resolve(WASI, gen)], { stdio: "inherit" });
    if (r.status === 0) return true;
    // Try the other interpreter both when this one is absent and when it lacks numpy. A machine can
    // legitimately have `python3` and `python` backed by different environments.
    continue;
  }
  return false;
}

const missing = EXPECTED.filter((f) => !existsSync(resolve(FX, f)));
if (missing.length) {
  console.log(`[wasi-fixtures] regenerating (missing ${missing.length}: ${missing.slice(0, 3).join(", ")}${missing.length > 3 ? ", …" : ""})`);
  for (const gen of ["gen_fixture.py", "gen_test_fixtures.py"]) {
    if (!runGen(gen)) {
      console.error(`[wasi-fixtures] '${gen}' failed. This gate needs python3 (or python) with numpy — 'pip install numpy'.`);
      process.exit(1);
    }
  }
}
const still = EXPECTED.filter((f) => !existsSync(resolve(FX, f)));
if (still.length) {
  console.error(`[wasi-fixtures] still missing after generation: ${still.join(", ")}`);
  process.exit(1);
}

// Record hashes; warn (do not fail) on drift from a prior manifest — numpy build differences can
// move the last ULP of the phantom data, which is not a suite failure.
const hashes: Record<string, string> = {};
for (const f of EXPECTED) hashes[f] = createHash("sha256").update(readFileSync(resolve(FX, f))).digest("hex");
if (existsSync(MANIFEST)) {
  const prev = JSON.parse(readFileSync(MANIFEST, "utf8")) as Record<string, string>;
  const drift = EXPECTED.filter((f) => prev[f] && prev[f] !== hashes[f]);
  if (drift.length) console.warn(`[wasi-fixtures] hash drift (numpy build change?): ${drift.join(", ")}`);
}
writeFileSync(MANIFEST, JSON.stringify(hashes, null, 2) + "\n");

const missingBuilds: string[] = [];
if (!existsSync(NII)) missingBuilds.push("native binary (build with `make -C ../src`)");
if (!existsSync(REACTOR)) missingBuilds.push("WASI reactor (build with `make -C ../src wasm-wasi`)");
if (missingBuilds.length) {
  const message = `[wasi-fixtures] missing test prerequisite(s): ${missingBuilds.join(", ")}`;
  if (requireBuilds) {
    console.error(message);
    process.exit(1);
  }
  console.warn(message);
}
console.log(`[wasi-fixtures] ${EXPECTED.length} fixtures present; hashes recorded.`);
