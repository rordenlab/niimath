// wasi-node-smoke.mjs — M4 Node 20+ functional smoke: real reactor run + gzip, under plain Node.
// Node 22.6+/26 strips TS types natively, so we can import the .ts runner directly.
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { pathToFileURL } from "node:url";
import { WasiRunner } from "../src/wasi/WasiRunner.ts";
import { isGzip, gzip, gunzip } from "../src/wasi/gzip.ts";

// Resolve against the repo root: NIIMATH_ROOT env, else the source-tree location (unbundled run).
const ROOT = process.env.NIIMATH_ROOT
  ? new URL("./", pathToFileURL(resolve(process.env.NIIMATH_ROOT) + "/"))
  : new URL("../../", import.meta.url);
const rd = (p) => new Uint8Array(readFileSync(new URL(p, ROOT)));

const runner = await WasiRunner.create(rd("js/src/niimath-wasi.wasm"));
const fixture = rd("js/src/wasi/fixtures/t2w_small_nifti2.nii");

// 1) raw run
let r = await runner.runFiles({ argv: ["in.nii", "-add", "1", "out.nii"], inputs: { "in.nii": fixture }, outputs: ["out.nii"] });
if (r.exitCode !== 0 || !r.files["out.nii"]) throw new Error("FAIL: raw run");
console.log("node raw run: exit", r.exitCode, "out bytes", r.files["out.nii"].byteLength);

// 2) gzip in + gzip out
const gzIn = await gzip(fixture);
r = await runner.runFilesGz({ argv: ["in.nii.gz", "-add", "1", "out.nii.gz"], inputs: { "in.nii.gz": gzIn }, outputs: ["out.nii.gz"] });
if (r.exitCode !== 0) throw new Error("FAIL: gz run exit " + r.exitCode);
const outGz = r.files["out.nii.gz"];
if (!outGz || !isGzip(outGz)) throw new Error("FAIL: gz output not gzip-framed");
const back = await gunzip(outGz);
console.log("node gzip run: exit", r.exitCode, "decompressed out bytes", back.byteLength);

// 3) failure recovery
r = await runner.runFiles({ argv: ["nope.nii", "-add", "1", "x.nii"], inputs: {}, outputs: ["x.nii"] });
if (r.exitCode === 0) throw new Error("FAIL: malformed returned 0");
r = await runner.runFiles({ argv: ["in.nii", "-add", "1", "ok.nii"], inputs: { "in.nii": fixture }, outputs: ["ok.nii"] });
if (r.exitCode !== 0) throw new Error("FAIL: recovery run");
console.log("node failure-recovery: ok");

// 4) the browser/Node DecompressionStream path enforces one aggregate decompressed-byte cap
const capped = await WasiRunner.create(rd("js/src/niimath-wasi.wasm"),
  { env: { FSLOUTPUTTYPE: "NIFTI" } }, { maxBytes: 12 });
const smallGz = await gzip(new Uint8Array(8));
let quotaRejected = false;
try {
  await capped.runFilesGz({
    argv: ["a.nii.gz", "b.nii.gz"],
    inputs: { "a.nii.gz": smallGz, "b.nii.gz": smallGz },
  });
} catch (error) {
  quotaRejected = error instanceof Error && error.message.includes("exceeds cap");
}
if (!quotaRejected) throw new Error("FAIL: decompression quota was not enforced");
console.log("node decompression quota: ok");

console.log("\nNODE SMOKE: ALL PASS");
