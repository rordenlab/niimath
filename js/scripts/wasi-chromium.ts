// wasi-chromium.ts — M6 Chromium decision run + Bun⇄Chromium agreement.
// Bundles the WASI-C browser entry, serves wasm + fixture + page over localhost, drives real
// Chromium via Playwright to run the pipelines, then reproduces the SAME payload hashes in Bun and
// checks that both runtimes reach the same correctness conclusion. Writes a temporary report by
// default; set WASI_CHROMIUM_REPORT to retain it at a chosen path.
// Run from js/: bun run scripts/wasi-chromium.ts
import { chromium } from "playwright";
import { readFileSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { fileURLToPath } from "node:url";
import { WasiRunner } from "../src/wasi/WasiRunner";

const ROOT = fileURLToPath(new URL("../../", import.meta.url));
const WASM = `${ROOT}js/src/niimath-wasi.wasm`;
const FIX = `${ROOT}js/src/wasi/fixtures/t2w_256x256x192_int16.nii`;
const PIPELINES = [
  { name: "add_mul_sqrt", argv: ["-add", "10", "-mul", "2", "-sqrt"] },
  { name: "gauss_s3", argv: ["-s", "3"] },
  { name: "edge", argv: ["-edge"] },
  { name: "fmean_boxv5", argv: ["-kernel", "boxv", "5", "-fmean"] },
];
const REPS = 3;

// --- 1) bundle the browser entry ---
const build = await Bun.build({ entrypoints: [`${ROOT}js/scripts/wasi-browser-entry.ts`], target: "browser", format: "iife", minify: true });
const bundle = await build.outputs[0].text();

// --- 2) serve bundle + wasm + fixture + page ---
const wasmBytes = readFileSync(WASM);
const fixBytes = readFileSync(FIX);
const page_html = `<!doctype html><html><head><meta charset="utf-8"></head><body><script>${bundle}</script></body></html>`;
const server = Bun.serve({
  port: 0,
  fetch(req) {
    const u = new URL(req.url);
    if (u.pathname === "/") return new Response(page_html, { headers: { "content-type": "text/html" } });
    if (u.pathname === "/niimath-wasi.wasm") return new Response(wasmBytes, { headers: { "content-type": "application/wasm" } });
    if (u.pathname === "/fixture.nii") return new Response(fixBytes);
    return new Response("not found", { status: 404 });
  },
});
const base = `http://localhost:${server.port}`;

// --- 3) drive Chromium ---
const browser = await chromium.launch();
const page = await browser.newPage();
const logs: string[] = [];
page.on("console", (m) => logs.push(m.text()));
page.on("pageerror", (e) => logs.push("PAGEERROR " + e.message));
await page.goto(base);
const chromeResults = await page.evaluate(async ({ base, pipelines, reps }) => {
  const wasm = await (await fetch(base + "/niimath-wasi.wasm")).arrayBuffer();
  const fixture = await (await fetch(base + "/fixture.nii")).arrayBuffer();
  return await (window as any).runWasiPipelines(wasm, fixture, pipelines, reps);
}, { base, pipelines: PIPELINES, reps: REPS });
await browser.close();
server.stop();

// --- 4) reproduce payload hashes in Bun for cross-runtime equality ---
function hashPayload(out: Uint8Array): number {
  let h = 2166136261 >>> 0;
  const step = Math.max(1, Math.floor(out.byteLength / 65536));
  for (let j = 0; j < out.byteLength; j += step) h = Math.imul(h ^ out[j], 16777619) >>> 0;
  return h;
}
const runner = await WasiRunner.create(new Uint8Array(wasmBytes));
const agreement: Record<string, any> = {};
let allAgree = true;
for (const { name, argv } of PIPELINES) {
  const res = await runner.runFiles({ argv: ["in.nii", ...argv, "out.nii"], inputs: { "in.nii": new Uint8Array(fixBytes) }, outputs: ["out.nii"] });
  const bunHash = res.files["out.nii"] ? hashPayload(res.files["out.nii"]) : -1;
  const c = chromeResults[name];
  const agree = c.exitCode === 0 && res.exitCode === 0 && c.outHash === bunHash && c.outLen === res.files["out.nii"]?.byteLength;
  if (!agree) allAgree = false;
  agreement[name] = { chromium_exit: c.exitCode, bun_exit: res.exitCode, chromium_hash: c.outHash, bun_hash: bunHash, hashes_match: c.outHash === bunHash, chromium_p50_ms: c.p50 };
  console.log(`${agree ? "AGREE" : "DIFFER"} ${name.padEnd(14)} chromium_exit=${c.exitCode} hashMatch=${c.outHash === bunHash} chromium_p50=${c.p50}ms`);
}

const out = { runtime_agreement_all: allAgree, note: "Bun⇄Chromium: same reactor, same host VFS; payload hashes must match and both must succeed.", pipelines: agreement, browser_logs: logs.slice(-10) };
const report = process.env.WASI_CHROMIUM_REPORT ?? join(tmpdir(), "niimath-wasi-chromium.json");
writeFileSync(report, JSON.stringify(out, null, 2) + "\n");
console.log(`\n${report} written; runtime_agreement_all=${allAgree}`);
if (!allAgree) process.exit(1);
