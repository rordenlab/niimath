// wasi-browser-entry.ts — bundled to a browser IIFE; exposes the WASI-C backend on window for the
// Playwright/Chromium decision run. Verifies the reactor + host VFS work in a real browser and
// records per-pipeline timings so Bun and Chromium can be checked for the same conclusion.
import { WasiRunner } from "../src/wasi/WasiRunner";

declare global {
  interface Window {
    runWasiPipelines: (
      wasm: ArrayBuffer,
      fixture: ArrayBuffer,
      pipelines: { name: string; argv: string[] }[],
      reps: number,
    ) => Promise<any>;
  }
}

window.runWasiPipelines = async (wasm, fixture, pipelines, reps) => {
  const fx = new Uint8Array(fixture);
  const results: Record<string, any> = {};
  for (const { name, argv } of pipelines) {
    const times: number[] = [];
    let exitCode = -1;
    let outLen = 0;
    let outHash = 0;
    for (let i = 0; i < reps + 1; i++) {
      const t0 = performance.now();
      const runner = await WasiRunner.create(wasm); // cold instance per rep
      const res = await runner.runFiles({
        argv: ["in.nii", ...argv, "out.nii"],
        inputs: { "in.nii": fx },
        outputs: ["out.nii"],
      });
      const t1 = performance.now();
      exitCode = res.exitCode;
      const out = res.files["out.nii"];
      if (out) {
        outLen = out.byteLength;
        // cheap FNV-ish hash of the payload for cross-runtime equality check
        let h = 2166136261 >>> 0;
        const step = Math.max(1, Math.floor(out.byteLength / 65536));
        for (let j = 0; j < out.byteLength; j += step) h = (Math.imul(h ^ out[j], 16777619)) >>> 0;
        outHash = h;
      }
      if (i > 0) times.push(t1 - t0); // drop warmup
    }
    times.sort((a, b) => a - b);
    results[name] = { exitCode, outLen, outHash, p50: +times[Math.floor(times.length / 2)].toFixed(2) };
  }
  return results;
};
