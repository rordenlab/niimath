# WASI-C niimath backend

A lean, zlib-free, single-threaded WASI reactor of niimath's **computational BSD** feature set: core math, multi-image operations, allineate/deface, dtifit, QC, conform and butterworth. It is the experimental WASI-C backend. The shipped default backend remains the Emscripten build (`../niimath.js` and `../niimath.wasm`).

## What it is and is not

- **Is:** the whole scoped compute engine compiled with `zig cc -target wasm32-wasi -mexec-model=reactor`. Zig bundles wasi-libc, so no wasi-sdk install is needed. File I/O flows through real wasi-libc stdio, serviced by a small host-side in-memory VFS. NIfTI gzip is done host-side (`CompressionStream` or `node:zlib`), so the module sees only raw `.nii` bytes.
- **Is not:** a drop-in replacement yet. **Bitmap (`-bitmap`) and mesh (`-mesh`) are not compiled.** They fail clearly with a nonzero exit and never silently do nothing. There is no GPL, zstd or OpenMP.

## Files

| File | Role |
|---|---|
| `../../../src/wasi_shim.c` | Reactor entry (`nii_run`). The only C added; no algorithm changes. |
| `WasiVfs.ts` | Narrow host-side in-memory VFS: flat relative namespace, quota, stdout/stderr capture. |
| `WasiImports.ts` | The 15 frozen WASI Preview-1 imports over the VFS. `proc_exit` throws `WasiExit`. |
| `WasiRunner.ts` | Reactor lifecycle, argv marshaling, instance recycling, `runFiles` and `runFilesGz`. |
| `gzip.ts` | Host-side magic-byte gzip (CompressionStream, with a fallback to `node:zlib` in Bun). |
| `featureParity.test.ts`, `compression.test.ts` | Bun conformance, lifecycle, quota and gzip tests. |

## Build

```bash
make -C ../../../src wasm-wasi        # -> ../niimath-wasi.wasm (+ import-surface check)
make -C ../../../src wasm-emcc-core   # feature-matched Emscripten comparison artifact

cd ../..
bun run test:wasi                     # fixtures + Bun conformance/lifecycle tests
bun run wasiNodeSmoke                 # plain-Node stream/gzip/quota smoke
bun run wasiChromium                  # optional real-browser correctness/agreement check
```

`wasm-wasi` verifies the toolchain (Zig ≥ 0.16), builds the reactor with `-O3 -flto -ffast-math --strip-all`, and fails if the WASI import surface grows beyond [`import-manifest.json`](./import-manifest.json).

## Use

```ts
import { WasiRunner } from "./WasiRunner";
const runner = await WasiRunner.create(wasmBytes);

// low-level multi-file API
const r = await runner.runFiles({
  argv: ["in.nii", "-add", "1", "out.nii"],
  inputs: { "in.nii": rawNiftiBytes },
  outputs: ["out.nii"],
});
// r.exitCode, r.files["out.nii"], r.stdout, r.stderr

// gzip-aware: gzip inputs decompressed + argv rewritten; .gz outputs recompressed
const g = await runner.runFilesGz({
  argv: ["in.nii.gz", "-s", "3", "out.nii.gz"],
  inputs: { "in.nii.gz": gzippedBytes },
  outputs: ["out.nii.gz"],
});
```

`create()` calls `_initialize` once. A normal run reuses the instance after `reset()`. A `proc_exit` (the error path) recreates the instance automatically. A multi-output command (dtifit writes 11 files) must list every output name in `outputs`. One runner is single-flight: await a run before you call another run, `reset()` or `addFile()`, and use one runner per concurrent job. Returned file buffers are detached copies, and the caller may mutate them. `runFilesGz()` clears its raw staging state before compression and return, so keep its returned buffers instead of expecting `readFile()` to expose them.

## Lifecycle contract

1. `_initialize` runs exactly once per instance.
2. A normal `nii_run` return permits reuse after `reset()`, which clears the VFS.
3. A `proc_exit` throws `WasiExit`. The exit code is captured, and the instance is discarded and recreated, because the C cleanup was skipped and the linear memory and libc state are untrusted.

## Status

Functionally complete and verified: feature parity and payload/header conformance against the native binary pass (`bun run test:wasi`). It is a correct, portable, glue-free backend kept for audit and portability. Performance evidence is deferred until the current source is committed. The benchmark lives in `niimath_tests/wasm_benchmark` and is not a release-quality claim.
