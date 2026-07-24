import * as esbuild from "esbuild";
import type { BuildOptions } from "esbuild";
import { copyFileSync, mkdirSync, rmSync } from "fs";
import { execSync } from "child_process";

rmSync("./dist", { recursive: true, force: true });
mkdirSync("./dist", { recursive: true });

// Common build options
const commonOptions: Partial<BuildOptions> = {
  bundle: true,
  format: "esm",
  target: ["es2020"],
  minify: false,
  define: {
    "process.env.NODE_ENV": '"production"',
  },
};

// Build main (BSD) index.ts
await esbuild.build({
  ...commonOptions,
  entryPoints: ["./src/index.ts"],
  outfile: "./dist/index.js",
  loader: {
    ".json": "json",
  },
});

// Build BSD worker.ts
await esbuild.build({
  ...commonOptions,
  entryPoints: ["./src/worker.ts"],
  outfile: "./dist/worker.js",
  external: ["./niimath.js"], // Keep niimath.js as external import
});

// The published npm package is BSD-only and does NOT ship the optional GPL
// spm_coreg build. The C sources (src/GPL submodule) and the `GPL=1 make` /
// `makeWasmGpl` opt-in remain for local/historical use, but esbuild never emits
// index-gpl/worker-gpl/niimath-gpl here, and package.json exposes no ./gpl export.

// Generate TypeScript declarations using tsc with tsconfig. Declaration-only is
// important: plain `tsc --project` emits JS into dist and overwrites esbuild's
// browser-ready bundles with extensionless imports (`./core`), which fail under
// native ESM in browsers.
execSync("bun x tsc --project tsconfig.json --emitDeclarationOnly", { stdio: "inherit" });

// copy BSD niimath.wasm, niimath.js, niimathOperators.json to dist folder
copyFileSync("./src/niimath.wasm", "./dist/niimath.wasm");
copyFileSync("./src/niimath.js", "./dist/niimath.js");
copyFileSync("./src/niimathOperators.json", "./dist/niimathOperators.json");

console.log("Build completed!");
