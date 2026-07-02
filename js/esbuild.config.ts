import * as esbuild from "esbuild";
import type { BuildOptions } from "esbuild";
import { copyFileSync, existsSync, mkdirSync, rmSync, writeFileSync } from "fs";
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

// The GPL build is optional: it requires the GPL submodule and a `GPL=1 make wasm`
// (see the makeWasmGpl npm script) to have produced src/niimath-gpl.js. Only build
// and ship the GPL entry points when that artifact exists, so a plain BSD build
// (without the GPL submodule) still succeeds.
const hasGpl =
  existsSync("./src/niimath-gpl.js") && existsSync("./src/niimath-gpl.wasm");

if (hasGpl) {
  // Build GPL index-gpl.ts
  await esbuild.build({
    ...commonOptions,
    entryPoints: ["./src/index-gpl.ts"],
    outfile: "./dist/index-gpl.js",
    loader: {
      ".json": "json",
    },
  });

  // Build GPL worker-gpl.ts
  await esbuild.build({
    ...commonOptions,
    entryPoints: ["./src/worker-gpl.ts"],
    outfile: "./dist/worker-gpl.js",
    external: ["./niimath-gpl.js"], // Keep niimath-gpl.js as external import
  });
} else {
  console.warn(
    "Warning: src/niimath-gpl.js not found — skipping GPL build. Run `bun run makeWasmGpl` (needs the GPL submodule) to produce it.",
  );
}

// Generate TypeScript declarations using tsc with tsconfig. Declaration-only is
// important: plain `tsc --project` emits JS into dist and overwrites esbuild's
// browser-ready bundles with extensionless imports (`./core`), which fail under
// native ESM in browsers.
try {
  execSync("bun x tsc --project tsconfig.json --emitDeclarationOnly", { stdio: "inherit" });
} catch (error) {
  const err = error as Error;
  console.warn("Warning: Could not generate .d.ts files:", err.message);
}

// copy BSD niimath.wasm, niimath.js, niimathOperators.json to dist folder
copyFileSync("./src/niimath.wasm", "./dist/niimath.wasm");
copyFileSync("./src/niimath.js", "./dist/niimath.js");
copyFileSync("./src/niimathOperators.json", "./dist/niimathOperators.json");

// copy GPL niimath-gpl.wasm / niimath-gpl.js when present
if (hasGpl) {
  copyFileSync("./src/niimath-gpl.wasm", "./dist/niimath-gpl.wasm");
  copyFileSync("./src/niimath-gpl.js", "./dist/niimath-gpl.js");
  // Assemble the complete corresponding source (GPL-2 requirement into the package.
  execSync("bun run scripts/collect-gpl-source.ts", { stdio: "inherit" });
} else {
  const unavailableMessage =
    "The @niivue/niimath GPL build was not generated. Run `bun run makeWasmGpl` before importing @niivue/niimath/gpl.";
  rmSync("./dist/niimath-gpl.wasm", { force: true });
  rmSync("./corresponding-source", { recursive: true, force: true });
  writeFileSync(
    "./dist/index-gpl.js",
    `export { dataTypes } from "./core.js";\nexport class Niimath {\n  constructor() {\n    throw new Error(${JSON.stringify(unavailableMessage)});\n  }\n}\n`,
  );
  writeFileSync(
    "./dist/niimath-gpl.js",
    `export default function niimathGplUnavailable() {\n  return Promise.reject(new Error(${JSON.stringify(unavailableMessage)}));\n}\n`,
  );
}

console.log("Build completed!");
