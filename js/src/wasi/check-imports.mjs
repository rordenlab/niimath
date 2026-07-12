// check_imports.mjs <wasm> <import-manifest.json>
// Fails (exit 1) if the reactor's WASI import surface GREW beyond the frozen manifest.
// A shrink is fine (fewer host syscalls to service); a new import must be reviewed first.
import { readFileSync } from "node:fs";

const [wasmPath, manifestPath] = process.argv.slice(2);
if (!wasmPath || !manifestPath) {
  console.error("usage: check_imports.mjs <wasm> <import-manifest.json>");
  process.exit(2);
}
const mod = new WebAssembly.Module(readFileSync(wasmPath));
const actual = WebAssembly.Module.imports(mod)
  .filter((i) => i.module === "wasi_snapshot_preview1")
  .map((i) => i.name)
  .sort();
const frozen = JSON.parse(readFileSync(manifestPath, "utf8")).imports.slice().sort();

const added = actual.filter((n) => !frozen.includes(n));
const removed = frozen.filter((n) => !actual.includes(n));

console.log(`imports: ${actual.length} actual, ${frozen.length} frozen`);
if (removed.length) console.log(`  (removed, ok): ${removed.join(", ")}`);
if (added.length) {
  console.error(`ERROR: new WASI imports not in the frozen manifest: ${added.join(", ")}`);
  console.error("Review + implement the host-side syscall, then update import-manifest.json.");
  process.exit(1);
}
console.log("import surface OK (no un-reviewed growth)");
