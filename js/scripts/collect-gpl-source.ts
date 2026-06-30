// Assemble the COMPLETE CORRESPONDING SOURCE for the GPL build of niimath
// (dist/niimath-gpl.wasm) into js/corresponding-source/, so it can be shipped
// inside the npm tarball. This discharges GPL-2 §3(a) — "accompany the binary
// with the complete corresponding source" — directly, instead of a §3(b) written
// offer.
//
// "Complete corresponding source" = every source file compiled/linked into
// niimath-gpl.wasm (the BSD niimath C sources AND the GPL spm_coreg sources) plus
// the scripts used to control its compilation (the src Makefile/CMake and the js
// build wiring). System libraries provided by the Emscripten toolchain (zlib) are
// excluded under GPL-2's System Library exception.
//
// Only generated build artifacts and VCS/dependency dirs are excluded; everything
// needed to rebuild the binary from scratch is included. Run after the GPL wasm
// has been built (see makeWasmGpl / esbuild.config.ts).
import { cpSync, rmSync, mkdirSync, writeFileSync, existsSync } from 'node:fs';
import path from 'node:path';

const repoRoot = path.resolve('..');
const dest = path.resolve('corresponding-source');

// Generated binaries / objects / VCS that are NOT corresponding source.
const EXCLUDE_NAMES = new Set([
  'niimath', 'niimath.js', 'niimath.wasm', 'niimath-gpl.js', 'niimath-gpl.wasm',
  'al_wasm.o', 'pn_wasm.o', 'allineate.o', 'powell_newuoa.o',
  'niimathOperators.json', // generated from CLI help at build time
  '.DS_Store', '.git', 'node_modules', 'dist', 'corresponding-source',
]);

function keep(src: string): boolean {
  const base = path.basename(src);
  if (EXCLUDE_NAMES.has(base)) return false;
  if (base.endsWith('.o')) return false;
  // generated Emscripten output and copied wasm in js/src
  if (/^niimath(-gpl)?\.(js|wasm)$/.test(base)) return false;
  return true;
}

console.log('Collecting GPL corresponding source into', dest);
rmSync(dest, { recursive: true, force: true });
mkdirSync(dest, { recursive: true });

// 1. The C sources + build scripts (src/, includes the GPL submodule at src/GPL).
cpSync(path.join(repoRoot, 'src'), path.join(dest, 'src'), {
  recursive: true,
  filter: (src) => keep(src),
});

// Sanity-check the GPL sources actually came along (submodule must be initialized).
if (!existsSync(path.join(dest, 'src', 'GPL', 'spm_coreg.c'))) {
  throw new Error('GPL sources missing from corresponding-source — is the src/GPL submodule initialized?');
}

// 2. Top-level build wiring / licenses.
for (const f of ['CMakeLists.txt', 'LICENSE', 'README.md']) {
  const p = path.join(repoRoot, f);
  if (existsSync(p)) cpSync(p, path.join(dest, f));
}

// 3. The js wrapper sources + build config (the scripts that produce the published
//    JS and that invoke the wasm build). Copy specific entries (not the whole js/
//    dir, which would contain this very output directory). Excludes generated
//    wasm/js and deps via the filter.
const jsRoot = path.resolve('.');
const JS_ENTRIES = [
  'package.json', 'tsconfig.json', 'esbuild.config.ts', 'server.ts', 'index.html',
  'README.md', 'LICENSE', 'LICENSE.GPL-2.0.txt', 'GPL-NOTICE.md',
  'src', 'scripts', 'tests',
];
mkdirSync(path.join(dest, 'js'), { recursive: true });
for (const entry of JS_ENTRIES) {
  const from = path.join(jsRoot, entry);
  if (!existsSync(from)) continue;
  cpSync(from, path.join(dest, 'js', entry), { recursive: true, filter: (src) => keep(src) });
}

writeFileSync(
  path.join(dest, 'README.SOURCE.md'),
  `# Complete corresponding source — @niivue/niimath GPL build

This directory is the complete corresponding source (GPL-2 §3(a)) for the GPL
build shipped in this package: \`dist/niimath-gpl.js\` and \`dist/niimath-gpl.wasm\`,
reached via \`import { Niimath } from '@niivue/niimath/gpl'\`.

It contains every source file linked into \`niimath-gpl.wasm\` — the BSD-2-Clause
niimath C sources (\`src/*.c\`, \`src/*.h\`) and the GPL-2 spm_coreg sources
(\`src/GPL/\`) — plus the scripts used to control compilation (\`src/Makefile\`,
\`src/CMakeLists.txt\`, and the \`js/\` build wiring). The zlib library used at link
time is supplied by the Emscripten toolchain and is excluded under GPL-2's System
Library exception.

## Rebuild

\`\`\`bash
cd js
bun install
bun run makeWasmGpl   # emits src/niimath-gpl.{js,wasm}; needs emscripten (emcc)
bun run build         # bundles dist/index-gpl.js + dist/niimath-gpl.*
\`\`\`

See \`js/LICENSE\`, \`js/LICENSE.GPL-2.0.txt\`, and \`js/GPL-NOTICE.md\` for licensing.
`,
);

console.log('GPL corresponding source collected.');
