import * as esbuild from 'esbuild';
import type { BuildOptions } from 'esbuild';
import { copyFileSync } from 'fs';
import { execSync } from 'child_process';

// Common build options
const commonOptions: Partial<BuildOptions> = {
  bundle: true,
  format: 'esm',
  target: ['es2020'],
  minify: false,
  define: {
    'process.env.NODE_ENV': '"production"',
  },
};

// Build main index.ts
await esbuild.build({
  ...commonOptions,
  entryPoints: ['./src/index.ts'],
  outfile: './dist/index.js',
  loader: {
    '.json': 'json',
  },
});

// Build worker.ts
await esbuild.build({
  ...commonOptions,
  entryPoints: ['./src/worker.ts'],
  outfile: './dist/worker.js',
  external: ['./niimath.js'], // Keep niimath.js as external import
});

// Generate TypeScript declarations using tsc with tsconfig
try {
  execSync('bun x tsc --project tsconfig.json', { stdio: 'inherit' });
} catch (error) {
  const err = error as Error;
  console.warn('Warning: Could not generate .d.ts files:', err.message);
}

// copy niimath.wasm, niimath.js, niimathOperators.json to dist folder
copyFileSync('./src/niimath.wasm', './dist/niimath.wasm');
copyFileSync('./src/niimath.js', './dist/niimath.js');
copyFileSync('./src/niimathOperators.json', './dist/niimathOperators.json');

console.log('Build completed!');