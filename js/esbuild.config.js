import * as esbuild from 'esbuild';
import { copyFileSync } from 'fs';

await esbuild.build({
  entryPoints: ['./src/index.js'],
  outfile: './dist/index.js',
  bundle: true,
  format: 'esm',
  target: ['es2020'],
  minify: false,
  define: {
    'process.env.NODE_ENV': '"production"',
  },
});

// copy worker.js, niimath.wasm, niimath.js to dist folder
// (they do not require any processing by esbuild).
// Technically, none of the files in the src folder require processing by esbuild,
// but it does allow minification (optional), and ES version target specification if needed.
// In the future, if we use Typescript, we can use esbuild to transpile the Typescript to JS.
copyFileSync('./src/worker.js', './dist/worker.js');
copyFileSync('./src/niimath.wasm', './dist/niimath.wasm');
copyFileSync('./src/niimath.js', './dist/niimath.js');
console.log('Build completed!');