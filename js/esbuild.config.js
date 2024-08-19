const esbuild = require('esbuild');
const fs = require('fs');

esbuild.build({
  entryPoints: ['./src/index.js'],
  outfile: './dist/index.js',
  bundle: true,
  format: 'esm',
  target: ['es2020'],
  minify: false,
  define: {
    'process.env.NODE_ENV': '"production"',
  },
}).then(() => {
  // copy worker.js, niimath.wasm, niimath.js to dist folder
  // (they do not require any processing by esbuild).
  // Technically, none of the files in the src folder require processing by esbuild, 
  // but it does allow minification (optional), and ES version target specification if needed.
  // In the future, if we use Typescript, we can use esbuild to transpile the Typescript to JS.
  fs.copyFileSync('./src/worker.js', './dist/worker.js');
  fs.copyFileSync('./src/niimath.wasm', './dist/niimath.wasm');
  fs.copyFileSync('./src/niimath.js', './dist/niimath.js');
  console.log('Build completed!');
}).catch(() => process.exit(1));