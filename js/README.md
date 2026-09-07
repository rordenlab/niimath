# @niivue/niimath

`@niivue/niimath` is a JavaScript and WASM library for mathematical operations on NIfTI files. It is intended for use **in the browser**, not in Node.js.

All image processing runs in the WASM build of [niimath](https://github.com/rordenlab/niimath), which is much faster than a pure JavaScript implementation. The processing runs in a separate worker thread, so it does not block the main thread of your application.

## Installation

```bash
npm install @niivue/niimath # or bun install @niivue/niimath
```

### Install a local build

Build the library, pack it, and install the package in your application:

```bash
# from the niimath root directory
cd js
bun run build
npm pack   # creates a .tgz file in the current directory
npm install /path/to/niivue-niimath.tgz
```

## Usage

The library offers an object-oriented API over the `niimath` CLI. Because `niimath` is a command-line tool, the API is a wrapper around the CLI options and arguments.

### Process a volume

The [difference of Gaussians](https://www.biorxiv.org/content/biorxiv/early/2022/09/17/2022.09.14.507937.full.pdf) command `niimath input.nii -dog 2 3.2 output.nii` becomes:

```javascript
import { Niimath } from '@niivue/niimath';

const niimath = new Niimath();
// call init() to load the WASM before you process images
await niimath.init();

// selectedFile is a browser File object.
// run() executes the command. It returns a promise that resolves to the output file when the command succeeds.
const outFile = await niimath.image(selectedFile).dog(2, 3.2).run();
```

### Register and deface

The default (BSD-2-Clause) build includes the affine registration and defacing operations `-allineate` and `-deface`, adapted from AFNI 3dAllineate (public domain). They take other browser `File` objects as arguments:

```javascript
import { Niimath } from '@niivue/niimath';
const niimath = new Niimath();
await niimath.init();

// affine-register selectedFile onto a base volume
const registered = await niimath.image(selectedFile).allineate(baseFile).run();

// deface with a template and mask pair
const defaced = await niimath.image(selectedFile).deface(templateFile, maskFile).run();
```

`allineate(base, opts?, weight?)` accepts an optional third `File`: a graded weight in the base's space. Its dimensions and world frame must match `base`. Values are normalized to `[0, 1]`. Zero excludes a voxel, and larger values contribute more. Both registration engines honor the weight (the fast engine at its finest stage). Keep the outer head attenuated but nonzero, so the whole-head boundary still anchors the scale.

```javascript
// register selectedFile onto baseFile, with a graded whole-head weight
const registered = await niimath.image(selectedFile).allineate(baseFile, [], weightFile).run();
```

### Use an image as an operand

Two more operations take a `File` operand. `resliceNN(refFile)` reslices the current image onto another image's grid with nearest-neighbor interpolation. `mulImage(imgFile)` multiplies the current image by another image (the generated `mul` takes only a scalar). For example, to reslice a brain mask onto a native grid and apply it:

```javascript
// reslice maskFile (conformed space) onto nativeFile's grid, binarize, save
const maskBlob = await niimath.image(maskFile).resliceNN(nativeFile).bin().run();
// run() returns a Blob; wrap it in a File so it can be used as an operand
const nativeMask = new File([maskBlob], 'nativeMask.nii.gz');
// keep only the masked region of the native image (now on the same grid)
const brain = await niimath.image(nativeFile).mulImage(nativeMask).run();
```

### Create a mesh

The library supports the `-mesh` options of the `niimath` CLI. The JavaScript API differs slightly from volume processing, because `-mesh` has sub-options. Pass them as an object whose keys are the CLI sub-option letters: `i` (isosurface), `a` (atlas file), `b` (fill bubbles), `l` (only largest), `o` (original marching cubes), `q` (quality), `s` (post smooth), `r` (reduce fraction) and `v` (verbose). See the `-mesh` reference in the [niimath README](https://github.com/rordenlab/niimath#-mesh-opts-output).

```javascript
import { Niimath } from '@niivue/niimath';
const niimath = new Niimath();
await niimath.init();
const outName = 'out.mz3'; // the output name must use a mesh format
const outMesh = await niimath.image(selectedFile)
  .mesh({
    i: 'm', // 'd'ark, 'm'edium, 'b'right or numeric (e.g. 128) isosurface
    b: 1, // fill bubbles
  })
  .run(outName);
```

## Worker lifecycle

- Call `await niimath.init()` once before you process images. It spawns one persistent Web Worker. The promise rejects if the worker fails to load or instantiate, for example when the WASM cannot be fetched.
- One `.run()` is in flight per instance at a time. A second overlapping `run()`, or a `run()` before `init()` has resolved, rejects immediately instead of interleaving. Serialize calls (await the previous `run()`), or use a separate instance for each concurrent stream.
- `niimath.dispose()` terminates the worker and releases its WASM heap. It is idempotent, and it rejects any in-flight `init()` or `run()`.
- A worker crash during a `run()` rejects that run and invalidates the worker. A later `image(...).run()` rejects with "Worker not initialized" until you call `init()` again.

## SPM coregistration

The published `@niivue/niimath` package is **BSD-2-Clause only**. It no longer ships the optional GPL-2 SPM coregistration WASM module, because the permissively licensed `-allineate` and `-deface` engine supersedes it. The `-spm_coreg` and `-spm_deface` C sources remain in the [`niimath_gpl`](https://github.com/rordenlab/niimath_gpl) submodule. You can still build them from source for local or historical use (`GPL=1 make`, or `bun run makeWasmGpl` to produce a GPL WASM), but they are not part of the npm package or its exports.

## Development

Install [Bun](https://bun.com/docs/installation). Then, from the `js` directory of the `niimath` repository:

```bash
cd js
bun install      # install the dependencies
bun run build    # build the library
bun run test     # run the tests
bun run dev      # start the development server
```

`src/niimathOperators.json` and `src/types.ts` are **generated** from the niimath CLI help text and are not checked into git. `bun run build` regenerates them in its `prebuild` step (`parseHelpText` + `generateTypes`). On a fresh clone, run `bun run prebuild` (or `bun run parseHelpText && bun run generateTypes`) once before you use your editor or `tsc`. Otherwise the imports in `src/index.ts` appear missing.

The tests in `tests/` load the built WASM module from `dist/` directly, through the in-memory filesystem and without a browser Worker, so run `bun run build` first. The package is BSD-only. The historical GPL binding and test are kept under `gpl-historical/` and are not part of the default suite.

`bun run dev` starts a development server at `http://localhost:3000` with automatic page reloading when source files change.
