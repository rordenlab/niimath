# @niivue/niimath

`@niivue/niimath` is a JavaScript + WASM library for performing mathemetical operations on NIFTI files. This library is intended to be **used in the browser**, not in a Node.js environment.

> All image processing operations are performed using the WASM build of [niimath](https://github.com/rordenlab/niimath), making it much faster than a pure JavaScript implementation. The image processing takes place in a separate worker thread, so it won't block the main thread in your application.

## Usage

The `@niivue/niimath` JavaScript library offers an object oriented API for working with the `niimath` CLI. Since `niimath` is a CLI tool, the API implemented in `@niivue/niimath` is just a wrapper around the CLI options and arguments. 

### Example: volumes

For example, the [difference of gaussian](https://www.biorxiv.org/content/biorxiv/early/2022/09/17/2022.09.14.507937.full.pdf) command `niimath input.nii -dog 2 3.2 output.nii` can be executed using the following `@niivue/niimath` JavaScript code:

```javascript
import { Niimath } from '@niivue/niimath';

const niimath = new Niimath();
// call the init() method to load the wasm before processing images
await niimath.init();

// 1. selectedFile is a browser File object
// 2. note the use of the final run() method to execute the command. 
// 3. note the use of await. The run method returns a promise that resolves to the output file if the command is successful.
const outFile = await niimath.image(selectedFile).dog(2, 3.2).run();
```

### Registration & defacing

The default (BSD-2-Clause) build includes the affine registration and defacing operations `-allineate` and `-deface` (adapted from AFNI 3dAllineate, public domain). These take other browser `File` objects as arguments:

```javascript
import { Niimath } from '@niivue/niimath';
const niimath = new Niimath();
await niimath.init();

// affine-register `selectedFile` onto a base volume
const registered = await niimath.image(selectedFile).allineate(baseFile).run();

// deface using a template + mask pair
const defaced = await niimath.image(selectedFile).deface(templateFile, maskFile).run();
```

### GPL build (`-spm_coreg`, `-spm_deface`)

A second, larger WASM module adds the optional **GPL-2** SPM coregistration operations (`-spm_coreg`, `-spm_deface`) on top of everything in the BSD build. It is exposed under a separate subpath export so you explicitly opt into the GPL licensing:

```javascript
// GPL-2 build — same API as the default import, plus SPM coregistration
import { Niimath } from '@niivue/niimath/gpl';

const niimath = new Niimath();
await niimath.init();

// rigid-body coregister `selectedFile` onto a reference volume
const coregistered = await niimath.image(selectedFile).spmcoreg(referenceFile).run();

// SPM rigid-body defacing with a template + mask pair
const defaced = await niimath.image(selectedFile).spmDeface(templateFile, maskFile).run();
```

> **Licensing:** importing from `@niivue/niimath/gpl` pulls in GPL-2 code, so a bundle that includes it becomes a GPL-2 combined work. Use the default `@niivue/niimath` import if your project must remain BSD-2-Clause — it still provides `-allineate`/`-deface`, just not the SPM operations. The GPL WASM binary is built from the [`niimath_gpl`](https://github.com/rordenlab/niimath_gpl) submodule; a plain clone without that submodule still builds the BSD package (the GPL entry point is simply omitted).

### Example: meshes

The `@niivue/niimath` library also supports the `-mesh` options available in the `niimath` CLI. However, the JavaScript API is slightly different from the volume processing due to the use of the `-mesh` suboptions. 

```javascript
import { Niimath } from '@niivue/niimath';
const niimath = new Niimath();
await niimath.init();
const outName = 'out.mz3'; // outname must be a mesh format!
const outMesh = await niimath.image(selectedFile)
  .mesh({
    i: 'm', // 'd'ark, 'm'edium, 'b'right or numeric (e.g. 128) isosurface
    b: 1, // fill bubbles
  })
  .run(outName);
/*
Here's the help from the niimath CLI program
The mesh option has multiple sub-options:
 -mesh                    : meshify requires 'd'ark, 'm'edium, 'b'right or numeric isosurface ('niimath bet -mesh -i d mesh.gii')
        -i <isovalue>            : 'd'ark, 'm'edium, 'b'right or numeric isosurface
        -a <atlasFile>           : roi based atlas to mesh
        -b <fillBubbles>         : fill bubbles
        -l <onlyLargest>         : only largest
        -o <originalMC>          : original marching cubes
        -q <quality>             : quality
        -s <postSmooth>          : post smooth
        -r <reduceFraction>      : reduce fraction
        -v <verbose>             : verbose
*/
```

## Installation

To install `@niivue/niimath` in your project, run the following command:

```bash
npm install @niivue/niimath # or bun install @niivue/niimath
```

### To install a local build of the library

Fist, `cd` into the `js` directory of the `niimath` repository.

```bash
# from niimath root directory
cd js
```

To install a local build of the library, run the following command:

```bash
bun run build
```

Then, install the library using the following command:

```bash
npm pack # will create a .tgz file in the root directory
```

Then, install the `@niivue/niimath` library in your application locally using the following command:

```bash
npm install /path/to/niivue-niimath.tgz
```

## Development

Install [Bun](https://bun.com/docs/installation)

First `cd` into the `js` directory of the `niimath` repository.

```bash
# from niimath root directory
cd js
```

To install the dependencies, run the following command:

```bash
bun install
```

To build the library, run the following command

```bash
bun run build
```

> **Note:** `src/niimathOperators.json` and `src/types.ts` are **generated** from the
> niimath CLI help text and are not checked into git. `bun run build` regenerates them
> via its `prebuild` step (`parseHelpText` + `generateTypes`). On a fresh clone, run
> `bun run prebuild` (or `bun run parseHelpText && bun run generateTypes`) once before
> using your editor / `tsc`, otherwise the imports in `src/index.ts` will appear missing.

To run the tests, run the following command:

```bash
bun run test
```

The tests in `tests/` load the built WASM modules from `dist/` directly (via the
in-memory filesystem, no browser Worker), so run `bun run build` first. The GPL
tests (`tests/gpl.test.ts`) automatically **skip** when `dist/niimath-gpl.js` was
not produced (e.g. a clone without the `niimath_gpl` submodule), so the suite still
passes on a BSD-only build.

### Development server with Hot Module Reloading

To start the development server with hot module reloading:

```bash
bun run dev
```

This will start a development server at `http://localhost:3000` with automatic page reloading when source files change.



