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
# TODO: publish to npm
# npm install @niivue/niimath
```

### To install a local build of the library

Fist, `cd` into the `js` directory of the `niimath` repository.

```bash
# from niimath root directory
cd js
```

To install a local build of the library, run the following command:

```bash
npm run build
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

First `cd` into the `js` directory of the `niimath` repository.

```bash
# from niimath root directory
cd js
```

To install the dependencies, run the following command:

```bash
npm install
```

To build the library, run the following command

```bash
npm run build
```

To run the tests, run the following command:

```bash
npm run test
```

### Test using a simple demo

To test that the `@niivue/niimath` library is working correctly, you can run the following command:

```bash
npm run demo
```



