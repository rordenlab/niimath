// Load the Emscripten-generated JavaScript, which will handle the WASM binary loading.
// The worker is of type "module" so that it can use ES6 module syntax.
// Importantly, the emscripten code must be compiled with: -s EXPORT_ES6=1 -s MODULARIZE=1.
// This allows proper module bundlers to import the worker and wasm properly with code splitting. 
import Module from './niimath.js';

// initialise an instance of the Emscripten Module so that
// it is ready when the worker receives a message.
// We keep a reference to the module so that the worker can reuse it 
// for all subsequent calls without having to reinitialise it (which could be slow due to the WASM loading)
let mod = null
Module().then((initializedMod) => {
  mod = initializedMod
  // Send a ready message once initialization is complete
  // so we can signal to the main thread that the worker is ready.
  // The Niimath.init() method will wait for this message before resolving the promise.
  self.postMessage({ type: 'ready' });
})

// error handler in the worker
self.onerror = function (message, error) {
  self.postMessage({ type: 'error', message: message, error: error ? error.stack : null });
};
// unhandled promise rejection handler in the worker
self.onunhandledrejection = function (event) {
  self.postMessage({ type: 'error', message: event.reason ? event.reason.message : 'Unhandled rejection', error: event.reason ? event.reason.stack : null });
};

const handleMessage = (e) => {
  try {
    const file = e.data.blob;
    const args = e.data.cmd;
    const outName = e.data.outName || 'out.nii';

    if (!file || args.length < 1) {
      throw new Error("Expected a file and at least one command");
    }

    const inName = file.name;
    const fr = new FileReader();

    fr.readAsArrayBuffer(file);

    fr.onloadend = function () {
      const data = new Uint8Array(fr.result);
      // if the module is loaded and not null,
      // then we can proceed with the image processing
      if (mod) {
        // Create a virtual file in the Emscripten filesystem
        mod.FS_createDataFile(".", inName, data, true, true);
        // Ensure args is an array
        if (!Array.isArray(args)) {
          throw new Error("Expected args to be an array");
        }
        // call the main niimath entry point with the arguments.
        // equivalent to calling `niimath <in_file.nii> <args> <out_file.nii>` from the command line
        const exitCode = mod.callMain(args);

        // read the output file from the Emscripten filesystem
        const out_bin = mod.FS_readFile(outName);

        // binary output file from niimath wasm: nii or mz3
        const outputFile = new Blob([out_bin], { type: 'application/sla' });

        // send a message back to the main thread with the output file, exit code and output file name
        self.postMessage({ blob: outputFile, outName: outName, exitCode: exitCode });

        // Free virtual files
        mod.FS_unlink(inName);

        // if the input and output file names are different, then delete the output file
        if (inName !== outName) {
          mod.FS_unlink(outName);
        }

      } else {
        throw new Error("WASM module not loaded yet!");
      }

    };
  } catch (err) {
    // Send error details back to the main thread
    self.postMessage({ type: 'error', message: err.message, error: err.stack });
  }
}

// Handle messages from the main thread
self.addEventListener('message', handleMessage, false);