// Load the Emscripten-generated JavaScript, which will handle the WASM binary loading.
// The worker is of type "module" so that it can use ES6 module syntax.
// Importantly, the emscripten code must be compiled with: -s EXPORT_ES6=1 -s MODULARIZE=1.
// This allows proper module bundlers to import the worker and wasm properly with code splitting.
import Module from './niimath.js';

// Emscripten module interface
interface EmscriptenModule {
  FS_createDataFile(
    parent: string,
    name: string,
    data: Uint8Array,
    canRead: boolean,
    canWrite: boolean
  ): void;
  callMain(args: string[]): number;
  FS_readFile(filename: string): Uint8Array;
  FS_unlink(filename: string): void;
}

interface WorkerInputMessage {
  blob: Blob;
  cmd: string[];
  outName?: string;
}

interface WorkerReadyMessage {
  type: 'ready';
}

interface WorkerErrorMessage {
  type: 'error';
  message: string;
  error: string | null;
}

interface WorkerSuccessMessage {
  blob: Blob;
  outName: string;
  exitCode: number;
}

// initialise an instance of the Emscripten Module so that
// it is ready when the worker receives a message.
// We keep a reference to the module so that the worker can reuse it
// for all subsequent calls without having to reinitialise it (which could be slow due to the WASM loading)
let mod: EmscriptenModule | null = null;

(Module() as Promise<EmscriptenModule>).then((initializedMod: EmscriptenModule) => {
  mod = initializedMod;
  // Send a ready message once initialization is complete
  // so we can signal to the main thread that the worker is ready.
  // The Niimath.init() method will wait for this message before resolving the promise.
  const readyMsg: WorkerReadyMessage = { type: 'ready' };
  self.postMessage(readyMsg);
});

// error handler in the worker
self.onerror = function (event: Event | string) {
  const errorMsg: WorkerErrorMessage = {
    type: 'error',
    message: typeof event === 'string' ? event : (event as ErrorEvent).message ?? 'Unknown error',
    error: (event as ErrorEvent).error?.stack ?? null
  };
  self.postMessage(errorMsg);
};

// unhandled promise rejection handler in the worker
self.onunhandledrejection = function (event: PromiseRejectionEvent) {
  const errorMsg: WorkerErrorMessage = {
    type: 'error',
    message: event.reason?.message ?? 'Unhandled rejection',
    error: event.reason?.stack ?? null
  };
  self.postMessage(errorMsg);
};

const handleMessage = (e: MessageEvent<WorkerInputMessage>) => {
  try {
    const file = e.data.blob;
    const args = e.data.cmd;
    const outName = e.data.outName ?? 'out.nii';

    if (!file || args.length < 1) {
      throw new Error("Expected a file and at least one command");
    }

    const inName = (file as File).name ?? 'input.nii';
    const fr = new FileReader();

    fr.readAsArrayBuffer(file);

    fr.onloadend = function () {
      const data = new Uint8Array(fr.result as ArrayBuffer);
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
        const outputFile = new Blob([out_bin.buffer as ArrayBuffer], { type: 'application/sla' });

        // send a message back to the main thread with the output file, exit code and output file name
        const successMsg: WorkerSuccessMessage = {
          blob: outputFile,
          outName: outName,
          exitCode: exitCode
        };
        self.postMessage(successMsg);

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
    const error = err as Error;
    const errorMsg: WorkerErrorMessage = {
      type: 'error',
      message: error.message,
      error: error.stack ?? null
    };
    self.postMessage(errorMsg);
  }
};

// Handle messages from the main thread
self.addEventListener('message', handleMessage, false);
