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
  // Additional files staged into MEMFS by name (e.g. template + mask for
  // -deface/-spm_deface, whose argv tokens are these filenames). Each is read,
  // staged before callMain, and unlinked afterward alongside the main input.
  extraFiles?: { name: string; data: Blob }[];
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

// niimath prints progress chatter to stdout/stderr; Emscripten routes those to
// console.log/console.error (the latter shows red in devtools). Capture both into
// a per-run buffer instead so a successful run is silent, and surface the captured
// text only when a run fails (see the nonzero-exit branch below). The buffer is
// cleared before each callMain.
let runLog: string[] = [];
const captureModule = { print: (s: string) => runLog.push(s), printErr: (s: string) => runLog.push(s) };

(Module(captureModule) as Promise<EmscriptenModule>).then((initializedMod: EmscriptenModule) => {
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
    const extraFiles = e.data.extraFiles ?? [];
    const fr = new FileReader();

    // Surface a read failure as a controlled error. The work runs in fr.onload
    // (success only), so a failed read never reaches the `fr.result as ArrayBuffer`
    // cast with a null, and nothing is staged yet so no cleanup is owed here.
    fr.onerror = function () {
      const errorMsg: WorkerErrorMessage = {
        type: 'error',
        message: `Failed to read input "${inName}": ${fr.error?.message ?? 'unknown read error'}`,
        error: fr.error?.stack ?? null
      };
      self.postMessage(errorMsg);
    };

    fr.readAsArrayBuffer(file);

    fr.onload = async function () {
      // Fires only on a successful read (errors go to fr.onerror), so fr.result is
      // a valid ArrayBuffer here. This callback runs asynchronously, outside the
      // outer try/catch, so it does its own error handling and cleans up staged
      // files in a finally block.
      if (!mod) {
        const errorMsg: WorkerErrorMessage = {
          type: 'error',
          message: "WASM module not loaded yet!",
          error: null
        };
        self.postMessage(errorMsg);
        return;
      }
      const data = new Uint8Array(fr.result as ArrayBuffer);
      let stagedInput = false;
      const stagedExtras: string[] = [];
      try {
        if (!Array.isArray(args)) {
          throw new Error("Expected args to be an array");
        }
        // Create a virtual file in the Emscripten filesystem
        mod.FS_createDataFile(".", inName, data, true, true);
        stagedInput = true;
        // Stage any extra inputs (template/mask) by name so the chain op can open them.
        for (const f of extraFiles) {
          const bytes = new Uint8Array(await f.data.arrayBuffer());
          mod.FS_createDataFile(".", f.name, bytes, true, true);
          stagedExtras.push(f.name);
        }
        // call the main niimath entry point with the arguments.
        // equivalent to calling `niimath <in_file.nii> <args> <out_file.nii>` from the command line
        runLog = []; // discard prior chatter; keep only this run's output for error reporting
        const exitCode = mod.callMain(args);
        // Check the exit code BEFORE reading the output: a failed run may not have
        // written the output file, so reading it would throw a less useful error.
        // On failure, attach the captured niimath stdout/stderr so the error is
        // actionable even though successful runs stay silent.
        if (exitCode !== 0) {
          const detail = runLog.join('\n').trim();
          throw new Error(`niimath exited with code ${exitCode}${detail ? `:\n${detail}` : ''}`);
        }
        // read the output file from the Emscripten filesystem
        const out_bin = mod.FS_readFile(outName);
        // binary output file from niimath wasm: nii or mz3. Copy into an exact-
        // length buffer so the Blob is precisely the file's bytes — the FS_readFile
        // view can in principle sit at an offset/short of its backing ArrayBuffer.
        const exact = new Uint8Array(out_bin.byteLength);
        exact.set(out_bin);
        const outputFile = new Blob([exact.buffer], { type: 'application/sla' });
        // send a message back to the main thread with the output file, exit code and output file name
        const successMsg: WorkerSuccessMessage = {
          blob: outputFile,
          outName: outName,
          exitCode: exitCode
        };
        self.postMessage(successMsg);
      } catch (err) {
        const error = err as Error;
        const errorMsg: WorkerErrorMessage = {
          type: 'error',
          message: error.message,
          error: error.stack ?? null
        };
        self.postMessage(errorMsg);
      } finally {
        // Always free staged files so a failed run cannot leave stale entries in the
        // MEMFS that would break later calls reusing the same filename. Guard each
        // unlink: the file may not exist (e.g. output never written on failure).
        if (stagedInput) {
          try { mod.FS_unlink(inName); } catch { /* already gone */ }
        }
        for (const name of stagedExtras) {
          try { mod.FS_unlink(name); } catch { /* already gone */ }
        }
        if (inName !== outName) {
          try { mod.FS_unlink(outName); } catch { /* never created */ }
        }
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
