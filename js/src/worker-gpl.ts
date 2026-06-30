// GPL WASM worker entry point. Identical to worker.ts except it loads the GPL
// build of niimath (./niimath-gpl.js), which additionally includes the GPL-2
// spm_coreg module (-spm_coreg/-spm_deface) and allineate (-allineate/-deface).
// A binary built from this module is a GPL-2 combined work — it is exposed via a
// separate package.json export ("@niivue/niimath/gpl") so consumers explicitly
// opt into the GPL licensing.
import Module from './niimath-gpl.js';
import { setupWorker, type EmscriptenModuleFactory } from './workerImpl';

setupWorker(Module as unknown as EmscriptenModuleFactory);
