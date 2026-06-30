// BSD WASM worker entry point. Loads the BSD-licensed niimath module and wires up
// the shared message handler. The worker is of type "module" so it can use ES6
// import syntax; the Emscripten code is compiled with -s EXPORT_ES6=1 -s
// MODULARIZE=1 so bundlers can import the worker and wasm with code splitting.
import Module from './niimath.js';
import { setupWorker, type EmscriptenModuleFactory } from './workerImpl';

setupWorker(Module as unknown as EmscriptenModuleFactory);
