// Threaded (OpenMP) worker, chosen by index.ts on cross-origin-isolated pages.
import Module from './niimath-mt.js';
import { setupWorker, type EmscriptenModuleFactory } from './workerImpl';

setupWorker(Module as unknown as EmscriptenModuleFactory);
