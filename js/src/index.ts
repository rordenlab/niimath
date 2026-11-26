import operators from './niimathOperators.json' with { type: 'json' };
import { NiimathBase, dataTypes, type Operators } from './core';

export { dataTypes } from './core';
export type {
  Operators,
  OperatorDefinition,
  ImageProcessorMethods,
  MeshOptions,
  BitmapOptions,
  DataType,
  ImageSource,
  QcTissues
} from './core';

/**
 * niimath (BSD-2-Clause). Loads the WASM module shipped with the `@niivue/niimath`
 * package, providing the full BSD feature set including `-allineate`/`-deface`.
 * (The optional GPL `-spm_coreg`/`-spm_deface` are no longer published; the
 * permissive allineate engine supersedes them.)
 */
export class Niimath extends NiimathBase {
  constructor() {
    // The `new Worker(new URL(...))` literals must stay here so bundlers can
    // statically discover both workers. Threads need SharedArrayBuffer, which only
    // cross-origin-isolated pages (COOP/COEP headers) have; others get the serial build.
    super(operators as Operators, () =>
      globalThis.crossOriginIsolated
        ? new Worker(new URL('./worker-mt.js', import.meta.url), { type: 'module' })
        : new Worker(new URL('./worker.js', import.meta.url), { type: 'module' })
    );
  }
}
