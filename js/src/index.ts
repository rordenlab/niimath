import operators from './niimathOperators.json' with { type: 'json' };
import { NiimathBase, dataTypes, type Operators } from './core';

export { dataTypes } from './core';
export type {
  Operators,
  OperatorDefinition,
  ImageProcessorMethods,
  MeshOptions,
  BitmapOptions,
  DataType
} from './core';

/**
 * niimath (BSD-2-Clause). Loads the WASM module shipped with the `@niivue/niimath`
 * package, providing the full BSD feature set including `-allineate`/`-deface`.
 * (The optional GPL `-spm_coreg`/`-spm_deface` are no longer published; the
 * permissive allineate engine supersedes them.)
 */
export class Niimath extends NiimathBase {
  constructor() {
    // The `new Worker(new URL(...))` literal must stay here so esbuild can
    // statically discover and bundle worker.js (and its niimath.wasm) into a
    // separate chunk.
    super(operators as Operators, () =>
      new Worker(new URL('./worker.js', import.meta.url), { type: 'module' })
    );
  }
}
