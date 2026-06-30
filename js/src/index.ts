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
 * BSD-2-Clause build of niimath. Loads the minimal WASM module that ships with
 * the published `@niivue/niimath` package. For the GPL-2 build that additionally
 * provides `-spm_coreg`/`-spm_deface`/`-allineate`/`-deface`, import from
 * `@niivue/niimath/gpl` instead.
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
