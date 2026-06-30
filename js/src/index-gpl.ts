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
 * GPL-2 build of niimath. Identical API to the BSD `Niimath`, but backed by the
 * WASM module compiled with `GPL=1` — which additionally includes the GPL-2
 * SPM coregistration operations (`-spm_coreg`, `-spm_deface`). The BSD allineate
 * registration/defacing (`-allineate`, `-deface`) is in both builds.
 *
 * **Licensing:** importing from `@niivue/niimath/gpl` pulls in GPL-2 code, so any
 * work that bundles this module becomes a GPL-2 combined work. Use the default
 * `@niivue/niimath` import if you need to stay BSD-2-Clause (it still provides
 * `-allineate`/`-deface`, just not the SPM operations).
 */
export class Niimath extends NiimathBase {
  constructor() {
    // The `new Worker(new URL(...))` literal must stay here so esbuild can
    // statically discover and bundle worker-gpl.js (and its niimath-gpl.wasm)
    // into a separate chunk.
    super(operators as Operators, () =>
      new Worker(new URL('./worker-gpl.js', import.meta.url), { type: 'module' })
    );
  }
}
