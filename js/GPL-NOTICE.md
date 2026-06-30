# GPL-2 notice and written offer for source (`@niivue/niimath/gpl`)

This file applies **only** to the GPL build of this package — the files
`dist/index-gpl.js`, `dist/niimath-gpl.js`, and `dist/niimath-gpl.wasm`, reached
via `import { Niimath } from '@niivue/niimath/gpl'`. The default
`@niivue/niimath` build (`dist/niimath.js` / `dist/niimath.wasm`) contains **no
GPL code** and is BSD-2-Clause; see `LICENSE`.

## What the GPL build is

`dist/niimath-gpl.wasm` is a **combined work**: the BSD-2-Clause niimath sources
compiled together with the optional GPL-2 SPM coregistration module
(`-spm_coreg` / `-spm_deface`). Because BSD-2-Clause is GPL-compatible, the
combined binary as a whole is licensed under the **GNU General Public License,
version 2** (the full text is in `LICENSE.GPL-2.0.txt`). Bundling or
redistributing this build, or a work that incorporates it, makes that work a
GPL-2 work and obliges you to honor GPL-2 §3 (source availability) downstream.

The GPL module is derived from SPM:

> SPM is free software; you can redistribute it and/or modify it under the terms
> of the GNU General Public Licence as published by the Free Software
> Foundation; either version 2 of the Licence, or (at your option) any later
> version. http://www.gnu.org/licenses/gpl-2.0.html

## Complete corresponding source

The complete corresponding source for the GPL build is **included in this
package**, in the `corresponding-source/` directory (GPL-2 §3(a)). It contains
every source file linked into `niimath-gpl.wasm` — the BSD niimath C sources and
the GPL `spm_coreg` sources (`corresponding-source/src/GPL/`) — plus the scripts
used to control its compilation (`src/Makefile`, `src/CMakeLists.txt`, and the
`js/` build wiring). The zlib library linked at build time is provided by the
Emscripten toolchain and is omitted under GPL-2's System Library exception. See
`corresponding-source/README.SOURCE.md` for rebuild instructions.

The same source is also maintained publicly at:

- **niimath (BSD sources + build scripts):** https://github.com/rordenlab/niimath
- **niimath_gpl (GPL `spm_coreg` sources), pinned submodule commit:**
  https://github.com/rordenlab/niimath_gpl — commit
  `d589203cb9a0b0bf8899aa365740c85ba00c825e`
  (the exact `src/GPL` submodule revision used to build the published
  `niimath-gpl.wasm`).

To reproduce the binary from the bundled source:

```bash
cd corresponding-source/js
bun install
bun run makeWasmGpl   # emits src/niimath-gpl.{js,wasm}; needs emscripten (emcc)
bun run build         # bundles dist/index-gpl.js + dist/niimath-gpl.*
```

## Written offer (GPL-2 §3(b))

The complete corresponding source is shipped with the binary under §3(a) (see
above), so no separate written offer is required. As a convenience, the same
source is also kept publicly at the repositories listed above.

## SPDX

```
SPDX-License-Identifier: GPL-2.0-only
```

(The package as a whole — both builds — is `BSD-2-Clause AND GPL-2.0-only`.)
