#!/usr/bin/env bash
# Regenerate the ROMEO parity oracle (romeo_plan.md M0).
#
#   test/romeo_oracle.sh [<ROMEO.jl-checkout>] [<outdir>]
#
# Defaults: checkout /Users/chris/src/ROMEO.jl, outdir test/romeo_ref (gitignored).
# JULIA_NUM_THREADS=1 is pinned: ROMEO's weight loop is `Threads.@threads for dim in 1:3`, and the
# oracle must be reproducible run-to-run.  The generated goldens are large and derived, so they are
# NOT committed; only this script, romeo_oracle.jl and the small scalar manifest under
# test/romeo_golden/ are.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
checkout="${1:-/Users/chris/src/ROMEO.jl}"
outdir="${2:-$here/romeo_ref}"

if ! command -v julia >/dev/null 2>&1; then
  echo "SKIP: julia not found — the ROMEO oracle needs the pinned Julia environment (romeo_plan.md §2.1)" >&2
  exit 77
fi
if [ ! -d "$checkout" ]; then
  echo "SKIP: ROMEO.jl checkout not found at '$checkout'" >&2
  exit 77
fi

mkdir -p "$outdir"
JULIA_NUM_THREADS=1 julia --startup-file=no --color=no \
  "$here/romeo_oracle.jl" "$checkout" "$outdir" "${@:3}"

echo "oracle written to $outdir"
