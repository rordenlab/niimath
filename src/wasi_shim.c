/* wasi_shim.c — reactor/lifecycle bridge for the zlib-free WASI-C niimath backend.
 *
 * This file contains ONLY the reactor entry points. It performs NO algorithm work and touches
 * none of the shared allineate sources. All file I/O flows through wasi-libc stdio, which the
 * host services via WASI Preview 1 imports backed by a host-side in-memory VFS (js/src/wasi).
 *
 * Design notes:
 *   - The niimath translation unit is compiled with -Dmain=niimath_entry (see the Makefile
 *     wasm-wasi target), renaming its `int main(argc,argv)` to `int niimath_entry(argc,argv)`
 *     without editing niimath.c. `main32`/`main64`/etc. are distinct tokens and unaffected.
 *   - The reactor exports memory, _initialize (emitted by -mexec-model=reactor), malloc/free,
 *     and nii_run. Input/output NIfTI + text files live in the host VFS; niimath fopen/fwrite
 *     reach them through path_open/fd_* imports. The JS adapter stages inputs and reads outputs
 *     directly from that host table (the production zero-copy path), so nii_add_file /
 *     nii_read_file are provided at the JS adapter layer, not as raw WASM exports.
 *   - nii_run takes a pointer to a NUL-separated packed argv buffer + its byte length, splits it
 *     into argv[], and invokes niimath_entry. A trailing NUL terminates the last token.
 *   - exit() on an error path becomes proc_exit -> a typed host exception; the adapter discards
 *     the instance. A successful run returns normally from niimath_entry, so the instance can be
 *     reused after the host resets the VFS.
 */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifndef NII_WASI_MAX_ARGS
#define NII_WASI_MAX_ARGS 256
#endif

extern int niimath_entry(int argc, char **argv);

/* The JS side writes the packed argv buffer using the reactor's exported malloc/free
 * (-Wl,--export=malloc,free) — no dedicated allocator wrappers here (they would collide with
 * core.c's nii_malloc datatype-buffer helper). */

/* Run niimath with an argv reconstructed from a packed NUL-separated buffer.
 * argvPtr points to argvLen bytes: "niimath\0in.nii\0-add\01\0out.nii\0" (a leading argv[0] of
 * "niimath" is expected; the JS adapter supplies it). Returns niimath's exit code, or -1 on a
 * malformed buffer (too many tokens / empty). Does not free argvPtr — the caller owns it. */
__attribute__((export_name("nii_run")))
int nii_run(const char *argvPtr, int argvLen) {
    if (!argvPtr || argvLen <= 0) return -1;
    static char *argv[NII_WASI_MAX_ARGS];
    int argc = 0;
    int i = 0;
    while (i < argvLen && argc < NII_WASI_MAX_ARGS) {
        /* start of a token */
        const char *tok = argvPtr + i;
        /* advance to the NUL that ends this token */
        while (i < argvLen && argvPtr[i] != '\0') i++;
        if (i >= argvLen) {
            /* unterminated final token: only valid if it is non-empty and we treat the
             * buffer end as an implicit terminator — but our contract requires a trailing
             * NUL, so a missing one is malformed. */
            if (argvPtr[argvLen - 1] != '\0') return -1;
        }
        argv[argc++] = (char *)tok;
        i++; /* skip the NUL */
    }
    if (argc == 0) return -1;
    if (argc >= NII_WASI_MAX_ARGS) return -1;
    argv[argc] = NULL;
    /* Flush any buffered stdout/stderr from a prior run before starting a fresh one. */
    fflush(stdout);
    fflush(stderr);
    int rc = niimath_entry(argc, argv);
    fflush(stdout);
    fflush(stderr);
    return rc;
}
