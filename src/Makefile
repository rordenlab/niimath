# Regular use, set compiler name, compiler flags, openmp flags zlib flags
CNAME=gcc
#universal flags always included
UFLAGS= niimath.c core.c core32.c nifti_io.c -DFSLSTYLE -DPIGZ -DREJECT_COMPLEX -lm
#experimental conform feature
UFLAGS+= conform.c -DHAVE_CONFORM
#Bias field correction (from AFNI 3dUnifize, public domain)
UFLAGS+= unifize.c
UFLAGS64= -DHAVE_64BITS core64.c

#OpenMP flags (platform-specific); used by both core ops and allineate
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
OMPFLAGS= -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include
OMPLINK= -L/opt/homebrew/opt/libomp/lib -lomp
else
OMPFLAGS= -fopenmp
OMPLINK= -fopenmp
endif

#Affine registration (from AFNI 3dAllineate, public domain)
AFASTMATH= -ffast-math -fno-finite-math-only
AFLAGS= -DHAVE_ALLINEATE
ALOBJS= allineate.o powell_newuoa.o coreg_fast.o reface.o

#ROMEO phase unwrapping (-romeo). Faithful MIT-licensed port of ROMEO.jl + the
#MriResearchTools.jl helpers its CLI app uses; see src/romeo.LICENSE. Enabled by default;
#run "ROMEO=0 make" to omit it. romeo.o is a PREBUILT OBJECT compiled STRICT-FP: the rest of
#niimath is whole-program -ffast-math, and a single reassociated/contracted expression can move
#an 8-bit edge weight across a rescale() bin boundary, change the spanning tree, and shift a
#whole connected region by 2*pi. These flags must come AFTER $(CFLAGS) so they win.
ROMEOFLAGS= -DHAVE_ROMEO
ROMEOOBJS= romeo.o
ROMEO_STRICT_FP= -fno-fast-math -ffp-contract=off

#Zlib flags
ZFLAGS= -DHAVE_ZLIB -lz
# Optional zlib-ng for native builds WITHOUT CMake. Point ZLIBNG_ROOT at a ZLIB_COMPAT build
# (drop-in libz.a + zlib.h). Build one using zlib-ng's own autotools configure — no CMake needed:
#   git clone https://github.com/zlib-ng/zlib-ng && cd zlib-ng && ./configure --zlib-compat && make
#   make -C src ZLIBNG_ROOT=/abs/path/to/zlib-ng
# NOTE: auto-detecting an *installed* zlib-ng does NOT work — the packaged builds (Homebrew, distros)
# are non-compat (zng_-prefixed zlib-ng.h), so they are not drop-in for niimath's + spng.c's <zlib.h>
# API. A drop-in requires the ZLIB_COMPAT build, which is why the CMake release path builds it itself.
ifdef ZLIBNG_ROOT
ZNG_LIB := $(firstword $(wildcard $(ZLIBNG_ROOT)/libz.a $(ZLIBNG_ROOT)/lib/libz.a))
ZNG_INC := $(firstword $(wildcard $(ZLIBNG_ROOT)/zlib.h $(ZLIBNG_ROOT)/include/zlib.h))
ifeq ($(ZNG_LIB),)
$(error ZLIBNG_ROOT set but no libz.a found under $(ZLIBNG_ROOT) or $(ZLIBNG_ROOT)/lib -- run ./configure --zlib-compat && make in the zlib-ng source)
endif
ifeq ($(ZNG_INC),)
$(error ZLIBNG_ROOT set but no zlib.h found under $(ZLIBNG_ROOT) or $(ZLIBNG_ROOT)/include -- use the same ZLIB_COMPAT build for headers and libz.a)
endif
ZFLAGS= -DHAVE_ZLIB -I$(dir $(ZNG_INC)) $(ZNG_LIB)
endif
# Emscripten zlib for the `wasm` target: defaults to emscripten's ported (madler) zlib via
# -sUSE_ZLIB. Override to link a custom static wasm zlib — e.g. zlib-ng (ZLIB_COMPAT), ~1.5x
# faster on the .nii.gz path — by pre-building it (emcmake cmake -DZLIB_COMPAT=ON ...) and:
#   make wasm WASM_ZLIB="-DHAVE_ZLIB -I<zng-wasm-build> <zng-wasm-build>/libz.a"
WASM_ZLIB= -DHAVE_ZLIB -lz -s USE_ZLIB=1
#Zstd flags (auto-detected; run "ZSTD=0 make" to disable)
ZSTD_FOUND := $(shell pkg-config --exists libzstd 2>/dev/null && echo 1 || (echo 'int main(){return 0;}' | $(CNAME) -xc - -lzstd -o /dev/null 2>/dev/null && echo 1 || echo 0))
ifeq "$(ZSTD_FOUND)" "1"
ZSTDFLAGS= -DHAVE_ZSTD $(shell pkg-config --cflags libzstd 2>/dev/null)
ZSTDLIB := $(shell pkg-config --libs --static libzstd 2>/dev/null || echo "-lzstd")
else
ZSTDFLAGS=
ZSTDLIB=
endif
#compiler flags
CFLAGS   = -O3 -ffunction-sections -fdata-sections
LDFLAGS += -Wl,--gc-sections
# CFLAGS= -O0 -Wall
# GIfTI, obj, vtk, stl, etc support
GFILES= -DHAVE_FORMATS base64.c
#Butterworth bandpass flags
BFLAGS= -DHAVE_BUTTERWORTH bw.c
#Tensor decomposition
TFILES= -DHAVE_TENSOR tensor.c
#Diffusion tensor fit (linear, emulates FSL dtifit; fit math from AFNI 3dDWItoDT, public domain; needs TFILES)
DTFLAGS= -DHAVE_DTIFIT dtifit.c
#Anatomical QC metrics (MRIQC-style IQMs from a T1 + segmentation; self-contained)
QCFLAGS= -DHAVE_QC qc.c

#MEDIC multi-echo distortion correction (--medic, -unwarp). Clean-room emulation of the workflow
#in Van et al. 2026 (doi:10.1162/IMAG.a.1262); see test/medic_reference_manifest.md in the medic_bench repository. REQUIRES
#ROMEO (it calls romeo_unwrap_frame), so "ROMEO=0" implies "MEDIC=0". Ordinary FP: unlike
#romeo.c, nothing here is reassociation-sensitive, so medic.c rides the normal source line.
MEDICFLAGS= -DHAVE_MEDIC medic.c
#dtifit's dependencies for the WASM target (which otherwise ships no tensor support);
#grouped so "DTIFIT=0 make wasm" restores the original tensor-free WASM build
DT_WASM= $(TFILES) $(DTFLAGS)
#QC metrics for the WASM target; grouped so "QC=0 make wasm" drops it (qc.c is
#self-contained and wasm-safe: only core.h helpers, comparator-free quickselect, no OpenMP)
QC_WASM= $(QCFLAGS)
# MEDIC in the default Emscripten build: --medic and -unwarp are useful in the browser and
# -unwarp is cheap, so they ship. DELIBERATE LIMIT, not an oversight: wasm32 is a 4 GiB address
# space and every wasm target sets -DFORCE_INT32_MAX, while --medic holds the whole run resident
# (n3 * T * (2*echoes + 3) * 4 bytes). The 170-frame two-echo demo needs ~1.9 GB, so long
# multi-echo runs will not fit in the browser -- estimate on a native binary and apply -unwarp in
# wasm. Omitted from tiny/nano (like every other optional feature) and from the WASI reactor,
# which is a lean zlib-free compute backend.
MEDIC_WASM= $(MEDICFLAGS)
#Optional GPL spm_coreg module (GNU GPL-2): "GPL=1 make" enables -spm_coreg. OFF by
#default so the repository builds remain BSD-2. The GPL sources live in src/GPL,
#provided by the niimath_gpl submodule (a binary built with GPL=1 is a GPL-2 work).
GPL_SRCS=
ifeq "$(GPL)" "1"
    ifeq ($(wildcard GPL/spm_coreg.c),)
        $(error GPL=1 set but src/GPL is empty -- initialize the GPL submodule: git submodule update --init src/GPL)
    endif
    # GPL estimate only (returns the rigid transform); reslice/interpolation use BSD code.
    GPL_SRCS= -DHAVE_GPL -IGPL -I. GPL/cost.c GPL/hist2.c GPL/loaduint8.c GPL/matrix.c GPL/powell.c GPL/smooth.c GPL/spm_coreg.c GPL/spmcoreg_niimath.c
endif
#WASM allineate inclusion. allineate (-allineate/-deface) is BSD/public-domain
#(adapted from AFNI 3dAllineate), so it ships in BOTH the default BSD wasm and the
#GPL wasm. The GPL build additionally REQUIRES it because the GPL reslice calls the
#BSD nii_reslice_affine in allineate.c. Link allineate/powell as pre-built objects
#(wasm_al_objs rule) so they get -ffast-math like the native build; -DHAVE_ALLINEATE
#flows to the other inline TUs (coreFLT deface dispatch, GPL reslice hook). Set
#AL=0 to omit it (e.g. for a smaller minimal wasm).
WASM_AL=
WASM_AL_PREREQ=
ifneq "$(AL)" "0"
    WASM_AL= $(filter-out -DHAVE_QWARP,$(AFLAGS)) al_wasm.o pn_wasm.o cf_wasm.o rf_wasm.o
    WASM_AL_PREREQ= wasm_al_objs
endif

#Marching Cubes flag: run "MC=1 make" for new MarchingCubes instead of classic oldcubes
MFLAGS= -DNII2MESH meshify.c quadric.c bwlabel.c radixsort.c fdr.c -DUSE_CLASSIC_CUBES oldcubes.c

#run "MC=1 make" for new MarchingCubes algorithm (handles ambiguities)
ifeq "$(MC)" "1"
	MFLAGS= -DNII2MESH meshify.c quadric.c bwlabel.c radixsort.c fdr.c MarchingCubes.c
endif
#Bitmap ceation: either stb_image or spng, the latter creates more compact PNGs, requires ZFLAGS
BMPFLAGS= -DHAVE_BMP filter.c bmp.c spng.c

#run "STB=1 make" for STB instead of libspng bitmaps
ifeq "$(STB)" "1"
	BMPFLAGS= -DHAVE_BMP filter.c bmp.c  -DHAVE_STB
endif
ifeq ($(CXX),g++)
	CFLAGS += -std=gnu99
endif

#OpenMP enabled by default; run "OMP=0 make -j" to disable
ifeq "$(OMP)" "0"
	# clear both so recipes that reference $(OMPLINK) explicitly do not link libomp
	OMPFLAGS=
	OMPLINK=
else
	CFLAGS+= $(OMPFLAGS)
	LDFLAGS+= $(OMPLINK)
endif

# Whole-program fast-math (release builds). The standalone allineate project builds ALL of
# its sources with -ffast-math, and the registration calls shared nifti_io matrix helpers
# (nifti_mat44_inverse / nifti_mat44_mul) DURING the fit. Compiling nifti_io strict-FP while
# allineate.o is fast-FP made those 4x4 ops round differently, so on the flat cross-modal cost
# surface NEWUOA settled at a slightly different optimum and niimath's -allineate diverged from
# the allineate standard. Matching the fast-math build makes the registration bit-identical.
# (The `debug` target re-assigns CFLAGS, so it stays strict.)
CFLAGS+= $(AFASTMATH)

#run "MESH=0 make" for minimal compile without nii2mesh functions
ifeq "$(MESH)" "0"
	MFLAGS=
endif

#run "AL=0 make" for compile without allineate registration
ifeq "$(AL)" "0"
	AFLAGS=
	ALOBJS=
endif

#run "ROMEO=0 make" for compile without ROMEO phase unwrapping (-romeo)
ifeq "$(ROMEO)" "0"
	ROMEOFLAGS=
	ROMEOOBJS=
	ROMEO_SRC=
	# MEDIC calls romeo_unwrap_frame(); without ROMEO there is nothing to call.
	MEDICFLAGS=
else
	ROMEO_SRC= romeo.c
endif

#run "MEDIC=0 make" for compile without MEDIC distortion correction (--medic, -unwarp)
ifeq "$(MEDIC)" "0"
	MEDICFLAGS=
endif

# run "QWARP=1 make" to include -qwarp nonlinear (deformable) registration, an attributed
# public-domain port of AFNI 3dQwarp (qwarp.c/.h). OFF by default: it is memory/CPU-heavy,
# benefits greatly from OpenMP, and is impractically slow in WASM. Requires allineate
# (uses powell_newuoa.o), so do not combine with AL=0.
ifeq "$(QWARP)" "1"
ifeq "$(AL)" "0"
$(error QWARP=1 requires allineate (qwarp uses powell_newuoa); do not combine with AL=0)
endif
ifneq "$(filter wasm wasm-wasi wasm-emcc-core tiny nano,$(MAKECMDGOALS))" ""
$(error QWARP=1 is native-only (qwarp benefits greatly from OpenMP and is impractically slow in serial/WASM builds); build a full native OpenMP target)
endif
	AFLAGS += -DHAVE_QWARP
	ALOBJS += qwarp.o
	QWARP_SRC= qwarp.c
endif

#run "ZSTD=0 make" for compile without zstd support
ifeq "$(ZSTD)" "0"
	ZSTDFLAGS=
	ZSTDLIB=
endif

#run "DTIFIT=0 make" for compile without diffusion tensor fit
ifeq "$(DTIFIT)" "0"
	DTFLAGS=
	DT_WASM=
endif

#run "QC=0 make" for compile without anatomical QC metrics (--qc)
ifeq "$(QC)" "0"
	QCFLAGS=
endif

ifeq ($(UNAME_S),Darwin)
	# macOS linker does not support --gc-sections
	LDFLAGS =
endif

all: $(ALOBJS) $(ROMEOOBJS)
	$(CNAME) $(CFLAGS) $(BFLAGS) $(TFILES) $(DTFLAGS) $(QCFLAGS) $(MEDICFLAGS) $(GPL_SRCS) $(GFILES) $(AFLAGS) $(ALOBJS) $(ROMEOFLAGS) $(ROMEOOBJS) $(MFLAGS) $(UFLAGS) $(UFLAGS64) $(ZFLAGS) $(ZSTDFLAGS) $(ZSTDLIB) $(OMPLINK) $(LDFLAGS) $(BMPFLAGS) -flto -o niimath

# Canonical local regression: build, then run the stdlib-only release smoke test
# (`.github/scripts/release_smoke.py`) against the freshly built binary. It synthesizes small
# NIfTI fixtures and exercises the packaging-sensitive paths — gzip read/write, -conform/-gz 0/
# -odt char, feature dispatch, optional zstd — with no external test data. Deeper coverage lives
# elsewhere: the corner-case/conformance suite in the separate `niimath_tests` repo, and the WASI
# feature-parity suite via `bun run test:wasi` in ../js (needs Zig).
.PHONY: test
test: all
	python3 ../.github/scripts/release_smoke.py ./niimath

# The prebuilt registration objects are compiled with flags that depend on command-line
# toggles (OMP=0, GPL=1, CNAME, feature flags), but make keys rebuilds on source timestamps
# only. Without a flag dependency, `make` then `make OMP=0` reuses the OpenMP-compiled objects
# and the link FAILS (unresolved __kmpc_*/omp_* in allineate.o). Fix: a flag-signature stamp.
# `.al_obj_flags` is regenerated every build but only its CONTENTS change (and thus its mtime)
# when the signature changes, so the three .o rebuild exactly when — and only when — the
# object-affecting flags change. No `make -B` needed.
AL_OBJ_SIG := $(CNAME)|$(CFLAGS)|$(AFLAGS)
.PHONY: force_al_flags
force_al_flags: ;
.al_obj_flags: force_al_flags
	@printf '%s' '$(AL_OBJ_SIG)' | cmp -s - $@ 2>/dev/null || printf '%s' '$(AL_OBJ_SIG)' > $@

# -MMD -MP emits a .d file listing every header each object actually includes (nifti_io.h,
# mat44/nifti_image definitions, core32.h, ...), so a changed struct ABI rebuilds the object
# — the manual `.h` prerequisites below only cover the direct headers. The .d files are
# -included at the bottom of this Makefile. Combined with the .al_obj_flags stamp (flag
# changes) this makes the prebuilt objects rebuild on both source/header AND flag changes.
allineate.o: allineate.c allineate.h .al_obj_flags
	$(CNAME) $(CFLAGS) -MMD -MP -c allineate.c -o allineate.o

powell_newuoa.o: powell_newuoa.c .al_obj_flags
	$(CNAME) $(CFLAGS) -MMD -MP -c powell_newuoa.c -o powell_newuoa.o

# Fast SPM/FLIRT-inspired coreg engine (-allineate ... -cost fast/fastcr). Calls
# nifti_smooth_gauss_f32 (generated by core32.c) and the allineate API.
coreg_fast.o: coreg_fast.c coreg_fast.h core32.h allineate.h .al_obj_flags
	$(CNAME) $(CFLAGS) -MMD -MP -c coreg_fast.c -o coreg_fast.o

reface.o: reface.c reface.h core32.h nifti_io.h .al_obj_flags
	$(CNAME) $(CFLAGS) -MMD -MP -c reface.c -o reface.o

# qwarp.o is only linked when QWARP=1 (added to ALOBJS above); the rule is always present.
qwarp.o: qwarp.c qwarp.h core32.h nifti_io.h .al_obj_flags
	$(CNAME) $(CFLAGS) -MMD -MP -c qwarp.c -o qwarp.o

# NOTE: the non-default targets below (debug/verbose/ubsan/sanitize/static/native-noomp) each
# build their OWN romeo_*.o rather than romeo.o. Writing to romeo.o would clobber the release
# object with e.g. UBSan-instrumented code, and `make all` would then find it NEWER than romeo.c
# with an unchanged flag stamp, skip the rebuild, and fail to link on __ubsan_handle_*.
#
# romeo.o has its OWN flag signature: it must rebuild when CNAME, the optimisation/debug/
# sanitizer flags, OMP or ROMEO change even in an AL=0 build (where .al_obj_flags is unused).
# ROMEO_STRICT_FP is appended AFTER $(CFLAGS) so -fno-fast-math wins over the repository-wide
# -ffast-math; that is the whole point of the separate object.
ROMEO_OBJ_SIG := $(CNAME)|$(CFLAGS)|$(ROMEOFLAGS)|$(ROMEO_STRICT_FP)
.PHONY: force_romeo_flags
force_romeo_flags: ;
.romeo_obj_flags: force_romeo_flags
	@printf '%s' '$(ROMEO_OBJ_SIG)' | cmp -s - $@ 2>/dev/null || printf '%s' '$(ROMEO_OBJ_SIG)' > $@

romeo.o: romeo.c romeo.h core.h nifti_io.h .romeo_obj_flags
	$(CNAME) $(CFLAGS) $(ROMEOFLAGS) $(ROMEO_STRICT_FP) -MMD -MP -c romeo.c -o romeo.o

# Auto-generated header dependencies for the prebuilt registration objects (from -MMD).
-include allineate.d powell_newuoa.d coreg_fast.d reface.d qwarp.d romeo.d

#allineate source files for non-default targets (compiled directly, not as separate .o)
ALSRC=
ifneq "$(AL)" "0"
# qwarp.c ($(QWARP_SRC)) is appended when QWARP=1, matching ALOBJS's qwarp.o — both the object
# and direct-source inventories derive from the same QWARP conditional so QWARP=1 links on every target.
ALSRC= allineate.c powell_newuoa.c coreg_fast.c reface.c $(QWARP_SRC)
endif

# True single-threaded release build for native performance comparisons. This compiles the same
# feature inventory and fast-math policy as `all`, but omits all OpenMP compile/link settings and
# emits a separate artifact so benchmarking it never overwrites the normal OpenMP binary.
.PHONY: native-noomp
native-noomp:
	$(if $(ROMEOOBJS),$(CNAME) $(filter-out $(OMPFLAGS),$(CFLAGS)) $(ROMEOFLAGS) $(ROMEO_STRICT_FP) -c romeo.c -o romeo_noomp.o)
	$(CNAME) $(filter-out $(OMPFLAGS),$(CFLAGS)) $(ROMEOFLAGS) $(ROMEOOBJS:romeo.o=romeo_noomp.o) $(BFLAGS) $(TFILES) $(DTFLAGS) $(QCFLAGS) $(MEDICFLAGS) $(GPL_SRCS) $(GFILES) $(AFLAGS) $(ALSRC) $(MFLAGS) $(UFLAGS) $(UFLAGS64) $(ZFLAGS) $(ZSTDFLAGS) $(ZSTDLIB) $(filter-out $(OMPLINK),$(LDFLAGS)) $(BMPFLAGS) -flto -o niimath-noomp

#issue30: static build
static:
	$(if $(ROMEOOBJS),gcc -O3 -static -std=gnu99 $(ROMEOFLAGS) $(ROMEO_STRICT_FP) $(OMPFLAGS) -c romeo.c -o romeo_static.o)
	gcc -O3 -static -std=gnu99 $(ROMEOFLAGS) $(ROMEOOBJS:romeo.o=romeo_static.o) $(BFLAGS) $(TFILES) $(DTFLAGS) $(QCFLAGS) $(MEDICFLAGS) $(GPL_SRCS) $(GFILES) $(MFLAGS) $(AFLAGS) $(AFASTMATH) $(OMPFLAGS) $(ALSRC) $(UFLAGS) $(UFLAGS64) $(ZFLAGS) $(ZSTDFLAGS) $(ZSTDLIB) $(OMPLINK) $(LDFLAGS) -flto -o niimath

# tiny: terminal executable to emulate WASM (no OpenMP)
tiny:
	$(CNAME) -O3 -ffunction-sections -fdata-sections -DFORCE_INT32_MAX $(AFASTMATH) $(MFLAGS) $(UFLAGS) $(LDFLAGS) $(ZFLAGS) $(BMPFLAGS) -o niimath

# tiny without nii2mesh (no OpenMP)
nano:
	$(CNAME) -O3 -ffunction-sections -fdata-sections $(AFASTMATH) $(UFLAGS) $(LDFLAGS) -o niimath

debug: CFLAGS   = -O0 -ffunction-sections -fdata-sections
debug:
	$(if $(ROMEOOBJS),$(CNAME) $(CFLAGS) $(OMPFLAGS) $(ROMEOFLAGS) $(ROMEO_STRICT_FP) -c romeo.c -o romeo_dbg.o)
	$(CNAME) $(CFLAGS) $(ROMEOFLAGS) $(ROMEOOBJS:romeo.o=romeo_dbg.o) $(BFLAGS) $(TFILES) $(DTFLAGS) $(QCFLAGS) $(MEDICFLAGS) $(GPL_SRCS) $(GFILES) $(MFLAGS) $(AFLAGS) $(OMPFLAGS) $(ALSRC) $(BMPFLAGS) $(UFLAGS) $(UFLAGS64) $(ZFLAGS) $(ZSTDFLAGS) $(ZSTDLIB) $(OMPLINK) $(LDFLAGS) -o niimath

verbose:
	$(if $(ROMEOOBJS),$(CNAME) -O0 -Wall -Wextra -Wno-sign-compare $(OMPFLAGS) $(ROMEOFLAGS) $(ROMEO_STRICT_FP) -c romeo.c -o romeo_verb.o)
	$(CNAME) -O0 -Wall -Wextra -Wno-sign-compare $(ROMEOFLAGS) $(ROMEOOBJS:romeo.o=romeo_verb.o) $(BFLAGS) $(TFILES) $(DTFLAGS) $(QCFLAGS) $(MEDICFLAGS) $(GPL_SRCS) $(GFILES) $(MFLAGS) $(AFLAGS) $(OMPFLAGS) $(ALSRC) $(UFLAGS) $(UFLAGS64) $(ZFLAGS) $(ZSTDFLAGS) $(ZSTDLIB) $(OMPLINK) $(LDFLAGS) $(BMPFLAGS) -o niimath

# Lightweight undefined-behavior diagnostics. Unlike the Apple-Clang ASan target below,
# UBSan is safe with Homebrew libomp on Apple Silicon and retains the normal OpenMP paths.
ubsan:
	$(if $(ROMEOOBJS),$(CNAME) -O1 -g -Wno-deprecated -fsanitize=undefined -fno-sanitize-recover=undefined -fno-omit-frame-pointer $(OMPFLAGS) $(ROMEOFLAGS) $(ROMEO_STRICT_FP) -c romeo.c -o romeo_ubsan.o)
	$(CNAME) -O1 -g -Wno-deprecated -fsanitize=undefined -fno-sanitize-recover=undefined -fno-omit-frame-pointer $(ROMEOFLAGS) $(ROMEOOBJS:romeo.o=romeo_ubsan.o) $(BFLAGS) $(TFILES) $(DTFLAGS) $(QCFLAGS) $(MEDICFLAGS) $(GPL_SRCS) $(GFILES) $(MFLAGS) $(AFLAGS) $(OMPFLAGS) $(ALSRC) $(UFLAGS) $(UFLAGS64) $(ZFLAGS) $(ZSTDFLAGS) $(ZSTDLIB) $(OMPLINK) $(LDFLAGS) $(BMPFLAGS) -o niimath

# sanitize checks memory errors - similar to tiny/wasm
# on MacOS consider: export MallocNanoZone=0
# macOS note: ASan here is very slow on real-size images (allineate's scattered 3D
# interpolation thrashes ASan shadow memory — >100x, not specific to these flags;
# OMP=0 and MallocNanoZone=0 do not help). For memory validation use a small image
# (the 64^3 CI fixture runs clean in <1s) plus `leaks`; full-size ASan is impractical.
# macOS/Apple-Clang note: OpenMP is intentionally OMITTED here. Apple Clang + Homebrew
# libomp + -fsanitize=address deadlocks at OpenMP runtime init (even a trivial parallel
# region hangs; OMP_NUM_THREADS=1 does NOT help — the deadlock is at libomp load). Dropping
# the OpenMP flags makes every `#pragma omp` a no-op, so ASan runs serial with no deadlock.
# LeakSanitizer is unavailable on Apple Silicon macOS: for leaks use a normal (non-ASan)
# build with `MallocNanoZone=0 leaks --atExit -- ./niimath ...`. For ASan+OpenMP together,
# build a separate target with Homebrew LLVM clang (/opt/homebrew/opt/llvm/bin/clang).
# Run with: MallocNanoZone=0 make sanitize
sanitize:
	$(if $(ROMEOOBJS),$(CNAME) -O1 -g -Wno-deprecated -fsanitize=address -fno-omit-frame-pointer $(ROMEOFLAGS) $(ROMEO_STRICT_FP) -c romeo.c -o romeo_asan.o)
	$(CNAME) -O1 -g -Wno-deprecated -fsanitize=address -fno-omit-frame-pointer $(ROMEOFLAGS) $(ROMEOOBJS:romeo.o=romeo_asan.o) $(BFLAGS) $(TFILES) $(DTFLAGS) $(QCFLAGS) $(MEDICFLAGS) $(GPL_SRCS) $(GFILES) $(MFLAGS) $(AFLAGS) $(ALSRC) $(UFLAGS) $(UFLAGS64) $(ZFLAGS) $(ZSTDFLAGS) $(ZSTDLIB) $(LDFLAGS) $(BMPFLAGS) -o niimath

# WASM build. The default (BSD) target is minimal. With GPL=1 it additionally links
# allineate (-deface/-allineate) and the spm_coreg module (-spm_coreg/-spm_deface,
# GPL-2) — both via $(WASM_AL)/$(GPL_SRCS), empty unless GPL=1. allineate.c/
# powell_newuoa.c are pre-built (wasm_al_objs) WITH -ffast-math to match the native
# build: that reassociation is load-bearing for the registration cost reduction —
# WITHOUT it the optimizer converges to a worse local minimum in WASM (verified:
# -deface masked 37% vs 60% native; fast-math restores 60%). As of the whole-program
# fast-math decision the MAIN wasm link line also compiles -ffast-math -fno-finite-math-only
# (matching Make/CMake/notarize), so every wasm TU — core ops, the GPL spm_coreg module, and
# the pre-built allineate objects — shares one FP contract; the pre-built objects remain (the
# link consumes al_wasm.o/pn_wasm.o/cf_wasm.o, not the .c). Every OpenMP construct in those
# files is #ifdef _OPENMP-guarded and _OPENMP is never set here, so the scalar fallbacks compile.
# STACK_SIZE is raised from emcc's 64 KB default because the GPL cost function keeps a
# 256x256 joint histogram (double H[65536] = 512 KB) on the stack (cost.c); the default
# overflows into a RuntimeError. Harmless for non-GPL builds (just reserves more).
# Pre-build allineate + powell with -ffast-math (matches the native allineate.o/
# powell_newuoa.o). Load-bearing: the cost-reduction reassociation it enables is what
# makes the WASM registration converge to the same optimum as native.
# The ROMEO object is compiled STRICT-FP even in WASM: unlike allineate (whose cost reduction
# NEEDS the fast-math reassociation to converge), ROMEO's fidelity depends on the exact bit
# pattern of every weight. Per-object FP policies coexist in one link.
WASM_ROMEO=
WASM_ROMEO_PREREQ=
ifneq "$(ROMEO)" "0"
    WASM_ROMEO= $(ROMEOFLAGS) rm_wasm.o
    WASM_ROMEO_PREREQ= wasm_romeo_obj
endif

wasm_romeo_obj:
	emcc -O3 $(ROMEO_STRICT_FP) -DEMSCRIPTEN -DFORCE_INT32_MAX $(ROMEOFLAGS) -c romeo.c -o rm_wasm.o

wasm_al_objs:
	emcc -O3 -ffast-math -fno-finite-math-only -DEMSCRIPTEN $(AFLAGS) -c allineate.c -o al_wasm.o
	emcc -O3 -ffast-math -fno-finite-math-only -DEMSCRIPTEN -c powell_newuoa.c -o pn_wasm.o
	emcc -O3 -ffast-math -fno-finite-math-only -DEMSCRIPTEN $(AFLAGS) -c coreg_fast.c -o cf_wasm.o
	emcc -O3 -ffast-math -fno-finite-math-only -DEMSCRIPTEN -c reface.c -o rf_wasm.o

# WASM_OUT names the emitted ES6 module (emcc derives the sibling .wasm from it).
# Default is the BSD ../js/src/niimath.js; the GPL build overrides it to
# niimath-gpl.js (see the makeWasmGpl npm script) so the two binaries can ship side
# by side under separate package.json exports.
WASM_OUT= ../js/src/niimath.js

wasm: $(WASM_AL_PREREQ) $(WASM_ROMEO_PREREQ)
	emcc -O3 -ffast-math -fno-finite-math-only -DEMSCRIPTEN -DFORCE_INT32_MAX $(MFLAGS) $(DT_WASM) $(QC_WASM) $(MEDIC_WASM) $(WASM_AL) $(WASM_ROMEO) $(GPL_SRCS) $(UFLAGS) $(BMPFLAGS) $(WASM_ZLIB) -s STACK_SIZE=4194304 -s EXPORTED_RUNTIME_METHODS='["callMain", "ccall", "cwrap", "FS_createDataFile", "FS_readFile", "FS_unlink", "allocateUTF8", "getValue", "stringToUTF8", "setValue"]' -s ALLOW_MEMORY_GROWTH=1 -s WASM=1 -s EXPORT_ES6=1 -s MODULARIZE=1 -s EXPORTED_FUNCTIONS='["_main", "_malloc", "_free"]' -s INVOKE_RUN=0 -o $(WASM_OUT)
	# hint: consider further optimizations:
	#   wasm-opt -O3 -o output.wasm niimath.wasm; rm niimath.wasm; mv output.wasm niimath.wasm

#WASM callMain details https://youtu.be/c8hZFtl8EuQ?si=NpAMjf7ka5XtsiYw

# ============================================================================
# WASI-C reactor backend (zlib-free reactor; import surface pinned in js/src/wasi/import-manifest.json)
# ----------------------------------------------------------------------------
# A lean, zlib-free, single-threaded WASI reactor of the computational BSD feature set
# (core math, allineate/deface, dtifit, QC, conform, butterworth). NO zlib/zstd/mesh/bmp/
# gpl/openmp: NIfTI gzip is done host-side (CompressionStream), so the module only sees raw
# .nii bytes. Built with `zig cc` because Zig bundles a complete wasi-libc sysroot + reactor
# linking — no separate wasi-sdk install is required. All file I/O flows through wasi-libc
# stdio, which the host services via 15 WASI Preview-1 imports over an in-memory VFS
# (js/src/wasi). niimath.c is compiled with -Dmain=niimath_entry so wasi_shim.c can own the
# reactor entry (nii_run) without editing niimath.c. Whole-program -ffast-math matches native
# (load-bearing for the registration golden — see the `wasm:` note above).
WASI_CC= zig cc -target wasm32-wasi
# -msimd128: lets clang auto-vectorize the contiguous Gaussian blur_row loops in coreFLT.c.
# Verified to preserve allineate registration
# quality (corr) and keep -deface byte-stable in conformance; every JS runtime (Bun/Node/Chromium)
# has shipped wasm SIMD since 2021. ~+11% wasm size for a large filter speedup.
WASI_CFLAGS= -O3 -flto -ffast-math -fno-finite-math-only -msimd128 -DFORCE_INT32_MAX -DNII_WASI -DFSLSTYLE -DREJECT_COMPLEX -DHAVE_CONFORM -DHAVE_64BITS -DHAVE_ALLINEATE -DHAVE_TENSOR -DHAVE_DTIFIT -DHAVE_QC -DHAVE_BUTTERWORTH
# Scoped compute sources. niimath.c is compiled separately (with -Dmain=niimath_entry).
WASI_SRCS= core.c core32.c core64.c nifti_io.c conform.c unifize.c tensor.c dtifit.c qc.c bw.c allineate.c powell_newuoa.c coreg_fast.c reface.c wasi_shim.c
# ROMEO is compiled separately (strict FP) and linked as an object; see wasm_romeo_obj.
WASI_ROMEO=
WASI_ROMEO_OBJ=
ifneq "$(ROMEO)" "0"
    WASI_ROMEO= $(ROMEOFLAGS)
    WASI_ROMEO_OBJ= romeo_wasi.o
endif
WASI_EXPORTS= -Wl,--export=nii_run -Wl,--export=malloc -Wl,--export=free -Wl,--export=_initialize -Wl,--export-memory
WASI_OUT= ../js/src/niimath-wasi.wasm

# Verify the WASI toolchain; print actionable install instructions, never auto-download.
wasi-toolcheck:
	@command -v zig >/dev/null 2>&1 || { echo "ERROR: 'zig' not found. The WASI-C backend needs Zig >= 0.16 (bundles wasi-libc)."; echo "  Install: brew install zig   (or https://ziglang.org/download/)"; exit 1; }
	@echo "WASI toolchain OK: zig $$(zig version)"

wasm-wasi: wasi-toolcheck
	$(if $(WASI_ROMEO_OBJ),$(WASI_CC) $(WASI_CFLAGS) $(ROMEO_STRICT_FP) $(ROMEOFLAGS) -c romeo.c -o $(WASI_ROMEO_OBJ))
	$(WASI_CC) $(WASI_CFLAGS) $(WASI_ROMEO) -Dmain=niimath_entry -c niimath.c -o niimath_wasi.o
	$(WASI_CC) $(WASI_CFLAGS) $(WASI_ROMEO) -mexec-model=reactor $(WASI_EXPORTS) -Wl,--gc-sections -Wl,--strip-all niimath_wasi.o $(WASI_SRCS) $(WASI_ROMEO_OBJ) -lm -o $(WASI_OUT)
	@rm -f niimath_wasi.o $(WASI_ROMEO_OBJ)
	@# Optional post-parity size pass (Binaryen). Skipped with a hint if wasm-opt is absent, so a
	@# plain checkout still builds. -all enables the wasm features the reactor uses (bulk-memory,
	@# sign-ext, etc.); re-run conformance after a toolchain bump to re-confirm parity.
	@command -v wasm-opt >/dev/null 2>&1 && { wasm-opt -all -O3 --strip-debug $(WASI_OUT) -o $(WASI_OUT).opt && mv $(WASI_OUT).opt $(WASI_OUT) && echo "wasm-opt applied"; } || echo "wasm-opt not found (brew install binaryen) — shipping un-optimized reactor"
	@echo "built $(WASI_OUT) ($$(wc -c < $(WASI_OUT)) bytes)"
	@[ -f ../js/src/wasi/check-imports.mjs ] && node ../js/src/wasi/check-imports.mjs $(WASI_OUT) ../js/src/wasi/import-manifest.json || echo "(WASI import-surface check skipped: js/src/wasi/check-imports.mjs not present)"

# Feature-matched, zlib-free Emscripten build that RETAINS emscripten glue/MEMFS so the
# benchmark isolates runtime+FS overhead (WASI VFS vs MEMFS) at an identical feature set and
# host-side gzip path. Same scoped sources/defines as wasm-wasi; NO -sUSE_ZLIB (raw bytes only).
EMCC_CORE_OUT= ../js/src/niimath-emcc-core.js
wasm-emcc-core:
	$(if $(ROMEOOBJS),emcc -O3 $(ROMEO_STRICT_FP) -DEMSCRIPTEN -DFORCE_INT32_MAX $(ROMEOFLAGS) -c romeo.c -o romeo_emcc.o)
	emcc -O3 -flto -ffast-math -fno-finite-math-only -msimd128 -DEMSCRIPTEN -DFORCE_INT32_MAX -DFSLSTYLE -DREJECT_COMPLEX -DHAVE_CONFORM -DHAVE_64BITS -DHAVE_ALLINEATE -DHAVE_TENSOR -DHAVE_DTIFIT -DHAVE_QC -DHAVE_BUTTERWORTH $(ROMEOFLAGS) \
	  niimath.c core.c core32.c core64.c nifti_io.c conform.c unifize.c tensor.c dtifit.c qc.c bw.c allineate.c powell_newuoa.c coreg_fast.c reface.c $(if $(ROMEOOBJS),romeo_emcc.o) \
	  -s STACK_SIZE=4194304 -s ALLOW_MEMORY_GROWTH=1 -s WASM=1 -s EXPORT_ES6=1 -s MODULARIZE=1 \
	  -s EXPORTED_RUNTIME_METHODS='["callMain","ccall","cwrap","FS_createDataFile","FS_readFile","FS_unlink","stringToUTF8","getValue","setValue"]' \
	  -s EXPORTED_FUNCTIONS='["_main","_malloc","_free"]' -s INVOKE_RUN=0 -o $(EMCC_CORE_OUT)
	@echo "built $(EMCC_CORE_OUT)"
