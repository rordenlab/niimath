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
ALOBJS= allineate.o powell_newuoa.o coreg_fast.o

#Zlib flags
ZFLAGS= -DHAVE_ZLIB -lz
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
#dtifit's dependencies for the WASM target (which otherwise ships no tensor support);
#grouped so "DTIFIT=0 make wasm" restores the original tensor-free WASM build
DT_WASM= $(TFILES) $(DTFLAGS)
#QC metrics for the WASM target; grouped so "QC=0 make wasm" drops it (qc.c is
#self-contained and wasm-safe: only core.h helpers, comparator-free quickselect, no OpenMP)
QC_WASM= $(QCFLAGS)
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
    WASM_AL= $(AFLAGS) al_wasm.o pn_wasm.o cf_wasm.o
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

all: $(ALOBJS)
	$(CNAME) $(CFLAGS) $(BFLAGS) $(TFILES) $(DTFLAGS) $(QCFLAGS) $(GPL_SRCS) $(GFILES) $(AFLAGS) $(ALOBJS) $(MFLAGS) $(UFLAGS) $(UFLAGS64) $(ZFLAGS) $(ZSTDFLAGS) $(ZSTDLIB) $(OMPLINK) $(LDFLAGS) $(BMPFLAGS) -flto -o niimath

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

# Auto-generated header dependencies for the prebuilt registration objects (from -MMD).
-include allineate.d powell_newuoa.d coreg_fast.d

#allineate source files for non-default targets (compiled directly, not as separate .o)
ALSRC=
ifneq "$(AL)" "0"
ALSRC= allineate.c powell_newuoa.c coreg_fast.c
endif

#issue30: static build
static:
	gcc -O3 -static -std=gnu99 $(BFLAGS) $(TFILES) $(DTFLAGS) $(QCFLAGS) $(GPL_SRCS) $(GFILES) $(MFLAGS) $(AFLAGS) $(AFASTMATH) $(OMPFLAGS) $(ALSRC) $(UFLAGS) $(UFLAGS64) $(ZFLAGS) $(ZSTDFLAGS) $(ZSTDLIB) $(OMPLINK) $(LDFLAGS) -flto -o niimath

# tiny: terminal executable to emulate WASM (no OpenMP)
tiny:
	$(CNAME) -O3 -ffunction-sections -fdata-sections $(AFASTMATH) $(MFLAGS) $(UFLAGS) $(LDFLAGS) $(ZFLAGS) $(BMPFLAGS) -o niimath

# tiny without nii2mesh (no OpenMP)
nano:
	$(CNAME) -O3 -ffunction-sections -fdata-sections $(AFASTMATH) $(UFLAGS) $(LDFLAGS) -o niimath

debug: CFLAGS   = -O0 -ffunction-sections -fdata-sections
debug:
	$(CNAME) $(CFLAGS) $(TFILES) $(DTFLAGS) $(QCFLAGS) $(GPL_SRCS) $(MFLAGS) $(AFLAGS) $(OMPFLAGS) $(ALSRC) $(BMPFLAGS) $(UFLAGS) $(ZFLAGS) $(ZSTDFLAGS) $(ZSTDLIB) $(OMPLINK) $(LDFLAGS) -o niimath

verbose:
	$(CNAME) -O0 -Wall -Wextra -Wno-sign-compare $(TFILES) $(DTFLAGS) $(QCFLAGS) $(GPL_SRCS) $(MFLAGS) $(AFLAGS) $(OMPFLAGS) $(ALSRC) $(UFLAGS) $(ZFLAGS) $(ZSTDFLAGS) $(ZSTDLIB) $(OMPLINK) $(LDFLAGS) -o niimath

# Lightweight undefined-behavior diagnostics. Unlike the Apple-Clang ASan target below,
# UBSan is safe with Homebrew libomp on Apple Silicon and retains the normal OpenMP paths.
ubsan:
	$(CNAME) -O1 -g -Wno-deprecated -fsanitize=undefined -fno-sanitize-recover=undefined -fno-omit-frame-pointer $(TFILES) $(DTFLAGS) $(QCFLAGS) $(GPL_SRCS) $(MFLAGS) $(AFLAGS) $(OMPFLAGS) $(ALSRC) $(UFLAGS) $(ZFLAGS) $(ZSTDFLAGS) $(ZSTDLIB) $(OMPLINK) $(LDFLAGS) -o niimath

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
	$(CNAME) -O1 -g -Wno-deprecated -fsanitize=address -fno-omit-frame-pointer $(TFILES) $(DTFLAGS) $(QCFLAGS) $(GPL_SRCS) $(MFLAGS) $(AFLAGS) $(ALSRC) $(UFLAGS) $(ZFLAGS) $(ZSTDFLAGS) $(ZSTDLIB) $(LDFLAGS) -o niimath

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
wasm_al_objs:
	emcc -O3 -ffast-math -fno-finite-math-only -DEMSCRIPTEN $(AFLAGS) -c allineate.c -o al_wasm.o
	emcc -O3 -ffast-math -fno-finite-math-only -DEMSCRIPTEN -c powell_newuoa.c -o pn_wasm.o
	emcc -O3 -ffast-math -fno-finite-math-only -DEMSCRIPTEN $(AFLAGS) -c coreg_fast.c -o cf_wasm.o

# WASM_OUT names the emitted ES6 module (emcc derives the sibling .wasm from it).
# Default is the BSD ../js/src/niimath.js; the GPL build overrides it to
# niimath-gpl.js (see the makeWasmGpl npm script) so the two binaries can ship side
# by side under separate package.json exports.
WASM_OUT= ../js/src/niimath.js

wasm: $(WASM_AL_PREREQ)
	emcc -O3 -ffast-math -fno-finite-math-only -DEMSCRIPTEN $(MFLAGS) $(DT_WASM) $(QC_WASM) $(WASM_AL) $(GPL_SRCS) $(UFLAGS) $(ZFLAGS) $(BMPFLAGS) -s USE_ZLIB=1 -s STACK_SIZE=4194304 -s EXPORTED_RUNTIME_METHODS='["callMain", "ccall", "cwrap", "FS_createDataFile", "FS_readFile", "FS_unlink", "allocateUTF8", "getValue", "stringToUTF8", "setValue"]' -s ALLOW_MEMORY_GROWTH=1 -s WASM=1 -s EXPORT_ES6=1 -s MODULARIZE=1 -s EXPORTED_FUNCTIONS='["_main", "_malloc", "_free"]' -s INVOKE_RUN=0 -o $(WASM_OUT)
	# hint: consider further optimizations:
	#   wasm-opt -O3 -o output.wasm niimath.wasm; rm niimath.wasm; mv output.wasm niimath.wasm

#WASM callMain details https://youtu.be/c8hZFtl8EuQ?si=NpAMjf7ka5XtsiYw
