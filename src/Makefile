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
ALOBJS= allineate.o powell_newuoa.o

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
#dtifit's dependencies for the WASM target (which otherwise ships no tensor support);
#grouped so "DTIFIT=0 make wasm" restores the original tensor-free WASM build
DT_WASM= $(TFILES) $(DTFLAGS)
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
#WASM allineate inclusion. The GPL build REQUIRES allineate (GPL reslice calls the
#BSD nii_reslice_affine in allineate.c), so it is bundled whenever GPL=1. The default
#BSD wasm omits it to keep the published @niivue/niimath package minimal — the full
#toolset (allineate + spm_coreg) ships in the separate GPL @niivue/niimath-gpl build.
WASM_AL=
WASM_AL_PREREQ=
ifeq "$(GPL)" "1"
    # Link allineate/powell as pre-built objects (wasm_al_objs rule) so they get
    # -ffast-math like the native build; -DHAVE_ALLINEATE still flows to the other
    # inline TUs (coreFLT deface dispatch, GPL reslice hook).
    WASM_AL= $(AFLAGS) al_wasm.o pn_wasm.o
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

ifeq ($(UNAME_S),Darwin)
	# macOS linker does not support --gc-sections
	LDFLAGS =
endif

all: $(ALOBJS)
	$(CNAME) $(CFLAGS) $(BFLAGS) $(TFILES) $(DTFLAGS) $(GPL_SRCS) $(GFILES) $(AFLAGS) $(ALOBJS) $(MFLAGS) $(UFLAGS) $(UFLAGS64) $(ZFLAGS) $(ZSTDFLAGS) $(ZSTDLIB) $(OMPLINK) $(LDFLAGS) $(BMPFLAGS) -flto -o niimath

# NOTE: these .o targets depend only on their sources, not on $(OMPFLAGS), so make
# does NOT rebuild them when you toggle OMP=0/OMP=1 on the command line. Switching
# OMP on an existing tree leaves a stale object (e.g. an OpenMP-compiled allineate.o
# fails to link under OMP=0, missing libomp). Force a rebuild when toggling: make -B
# (or rm -f allineate.o powell_newuoa.o). Only the .o-based `all`/`static` targets are
# affected; the inline-$(ALSRC) targets (wasm/debug/sanitize/verbose) recompile anyway.
allineate.o: allineate.c allineate.h
	$(CNAME) $(CFLAGS) $(AFASTMATH) $(OMPFLAGS) -c allineate.c -o allineate.o

powell_newuoa.o: powell_newuoa.c
	$(CNAME) $(CFLAGS) $(AFASTMATH) $(OMPFLAGS) -c powell_newuoa.c -o powell_newuoa.o

#allineate source files for non-default targets (compiled directly, not as separate .o)
ALSRC=
ifneq "$(AL)" "0"
ALSRC= allineate.c powell_newuoa.c
endif

#issue30: static build
static:
	gcc -O3 -static -std=gnu99 $(BFLAGS) $(TFILES) $(DTFLAGS) $(GPL_SRCS) $(GFILES) $(MFLAGS) $(AFLAGS) $(AFASTMATH) $(OMPFLAGS) $(ALSRC) $(UFLAGS) $(UFLAGS64) $(ZFLAGS) $(ZSTDFLAGS) $(ZSTDLIB) $(OMPLINK) $(LDFLAGS) -flto -o niimath

# tiny: terminal executable to emulate WASM (no OpenMP)
tiny:
	$(CNAME) -O3 -ffunction-sections -fdata-sections $(MFLAGS) $(UFLAGS) $(LDFLAGS) $(ZFLAGS) $(BMPFLAGS) -o niimath

# tiny without nii2mesh (no OpenMP)
nano:
	$(CNAME) -O3 -ffunction-sections -fdata-sections $(UFLAGS) $(LDFLAGS) -o niimath

debug: CFLAGS   = -O0 -ffunction-sections -fdata-sections
debug:
	$(CNAME) $(CFLAGS) $(TFILES) $(DTFLAGS) $(GPL_SRCS) $(MFLAGS) $(AFLAGS) $(OMPFLAGS) $(ALSRC) $(BMPFLAGS) $(UFLAGS) $(ZFLAGS) $(ZSTDFLAGS) $(ZSTDLIB) $(OMPLINK) $(LDFLAGS) -o niimath

verbose:
	$(CNAME) -O0 -Wall -Wextra -Wno-sign-compare $(TFILES) $(DTFLAGS) $(GPL_SRCS) $(MFLAGS) $(AFLAGS) $(OMPFLAGS) $(ALSRC) $(UFLAGS) $(ZFLAGS) $(ZSTDFLAGS) $(ZSTDLIB) $(OMPLINK) $(LDFLAGS) -o niimath

# sanitize checks memory errors - similar to tiny/wasm
# on MacOS consider: export MallocNanoZone=0
# $(AFASTMATH) matches the shipped allineate FP build: without it NEWUOA converges
# poorly and runs near its iteration cap, which ASan instrumentation amplifies into a
# pseudo-hang (a 0.2s registration took 30+ min). -ffast-math is orthogonal to ASan's
# memory checking, so this does not weaken leak/overflow detection.
sanitize:
	$(CNAME) -O1 -g -Wno-deprecated -fsanitize=address -fno-omit-frame-pointer $(AFASTMATH) $(TFILES) $(DTFLAGS) $(GPL_SRCS) $(MFLAGS) $(AFLAGS) $(OMPFLAGS) $(ALSRC) $(UFLAGS) $(ZFLAGS) $(ZSTDFLAGS) $(ZSTDLIB) $(OMPLINK) $(LDFLAGS) -o niimath

# WASM build. The default (BSD) target is minimal. With GPL=1 it additionally links
# allineate (-deface/-allineate) and the spm_coreg module (-spm_coreg/-spm_deface,
# GPL-2) — both via $(WASM_AL)/$(GPL_SRCS), empty unless GPL=1. allineate.c/
# powell_newuoa.c are pre-built (wasm_al_objs) WITH -ffast-math to match the native
# build: that reassociation is load-bearing for the registration cost reduction —
# WITHOUT it the optimizer converges to a worse local minimum in WASM (verified:
# -deface masked 37% vs 60% native; fast-math restores 60%). Scoped to allineate/
# powell only so the GPL spm_coreg golden match and the BSD-default wasm are untouched.
# Every OpenMP construct in those files is #ifdef _OPENMP-guarded and _OPENMP is never
# set here, so the scalar fallbacks compile.
# STACK_SIZE is raised from emcc's 64 KB default because the GPL cost function keeps a
# 256x256 joint histogram (double H[65536] = 512 KB) on the stack (cost.c); the default
# overflows into a RuntimeError. Harmless for non-GPL builds (just reserves more).
# Pre-build allineate + powell with -ffast-math (matches the native allineate.o/
# powell_newuoa.o). Load-bearing: the cost-reduction reassociation it enables is what
# makes the WASM registration converge to the same optimum as native.
wasm_al_objs:
	emcc -O3 -ffast-math -fno-finite-math-only -DEMSCRIPTEN $(AFLAGS) -c allineate.c -o al_wasm.o
	emcc -O3 -ffast-math -fno-finite-math-only -DEMSCRIPTEN -c powell_newuoa.c -o pn_wasm.o

wasm: $(WASM_AL_PREREQ)
	emcc -O3 -DEMSCRIPTEN $(MFLAGS) $(DT_WASM) $(WASM_AL) $(GPL_SRCS) $(UFLAGS) $(ZFLAGS) $(BMPFLAGS) -s USE_ZLIB=1 -s STACK_SIZE=4194304 -s EXPORTED_RUNTIME_METHODS='["callMain", "ccall", "cwrap", "FS_createDataFile", "FS_readFile", "FS_unlink", "allocateUTF8", "getValue", "stringToUTF8", "setValue"]' -s ALLOW_MEMORY_GROWTH=1 -s WASM=1 -s EXPORT_ES6=1 -s MODULARIZE=1 -s EXPORTED_FUNCTIONS='["_main", "_malloc", "_free"]' -s INVOKE_RUN=0 -o ../js/src/niimath.js
	# hint: consider further optimizations:
	#   wasm-opt -O3 -o output.wasm niimath.wasm; rm niimath.wasm; mv output.wasm niimath.wasm

#WASM callMain details https://youtu.be/c8hZFtl8EuQ?si=NpAMjf7ka5XtsiYw
