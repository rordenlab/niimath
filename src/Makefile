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
ifneq "$(OMP)" "0"
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

ifeq ($(UNAME_S),Darwin)
	# macOS linker does not support --gc-sections
	LDFLAGS =
endif

all: $(ALOBJS)
	$(CNAME) $(CFLAGS) $(BFLAGS) $(TFILES) $(GFILES) $(AFLAGS) $(ALOBJS) $(MFLAGS) $(UFLAGS) $(UFLAGS64) $(ZFLAGS) $(ZSTDFLAGS) $(ZSTDLIB) $(OMPLINK) $(LDFLAGS) $(BMPFLAGS) -flto -o niimath

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
	gcc -O3 -static -std=gnu99 $(BFLAGS) $(TFILES) $(GFILES) $(MFLAGS) $(AFLAGS) $(AFASTMATH) $(OMPFLAGS) $(ALSRC) $(UFLAGS) $(UFLAGS64) $(ZFLAGS) $(ZSTDFLAGS) $(ZSTDLIB) $(OMPLINK) $(LDFLAGS) -flto -o niimath

# tiny: terminal executable to emulate WASM (no OpenMP)
tiny:
	$(CNAME) -O3 -ffunction-sections -fdata-sections $(MFLAGS) $(UFLAGS) $(LDFLAGS) $(ZFLAGS) $(BMPFLAGS) -o niimath

# tiny without nii2mesh (no OpenMP)
nano:
	$(CNAME) -O3 -ffunction-sections -fdata-sections $(UFLAGS) $(LDFLAGS) -o niimath

debug: CFLAGS   = -O0 -ffunction-sections -fdata-sections
debug:
	$(CNAME) $(CFLAGS) $(MFLAGS) $(AFLAGS) $(OMPFLAGS) $(ALSRC) $(BMPFLAGS) $(UFLAGS) $(ZFLAGS) $(ZSTDFLAGS) $(ZSTDLIB) $(OMPLINK) $(LDFLAGS) -o niimath

verbose:
	$(CNAME) -O0 -Wall -Wextra -Wno-sign-compare $(MFLAGS) $(AFLAGS) $(OMPFLAGS) $(ALSRC) $(UFLAGS) $(ZFLAGS) $(ZSTDFLAGS) $(ZSTDLIB) $(OMPLINK) $(LDFLAGS) -o niimath

# sanitize checks memory errors - similar to tiny/wasm
# on MacOS consider: export MallocNanoZone=0
sanitize:
	$(CNAME) -O1 -g -Wno-deprecated -fsanitize=address -fno-omit-frame-pointer $(MFLAGS) $(AFLAGS) $(OMPFLAGS) $(ALSRC) $(UFLAGS) $(ZFLAGS) $(ZSTDFLAGS) $(ZSTDLIB) $(OMPLINK) $(LDFLAGS) -o niimath

wasm:
	emcc -O3 $(MFLAGS) $(UFLAGS) $(ZFLAGS) $(BMPFLAGS) -s USE_ZLIB=1 -s EXPORTED_RUNTIME_METHODS='["callMain", "ccall", "cwrap", "FS_createDataFile", "FS_readFile", "FS_unlink", "allocateUTF8", "getValue", "stringToUTF8", "setValue"]' -s ALLOW_MEMORY_GROWTH=1 -s WASM=1 -s EXPORT_ES6=1 -s MODULARIZE=1 -s EXPORTED_FUNCTIONS='["_main", "_malloc", "_free"]' -s INVOKE_RUN=0 -o ../js/src/niimath.js
	# hint: consider further optimizations:
	#   wasm-opt -O3 -o output.wasm niimath.wasm; rm niimath.wasm; mv output.wasm niimath.wasm

#WASM callMain details https://youtu.be/c8hZFtl8EuQ?si=NpAMjf7ka5XtsiYw
