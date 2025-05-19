# Regular use, set compiler name, compiler flags, openmp flags zlib flags
CNAME=gcc
#universal flags always included
UFLAGS= niimath.c core.c core32.c niftilib/nifti2_io.c znzlib/znzlib.c -I./niftilib -I./znzlib -DFSLSTYLE -DPIGZ -DREJECT_COMPLEX -lm
#experimental conform feature
UFLAGS+= conform.c -DHAVE_CONFORM
UFLAGS64= -DHAVE_64BITS core64.c


#Zlib flags
ZFLAGS= -DHAVE_ZLIB -lz
#compiler flags
CFLAGS= -O3
# CFLAGS= -O0 -Wall
# GIfTI, obj, vtk, stl, etc support
GFILES= -DHAVE_FORMATS base64.c
#Butterworth bandpass flags
BFLAGS= -DHAVE_BUTTERWORTH bw.c
#Tensor decomposition
TFILES= -DHAVE_TENSOR tensor.c
#Marching Cubes flag
# MFLAGS= -DNII2MESH meshify.c quadric.c bwlabel.c MarchingCubes.c radixsort.c fdr.c
MFLAGS= -DNII2MESH meshify.c quadric.c bwlabel.c radixsort.c fdr.c -DUSE_CLASSIC_CUBES oldcubes.c

ifeq ($(CXX),g++)
	CFLAGS += -std=gnu99
endif

#run "OMP=1 make -j" for both OpenMP and CloudFlare
ifeq "$(OMP)" "1"
	CNAME= gcc-9
	CFLAGS+=-fopenmp -flto -static-libgcc -static-libstdc++ -L.
endif

#run "MESH=0 make" for minimal compile without nii2mesh functions
ifeq "$(MESH)" "0"
	MFLAGS=
endif

all:
	$(CNAME) $(CFLAGS) $(BFLAGS) $(TFILES) $(GFILES) $(MFLAGS) $(UFLAGS) $(UFLAGS64) $(ZFLAGS) -flto -o niimath

#issue30: static build 
static:
	gcc -O3 -static -std=gnu99 $(BFLAGS) $(TFILES) $(GFILES) $(MFLAGS) $(UFLAGS) $(UFLAGS64) $(ZFLAGS) -flto -o niimath
	
# tiny: terminal executable to emulate WASM
tiny:
	$(CNAME) $(CFLAGS) $(MFLAGS) $(UFLAGS) -o niimath

# tiny without nii2mesh
nano:
	$(CNAME) $(CFLAGS) $(UFLAGS) -o niimath

# quick and dirty
debug:
	$(CNAME) -O0 $(MFLAGS) $(UFLAGS) $(ZFLAGS) -o niimath

verbose:
	$(CNAME) -O0 -Wall -Wextra -Wno-sign-compare $(MFLAGS) $(UFLAGS) $(ZFLAGS) -o niimath

# sanitize checks memory errors - similar to tiny/wasm
# on MacOS consider: export MallocNanoZone=0
sanitize:
	$(CNAME) -O1 -g -Wno-deprecated -fsanitize=address -fno-omit-frame-pointer $(MFLAGS) $(UFLAGS) $(ZFLAGS) -o niimath

wasm:
	emcc -O3 $(MFLAGS) $(UFLAGS) $(ZFLAGS) -s USE_ZLIB=1 -s DEMANGLE_SUPPORT=1 -s EXPORTED_RUNTIME_METHODS='["callMain", "ccall", "cwrap", "FS_createDataFile", "FS_readFile", "FS_unlink", "allocateUTF8", "getValue", "stringToUTF8", "setValue"]' -s ALLOW_MEMORY_GROWTH=1 -s WASM=1 -s EXPORT_ES6=1 -s MODULARIZE=1 -s EXPORTED_FUNCTIONS='["_main", "_malloc", "_free"]' -s INVOKE_RUN=0 -o ../js/src/niimath.js
	# hint: consider further optimizations:
	#   wasm-opt -O3 -o output.wasm niimath.wasm; rm niimath.wasm; mv output.wasm niimath.wasm
	
#WASM callMain details https://youtu.be/c8hZFtl8EuQ?si=NpAMjf7ka5XtsiYw
