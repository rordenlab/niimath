# Regular use, set compiler name, compiler flags, openmp flags zlib flags
CNAME=gcc
#linker flags
LFLAGS=-lm -lz
#c flags
CFLAGS=-O3 -DHAVE_ZLIB -DFSLSTYLE -DPIGZ -DREJECT_COMPLEX
#CFLAGS=-O3 -ffast-math -DHAVE_ZLIB -DFSLSTYLE -DPIGZ -DREJECT_COMPLEX -DNII2MESH
# Westmere required for fast gz https://github.com/tikv/tikv/issues/4843
#CFLAGS=-O3 -DHAVE_ZLIB -march=westmere

ifeq ($(CXX),g++)
	CFLAGS += -std=gnu99
endif

#run "OMP=1 make -j" for both OpenMP and CloudFlare
ifeq "$(OMP)" "1"
	CNAME= gcc-9
	CFLAGS+=-fopenmp -flto
	#to get full file path:
	#gcc-9 --print-file-name=libgomp.a
	#LFLAGS+= /usr/local/Cellar/gcc/9.2.0_3/lib/gcc/9/gcc/x86_64-apple-darwin18/9.2.0/../../../libgomp.a
	LFLAGS+= -static-libgcc -static-libstdc++ -L.
endif

MFLAGS= -DNII2MESH MarchingCubes.c meshify.c quadric.c
#run "MESH=0 make" for minimal compile without nii2mesh functions
ifeq "$(MESH)" "0"
	MFLAGS=
endif

all:
	$(CNAME) $(CFLAGS) $(MFLAGS) niimath.c base64.c radixsort.c fdr.c bwlabel.c bw.c core.c tensor.c core32.c core64.c niftilib/nifti2_io.c znzlib/znzlib.c -I./niftilib -I./znzlib $(LFLAGS) -o niimath
	