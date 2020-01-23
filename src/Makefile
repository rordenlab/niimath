# Regular use, set compiler name, compiler flags, openmp flags zlib flags
CNAME=gcc
#linker flags
LFLAGS=-lm -lz
#c flags
CFLAGS=-O3 -DHAVE_ZLIB

#run "OMP=1 make -j" for both OpenMP and CloudFlare
ifeq "$(OMP)" "1"
	CNAME= gcc-9
	CFLAGS+=-fopenmp -flto
endif

all:
	$(CNAME) $(CFLAGS) $(CFLAGS) -o niimath niimath.c core.c tensor.c core32.c core64.c niftilib/nifti2_io.c znzlib/znzlib.c -I./niftilib -I./znzlib $(LFLAGS)
	
