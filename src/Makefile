# Regular use, set compiler name, compiler flags, openmp flags zlib flags
CNAME=gcc
#linker flags
LFLAGS=-lm -lz
#c flags
CFLAGS=-O3 -DHAVE_ZLIB -DFSLSTYLE -DPIGZ -DREJECT_COMPLEX
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

all:
	$(CNAME) $(CFLAGS) -o niimath niimath.c  bwlabel.c bw.c core.c tensor.c core32.c core64.c niftilib/nifti2_io.c znzlib/znzlib.c -I./niftilib -I./znzlib $(LFLAGS)
	
