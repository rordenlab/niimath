#ifndef _PRINT_H_
#define _PRINT_H_

#ifdef USING_R
	#define R_NO_REMAP
	#define R_USE_C99_IN_CXX
	#include <R_ext/Print.h>
	#include <R_ext/Error.h>
	#define printfx(...) REprintf(__VA_ARGS__)
	int niimath_rand ();
#else
	#ifdef __EMSCRIPTEN__
		#define printfx(...) printf( __VA_ARGS__)
	#else
		#include <stdio.h>
		#define printfx(...) fprintf(stderr, __VA_ARGS__)
	#endif
#endif // USING_R

#endif
