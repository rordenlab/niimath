// nim2png_refactored.c
// Refactored nim2png implementation: single mapping for axes, smaller inner loop,
// normalized tile->value to signed index, and extracted helpers for maintainability.
//
// Based on your file. Reference version: :contentReference[oaicite:1]{index=1}

#include "filter.h"
#include <math.h>
#include "nifti_io.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
  #define strcasecmp  _stricmp
  #define strncasecmp _strnicmp
#endif

// The tool can either of two PNG ligraries:
// STB stands alone, you must provide stb_image_write.h/stb_image.h
// SPNG uses zlib with better compression
#ifdef HAVE_STB
	#define STB_IMAGE_WRITE_IMPLEMENTATION
	#include "stb_image_write.h"
#else

#include "spng.h"

int stbi_write_png(const char *filename, int x, int y, int comp, const void *data, int stride_bytes) {
	// Basic argument validation (same semantics as stb)
	if (!filename || !data || x <= 0 || y <= 0 || comp < 1 || comp > 4)
		return 0;
	// Require tightly-packed rows; minimal shim does not copy/pad rows.
	if ((size_t)stride_bytes != (size_t)x * (size_t)comp) {
		fprintf(stderr, "stbi_write_png shim (minimal): requires stride_bytes == x * comp\n");
		return 0;
	}
	// Create encoder context
	spng_ctx *enc = spng_ctx_new(SPNG_CTX_ENCODER);
	if (!enc) {
		fprintf(stderr, "stbi_write_png shim: spng_ctx_new failed\n");
		return 0;
	}
	spng_set_option(enc, SPNG_ENCODE_TO_BUFFER, 1);
	struct spng_ihdr ihdr;
	memset(&ihdr, 0, sizeof(ihdr));
	ihdr.width = (uint32_t)x;
	ihdr.height = (uint32_t)y;
	ihdr.bit_depth = 8;
	switch (comp) {
	case 1:
		ihdr.color_type = SPNG_COLOR_TYPE_GRAYSCALE;
		break;
	case 2:
		ihdr.color_type = SPNG_COLOR_TYPE_GRAYSCALE_ALPHA;
		break;
	case 3:
		ihdr.color_type = SPNG_COLOR_TYPE_TRUECOLOR;
		break;
	case 4:
	default:
		ihdr.color_type = SPNG_COLOR_TYPE_TRUECOLOR_ALPHA;
		break;
	}
	int ret = spng_set_ihdr(enc, &ihdr);
	if (ret != 0) {
		fprintf(stderr, "stbi_write_png shim: spng_set_ihdr failed (%d): %s\n", ret, spng_strerror(ret));
		spng_ctx_free(enc);
		return 0;
	}
	size_t img_size = (size_t)stride_bytes * (size_t)y;
	// Encode once according to IHDR
	ret = spng_encode_image(enc, data, img_size, SPNG_FMT_PNG, SPNG_ENCODE_FINALIZE);
	if (ret != 0) {
		fprintf(stderr, "stbi_write_png shim: spng_encode_image failed (%d): %s\n", ret, spng_strerror(ret));
		spng_ctx_free(enc);
		return 0;
	}
	size_t png_size = 0;
	void *png_buf = spng_get_png_buffer(enc, &png_size, &ret);
	if (!png_buf || ret != 0) {
		fprintf(stderr, "stbi_write_png shim: spng_get_png_buffer failed (%d): %s\n", ret, spng_strerror(ret));
		spng_ctx_free(enc);
		return 0;
	}
	FILE *f = fopen(filename, "wb");
	if (!f) {
		fprintf(stderr, "stbi_write_png shim: fopen failed for %s\n", filename);
		free(png_buf);
		spng_ctx_free(enc);
		return 0;
	}
	size_t wrote = fwrite(png_buf, 1, png_size, f);
	fclose(f);
	// free libspng buffer and context
	free(png_buf);
	spng_ctx_free(enc);
	return (wrote == png_size) ? 1 : 0;
}

#endif

// Helper to write scaled PNG
int write_png_scale(const char *filename, int x, int y, int comp,
					const void *data, int stride_bytes,
					float scale, int isLinear) {
	if (!filename || !data || x <= 0 || y <= 0 || comp <= 0 || scale <= 0.0f)
		return -1;
	/* If scale == 1, write directly */
	if (fabsf(scale - 1.0f) < 1e-6f) {
		int ok = stbi_write_png(filename, x, y, comp, data, stride_bytes);
		return ok ? 0 : -1;
	}

	/* Compute new size */
	int dst_w = (int)floor((double)x * (double)scale + 0.5);
	int dst_h = (int)floor((double)y * (double)scale + 0.5);
	if (dst_w <= 0 || dst_h <= 0) {
		fprintf(stderr, "write_png_scale: invalid scaled size %d×%d\n", dst_w, dst_h);
		return -1;
	}

	/* Choose filter */
	double (*filterf)(double) = isLinear ? triangle_filter : box_filter;
	double fwidth = isLinear ? 1.0 : 0.5;

	/* Wrap source buffer directly without copying */
	Image src;
	src.xsize = x;
	src.ysize = y;
	src.components = comp;
	src.span = stride_bytes > 0 ? stride_bytes : x * comp;
	src.data = (Pixel *)data; /* non-owning pointer */

	/* Allocate destination image */
	Image *dst = new_image(dst_w, dst_h, comp);
	if (!dst) {
		fprintf(stderr, "write_png_scale: out of memory allocating dst image\n");
		return -1;
	}

	/* Perform scaling */
	if (zoom(dst, &src, filterf, fwidth) != 0) {
		fprintf(stderr, "write_png_scale: zoom() failed\n");
		free_image(dst);
		return -1;
	}

	/* Write scaled PNG */
	int ret = stbi_write_png(filename, dst_w, dst_h, comp, dst->data, dst->span);
	free_image(dst);
	return ret;
}

// Enumeration of predefined color lookup tables
typedef enum {
	LUT_NONE = 0,
	LUT_GRAY,
	LUT_RED,
	LUT_GREEN,
	LUT_BLUE,
	LUT_CYAN,
	LUT_YELLOW,
	LUT_BLUEGREEN,
	LUT_REDYELLOW,
	LUT_VIRIDIS,
	LUT_INFERNO,
	LUT_MAGMA,
	LUT_PLASMA,
	LUT_MAX
} color_lut_t;

// Optional: LUT names for parsing from strings or CLI
static const char *LUT_NAMES[LUT_MAX] = {
	"none", "gray", "red", "green", "blue",
	"cyan", "yellow", "bluegreen", "redyellow", "viridis", "inferno", "magma", "plasma"};

// --------------------------- Types -------------------------------------
// Image-wide options (defaults applied unless per-tile flags override)
typedef struct {
	int isRadiological; // 1 = radiological (default), 0 = neurological
	int isLabel;		// 1 = draw labels by default
	int isNegativeColorMap;
	int interpolationOrder; // 0=nearest, 1= linear, [future 2=Mitchell]
	float scale;			// zoom factor
	float cal_min, cal_max, cal_min2, cal_max2;
	float overlayRGBA[4]; // overlay color + opacity (0..1) default {1,0,0,0.5}
	float RGBA[4];		  // base tint + opacity (alpha currently unused for base) default {1,1,1,1}
	float RGBA2[4];		  // background color + alpha default {0,0,0,0.5}
	color_lut_t LUT2;	  // new overlay color lookup table
	color_lut_t LUT;	  // new base color lookup table
} image_options_t;

// tile description produced by parsing the argv
typedef struct {
	char axis;		 // 'x','y','z' or 'N' for newline/separator
	long long value; // numeric slice index (signed integer index)
	size_t w;		 // tile pixel width (computed from nim dims)
	size_t h;		 // tile pixel height (computed from nim dims)
	size_t x_off;	 // x offset in row (filled during layout)
	size_t row;		 // row index (filled during layout)
	int isCrossLines;
} tile_t;

// layout summary returned from layout_tiles
typedef struct {
	size_t rows;		   // number of rows
	size_t total_width;	   // overall image width in pixels
	size_t total_height;   // overall image height in pixels
	size_t *row_widths;	   // length rows: width of each row
	size_t *row_heights;   // length rows: height of each row
	size_t *row_y_offsets; // length rows: vertical offset (y) of each row within final image
} layout_t;

// ------------------------ Small helpers --------------------------------

// check if token is parseable as a number (consumes entire string)
static int is_number_token(const char *s) {
	if (!s || *s == '\0')
		return 0;
	char *end = NULL;
	strtod(s, &end);
	return end && *end == '\0';
}

// clamp integer to [a,b]
static long long clamp_ll(long long v, long long a, long long b) {
	if (v < a)
		return a;
	if (v > b)
		return b;
	return v;
}

// convert a value in [cal_min, cal_max] to u8 [0..255]
static uint8_t float_to_u8(double val, double cal_min, double cal_max) {
	double denom = cal_max - cal_min;
	if (denom == 0.0)
		return 0;
	double t = (val - cal_min) / denom;
	if (t < 0.0)
		t = 0.0;
	if (t > 1.0)
		t = 1.0;
	return (uint8_t)(t * 255.0 + 0.5);
}

// burn a small 5x4 'L' or 'R' label into the RGB image buffer
// pixels: pointer to row-major RGB buffer (height*row_bytes bytes). row_bytes = width*3
// px,py coordinates of top-left of label
static void burn_label_rgb(uint8_t *pixels, size_t img_width, size_t img_height, size_t px, size_t py, char label) {
	// 5x4 bitmaps for L and R (last column removed)
	static const uint8_t L[5][4] = {
		{1, 0, 0, 0},
		{1, 0, 0, 0},
		{1, 0, 0, 0},
		{1, 0, 0, 0},
		{1, 1, 1, 1}};
	static const uint8_t R[5][4] = {
		{1, 1, 1, 0},
		{1, 0, 0, 1},
		{1, 1, 1, 0},
		{1, 0, 1, 0},
		{1, 0, 0, 1}};

	const uint8_t (*bmp)[4] = (label == 'R') ? R : L;
	const int w = 4; // width in pixels (columns)
	const int h = 5; // height in pixels (rows)
	const int channels = 3;

	if (px + (size_t)w > img_width)
		return;
	if (py + (size_t)h > img_height)
		return;

	for (int y = 0; y < h; ++y) {
		size_t yy = py + y;
		uint8_t *row_ptr = pixels + yy * (img_width * channels);
		for (int x = 0; x < w; ++x) {
			if (bmp[y][x]) {
				size_t xx = px + x;
				uint8_t *p = row_ptr + xx * channels;
				p[0] = 255;
				p[1] = 255;
				p[2] = 255;
			}
		}
	}
}

// write a pixel into RGB buffer (C replacement for previous lambda)
// pixels: RGB buffer, row_bytes = width*channels
static inline void write_pixel_rgb(uint8_t *pixels, size_t row_bytes, size_t W, size_t H,
								   size_t x, size_t y, uint8_t r, uint8_t g, uint8_t b) {
	const int channels = 3;
	if (x >= W || y >= H)
		return;
	size_t idx = y * row_bytes + x * channels;
	pixels[idx + 0] = r;
	pixels[idx + 1] = g;
	pixels[idx + 2] = b;
}

// ------------------------ Helpers for alias expansion --------------------

// helper to ensure capacity for tiles array
static int ensure_tile_capacity(tile_t **tiles, size_t *cap, size_t need) {
	if (*cap >= need)
		return 0;
	size_t newcap = *cap ? *cap * 2 : 8;
	while (newcap < need)
		newcap *= 2;
	tile_t *tmp = (tile_t *)realloc(*tiles, newcap * sizeof(tile_t));
	if (!tmp)
		return -1;
	*tiles = tmp;
	*cap = newcap;
	return 0;
}

// append a single tile for given axis and fractional value (in [0,1])
static int append_tile_for_axis_value(tile_t **tiles_ptr, size_t *n_ptr, size_t *cap_ptr,
									  const nifti_image *nim, char axis, double frac, int isCrossLines) {
	tile_t *tiles = *tiles_ptr;
	size_t n = *n_ptr;
	size_t cap = *cap_ptr;

	if (ensure_tile_capacity(&tiles, &cap, n + 1) != 0)
		return -1;

	size_t w = 0, h = 0;
	long long nx = nim->nx, ny = nim->ny, nz = nim->nz;
	long long value = 0;

	if (axis == 'x') {
		w = (size_t)ny;
		h = (size_t)nz;
		value = (long long)llround(frac * (nx - 1));
	} else if (axis == 'y') {
		w = (size_t)nx;
		h = (size_t)nz;
		value = (long long)llround(frac * (ny - 1));
	} else if (axis == 'z') {
		w = (size_t)nx;
		h = (size_t)ny;
		value = (long long)llround(frac * (nz - 1));
	} else
		return -1;

	// FSL-style: negative fraction means absolute slice number
	if (frac < 0.0)
		value = (long long)llround(fabs(frac));

	// Clamp to valid range
	if (axis == 'x')
		value = clamp_ll(value, 0, nx - 1);
	if (axis == 'y')
		value = clamp_ll(value, 0, ny - 1);
	if (axis == 'z')
		value = clamp_ll(value, 0, nz - 1);

	tiles[n].axis = axis;
	tiles[n].value = value;
	tiles[n].isCrossLines = isCrossLines;
	tiles[n].w = w;
	tiles[n].h = h;
	tiles[n].x_off = 0;
	tiles[n].row = 0;
	*tiles_ptr = tiles;
	*n_ptr = n + 1;
	*cap_ptr = cap;
	return 0;
}

// append an 'N' separator tile
static int append_separator(tile_t **tiles_ptr, size_t *n_ptr, size_t *cap_ptr) {
	tile_t *tiles = *tiles_ptr;
	size_t n = *n_ptr;
	size_t cap = *cap_ptr;

	if (ensure_tile_capacity(&tiles, &cap, n + 1) != 0)
		return -1;

	tiles[n].axis = 'N';
	tiles[n].value = 0;
	tiles[n].w = 0;
	tiles[n].h = 0;
	tiles[n].x_off = 0;
	tiles[n].row = 0;

	*tiles_ptr = tiles;
	*n_ptr = n + 1;
	*cap_ptr = cap;
	return 0;
}

/* helper: find LUT enum by name (case-insensitive). Returns -1 if not found. */
static int find_lut_by_name(const char *name) {
	if (!name)
		return -1;
	for (int k = 0; k < LUT_MAX; ++k) {
		if (strcasecmp(name, LUT_NAMES[k]) == 0)
			return k;
	}
	return -1;
}

// ------------------- Stage 1: parse_args_build_tiles --------------------
// parse_args_build_tiles: builds tiles AND updates image-wide options in *opts
static int parse_args_build_tiles(const nifti_image *nim,
								  int argc,
								  const char *argv[],
								  tile_t **out_tiles,
								  size_t *out_n_tiles,
								  image_options_t *opts) {
	if (!nim || argc <= 0 || !argv || !out_tiles || !out_n_tiles || !opts)
		return -1;

	tile_t *tiles = NULL;
	size_t n = 0, cap = 0;

	for (int i = 0; i < argc;) {
		const char *tok = argv[i];
		if (!tok) {
			++i;
			continue;
		}

		// ---------- Global color options (can appear anywhere) ----------
		// -N [0|1] or just -N (no arg -> set to 1)
		if (!strcmp(tok, "-N")) {
			if (i + 1 < argc && is_number_token(argv[i + 1])) {
				opts->isNegativeColorMap = atoi(argv[i + 1]) ? 1 : 0;
				i += 2;
			} else {
				// shorthand: -N alone -> enable negative color map
				opts->isNegativeColorMap = 1;
				++i;
			}
			continue;
		}
		/* Parse -c / -C: base or overlay LUT/RGBA */
		if ((!strcmp(tok, "-c")) || (!strcmp(tok, "-C"))) {
			int isOverlay = (tok[1] == 'C'); // uppercase = overlay (LUT2/RGBA2)
			color_lut_t *lut_field = isOverlay ? &opts->LUT2 : &opts->LUT;
			float *rgba_field = isOverlay ? opts->RGBA2 : opts->RGBA;

			if (i + 1 >= argc) {
				fprintf(stderr, "parse_args: %s requires arguments\n", tok);
				free(tiles);
				return -1;
			}
			/* numeric RGBA: expect 4 numbers */
			if (is_number_token(argv[i + 1])) {
				if (i + 4 >= argc) {
					fprintf(stderr, "parse_args: %s requires 4 numeric args\n", tok);
					free(tiles);
					return -1;
				}
				if (!is_number_token(argv[i + 1]) || !is_number_token(argv[i + 2]) ||
					!is_number_token(argv[i + 3]) || !is_number_token(argv[i + 4])) {
					fprintf(stderr, "parse_args: %s expects numeric args\n", tok);
					free(tiles);
					return -1;
				}

				rgba_field[0] = fminf(fmaxf((float)strtod(argv[i + 1], NULL), 0.0f), 1.0f);
				rgba_field[1] = fminf(fmaxf((float)strtod(argv[i + 2], NULL), 0.0f), 1.0f);
				rgba_field[2] = fminf(fmaxf((float)strtod(argv[i + 3], NULL), 0.0f), 1.0f);
				rgba_field[3] = fminf(fmaxf((float)strtod(argv[i + 4], NULL), 0.0f), 1.0f);
				*lut_field = LUT_NONE;

				i += 5;
				continue;
			}

			/* named LUT form: -c <name> [alpha] */
			int lut_id = find_lut_by_name(argv[i + 1]);
			if (lut_id < 0) {
				fprintf(stderr, "parse_args: unknown LUT name '%s'\n", argv[i + 1]);
				free(tiles);
				return -1;
			}
			*lut_field = (color_lut_t)lut_id;

			/* optional alpha after name */
			float alpha = 1.0f;
			if (i + 2 < argc && is_number_token(argv[i + 2])) {
				alpha = fminf(fmaxf((float)strtod(argv[i + 2], NULL), 0.0f), 1.0f);
				i += 3;
			} else {
				i += 2;
			}
			rgba_field[3] = alpha;
			continue;
		}
		if (!strcmp(tok, "-b")) { // RGBA2: -b R G B A
			if (i + 4 >= argc) {
				fprintf(stderr, "parse_args: -b requires 4 numeric args\n");
				free(tiles);
				return -1;
			}
			if (!is_number_token(argv[i + 1]) || !is_number_token(argv[i + 2]) || !is_number_token(argv[i + 3]) || !is_number_token(argv[i + 4])) {
				fprintf(stderr, "parse_args: -b expects numeric args\n");
				free(tiles);
				return -1;
			}
			opts->RGBA2[0] = fminf(fmaxf((float)strtod(argv[i + 1], NULL), 0.0f), 1.0f);
			opts->RGBA2[1] = fminf(fmaxf((float)strtod(argv[i + 2], NULL), 0.0f), 1.0f);
			opts->RGBA2[2] = fminf(fmaxf((float)strtod(argv[i + 3], NULL), 0.0f), 1.0f);
			opts->RGBA2[3] = fminf(fmaxf((float)strtod(argv[i + 4], NULL), 0.0f), 1.0f);
			i += 5;
			continue;
		}
		/* Calibration range: -t <min> <max> */
		if ((!strcmp(tok, "-t")) || (!strcmp(tok, "-T"))) {
			/* need two more args: argv[i+1], argv[i+2] */
			if (i + 2 >= argc) {
				fprintf(stderr, "parse_args: %s requires two numeric arguments (min max)\n", tok);
				free(tiles);
				return -1;
			}
			if (!is_number_token(argv[i + 1]) || !is_number_token(argv[i + 2])) {
				fprintf(stderr, "parse_args: %s expects two numeric args\n", tok);
				free(tiles);
				return -1;
			}

			float mn = (float)strtod(argv[i + 1], NULL);
			float mx = (float)strtod(argv[i + 2], NULL);

			if (mn >= mx) {
				fprintf(stderr, "parse_args: invalid range for %s: min >= max (%.3f >= %.3f)\n", tok, mn, mx);
				free(tiles);
				return -1;
			}

			if (!strcmp(tok, "-T")) {
				opts->cal_min2 = mn;
				opts->cal_max2 = mx;
			} else {
				opts->cal_min = mn;
				opts->cal_max = mx;
			}

			i += 3;
			continue;
		}

		// ---------- Global small flags: -f (radiological/neurological) and -u (labels) ----------
		// -f [0|1] or just -f (no arg -> set to 0)
		if (!strcmp(tok, "-f")) {
			if (i + 1 < argc && is_number_token(argv[i + 1])) {
				opts->isRadiological = atoi(argv[i + 1]) ? 1 : 0;
				i += 2;
			} else {
				// shorthand: -f alone -> set to 0 (neurological)
				opts->isRadiological = 0;
				++i;
			}
			continue;
		}
		if (!strcmp(tok, "-n")) {
			if (i + 1 < argc && is_number_token(argv[i + 1])) {
				opts->interpolationOrder = atoi(argv[i + 1]) ? 1 : 0;
				i += 2;
			} else {
				// shorthand: -f alone -> set to 0 (nearest)
				opts->interpolationOrder = 0;
				++i;
			}
			continue;
		}
		// -s <float> : global scale factor (must be positive)
		if (!strcmp(tok, "-s")) {
			if (i + 1 < argc && is_number_token(argv[i + 1])) {
				double v = strtod(argv[i + 1], NULL);
				if (!(isfinite(v) && v > 0.0)) {
					fprintf(stderr, "parse_args: -s requires a positive finite number\n");
					free(tiles);
					return -1;
				}
				opts->scale = (float)v;
				i += 2;
			} else {
				fprintf(stderr, "parse_args: -s requires a numeric argument\n");
				free(tiles);
				return -1;
			}
			continue;
		}
		// -u [0|1] or just -u (no arg -> set to 0)
		if (!strcmp(tok, "-u")) {
			if (i + 1 < argc && is_number_token(argv[i + 1])) {
				opts->isLabel = atoi(argv[i + 1]) ? 1 : 0;
				i += 2;
			} else {
				// shorthand: -u alone -> disable labels
				opts->isLabel = 0;
				++i;
			}
			continue;
		}
		int isCrossLines = 0;
		// ---------- Tile-related tokens ----------
		if (tok[0] == '-' && tok[1] != '\0') {
			char flag = tok[1];
			if (flag == 'X' || flag == 'Y' || flag == 'Z') {
				isCrossLines = 1;
				flag = tolower(flag);
			}
			if (flag == 'r') {
				// row separator
				if (append_separator(&tiles, &n, &cap) != 0) {
					free(tiles);
					return -1;
				}
				++i;
				continue;
			}

			// ALIASES: -a and -m expand into tile lists
			if (flag == 'a') {
				if (append_tile_for_axis_value(&tiles, &n, &cap, nim, 'x', 0.5, isCrossLines) != 0) {
					free(tiles);
					return -1;
				}
				if (append_tile_for_axis_value(&tiles, &n, &cap, nim, 'y', 0.5, isCrossLines) != 0) {
					free(tiles);
					return -1;
				}
				if (append_tile_for_axis_value(&tiles, &n, &cap, nim, 'z', 0.5, isCrossLines) != 0) {
					free(tiles);
					return -1;
				}
				++i;
				continue;
			}
			if (flag == 'm') {
				double vals[3] = {0.25, 0.5, 0.75};
				for (int k = 0; k < 3; ++k) {
					if (append_tile_for_axis_value(&tiles, &n, &cap, nim, 'x', vals[k], isCrossLines) != 0) {
						free(tiles);
						return -1;
					}
				}
				if (append_separator(&tiles, &n, &cap) != 0) {
					free(tiles);
					return -1;
				}
				for (int k = 0; k < 3; ++k) {
					if (append_tile_for_axis_value(&tiles, &n, &cap, nim, 'y', vals[k], isCrossLines) != 0) {
						free(tiles);
						return -1;
					}
				}
				if (append_separator(&tiles, &n, &cap) != 0) {
					free(tiles);
					return -1;
				}
				for (int k = 0; k < 3; ++k) {
					if (append_tile_for_axis_value(&tiles, &n, &cap, nim, 'z', vals[k], isCrossLines) != 0) {
						free(tiles);
						return -1;
					}
				}
				++i;
				continue;
			}
			if (flag == 'o') {
				// Determine voxel dimensions
				double nx = nim->nx;
				double ny = nim->ny;
				double nz = nim->nz;

				// Compute in-plane areas for each axis
				double area_x = fabs(ny * nz); // sagittal: y–z plane
				double area_y = fabs(nx * nz); // coronal: x–z plane
				double area_z = fabs(nx * ny); // axial: x–y plane
				flag = 'x';
				if ((area_y > area_x) && (area_y > area_z)) {
					flag = 'y';
				} else if ((area_z > area_x) && (area_z > area_y)) {
					flag = 'z';
				}
			}
			// axis flags: x,y,z collect numeric args
			if (flag != 'x' && flag != 'y' && flag != 'z') {
				++i;
				continue;
			}

			int j = i + 1;
			int found_any = 0;
			while (j < argc && is_number_token(argv[j])) {
				double val = strtod(argv[j], NULL);
				if (append_tile_for_axis_value(&tiles, &n, &cap, nim, flag, val, isCrossLines) != 0) {
					free(tiles);
					return -1;
				}
				found_any = 1;
				++j;
			}

			if (!found_any) {
				fprintf(stderr, "parse_args: axis -%c requires at least one numeric argument\n", flag);
				free(tiles);
				return -1;
			}
			i = j;
		} else {
			// non-flag token: skip (e.g., filename)
			++i;
		}
	}

	if (n == 0) {
		fprintf(stderr, "parse_args: no tiles requested\n");
		free(tiles);
		return -1;
	}

	*out_tiles = tiles;
	*out_n_tiles = n;
	return 0;
}

// ------------------- Stage 2: layout_tiles ------------------------------
// Arrange tiles into rows separated by 'N' tiles. Computes each tile's x_off and row.
// Produces layout_t with row metrics and total image size.
// Returns 0 on success, non-zero on failure.
static int layout_tiles(tile_t *tiles, size_t n_tiles, layout_t *out) {
	if (!tiles || n_tiles == 0 || !out)
		return -1;

	// temporary arrays for row widths/heights; start with reasonable capacity
	size_t rows_cap = 8;
	size_t rows = 0;
	size_t *row_widths = (size_t *)calloc(rows_cap, sizeof(size_t));
	size_t *row_heights = (size_t *)calloc(rows_cap, sizeof(size_t));
	if (!row_widths || !row_heights) {
		free(row_widths);
		free(row_heights);
		return -1;
	}

	size_t cur_row = 0;
	size_t cur_width = 0;
	size_t cur_height = 0;

	// First pass: assign tile->row and tile->x_off; compute row widths/heights
	for (size_t ti = 0; ti < n_tiles; ++ti) {
		tile_t *t = &tiles[ti];
		if (t->axis == 'N') {
			// finish current row
			if (rows + 1 > rows_cap) {
				size_t newcap = rows_cap * 2;
				size_t *rw = (size_t *)realloc(row_widths, newcap * sizeof(size_t));
				size_t *rh = (size_t *)realloc(row_heights, newcap * sizeof(size_t));
				if (!rw || !rh) {
					free(row_widths);
					free(row_heights);
					return -1;
				}
				row_widths = rw;
				row_heights = rh;
				rows_cap = newcap;
			}
			row_widths[rows] = cur_width;
			row_heights[rows] = cur_height;
			++rows;
			// start a new row
			cur_row++;
			cur_width = 0;
			cur_height = 0;
			// mark separator tile with special row
			t->row = (size_t)-1;
			t->x_off = 0;
		} else {
			t->row = cur_row;
			t->x_off = cur_width;
			cur_width += t->w;
			if (t->h > cur_height)
				cur_height = t->h;
		}
	}

	// finalize last row
	if (rows + 1 > rows_cap) {
		size_t newcap = rows_cap * 2;
		size_t *rw = (size_t *)realloc(row_widths, newcap * sizeof(size_t));
		size_t *rh = (size_t *)realloc(row_heights, newcap * sizeof(size_t));
		if (!rw || !rh) {
			free(row_widths);
			free(row_heights);
			return -1;
		}
		row_widths = rw;
		row_heights = rh;
		rows_cap = newcap;
	}
	row_widths[rows] = cur_width;
	row_heights[rows] = cur_height;
	++rows;

	// compute total width and total height and row_y_offsets
	size_t total_width = 0;
	for (size_t r = 0; r < rows; ++r) {
		if (row_widths[r] > total_width)
			total_width = row_widths[r];
	}
	size_t total_height = 0;
	size_t *row_y_offsets = (size_t *)malloc(rows * sizeof(size_t));
	if (!row_y_offsets) {
		free(row_widths);
		free(row_heights);
		return -1;
	}
	for (size_t r = 0; r < rows; ++r) {
		row_y_offsets[r] = total_height;
		total_height += row_heights[r];
	}

	// fill out result
	out->rows = rows;
	out->total_width = total_width;
	out->total_height = total_height;
	out->row_widths = row_widths;
	out->row_heights = row_heights;
	out->row_y_offsets = row_y_offsets;

	return 0;
}

// ------------------- Small mapping & pixel helper -----------------------

// Map local tile coordinates (xx,yy) into image voxel indices (xi,yi,zi)
// Caller ensures axis is 'x','y' or 'z'
static inline void map_tile_coords(char axis, long long slice_idx, size_t xx, size_t yy,
								   long long nx, long long ny, long long nz,
								   long long *xi, long long *yi, long long *zi) {
	if (axis == 'x') {
		*xi = slice_idx;
		*yi = (long long)xx;
		*zi = (long long)yy;
	} else if (axis == 'y') {
		*xi = (long long)xx;
		*yi = slice_idx;
		*zi = (long long)yy;
	} else { // 'z'
		*xi = (long long)xx;
		*yi = (long long)yy;
		*zi = slice_idx;
	}
}

// ------------------- Stage 3: render_and_write_png ----------------------
// ---------- LUT types + helpers (drop-in) ----------

typedef uint32_t lut_entry_t;	// packed 0xRRGGBBAA
typedef lut_entry_t lut_t[256]; // LUT table type

// pack/unpack macros: store bytes in order R(31..24) G(23..16) B(15..8) A(7..0)
#define LUT_PACK(r, g, b, a) ((uint32_t)(((uint32_t)(r) << 24) | ((uint32_t)(g) << 16) | ((uint32_t)(b) << 8) | ((uint32_t)(a))))
#define LUT_R(e) ((uint8_t)(((e) >> 24) & 0xFF))
#define LUT_G(e) ((uint8_t)(((e) >> 16) & 0xFF))
#define LUT_B(e) ((uint8_t)(((e) >> 8) & 0xFF))
#define LUT_A(e) ((uint8_t)(((e)) & 0xFF))

static inline uint8_t clamp_u8_int(int v) {
	if (v < 0)
		return 0;
	if (v > 255)
		return 255;
	return (uint8_t)v;
}

static void build_lut(lut_t lut,
					  const float RGBA0_in[4],
					  const float RGBA255_in[4],
					  uint8_t a,
					  color_lut_t whichLUT) {
	if (whichLUT == LUT_PLASMA || whichLUT == LUT_MAGMA || whichLUT == LUT_INFERNO || whichLUT == LUT_VIRIDIS) {
		/* --- PLASMA (fitted) --- */
		static const float plasma_R_coeffs[5] = {-4.82727553e-01f, 1.01534158e+00f, -9.08892060e-01f, 1.08043146e+00f, 6.12119843e-02f};
		static const float plasma_G_coeffs[5] = {-1.58245180e+00f, 3.26813002e+00f, -1.65150322e+00f, 1.86627404e-01f, 1.96796774e-02f};
		static const float plasma_B_coeffs[5] = {5.36338814e-01f, -1.02081633e+00f, 1.80382334e-02f, 4.30347618e-01f, 5.33572002e-01f};

		/* --- MAGMA (fitted) --- */
		static const float magma_R_coeffs[5] = {1.689417474714e+00f, -5.471495047641e+00f, 4.282408899801e+00f, 4.656519398970e-01f, -1.189895149846e-03f};
		static const float magma_G_coeffs[5] = {-2.336446198991e+00f, 5.430202008705e+00f, -2.896679494516e+00f, 8.157159254022e-01f, -9.312217129325e-03f};
		static const float magma_B_coeffs[5] = {7.454982027343e-01f, 3.885001753071e+00f, -7.488009034435e+00f, 3.709778066293e+00f, -5.099515131401e-02f};

		/* --- INFERNO (fitted) --- */
		static const float inferno_R_coeffs[5] = {9.283142715074e-01f, -2.319138388432e+00f, 1.913749123488e+00f, 3.403179719649e-01f, -8.344912439033e-03f};
		static const float inferno_G_coeffs[5] = {-6.861301782570e-01f, 1.485422391111e+00f, -7.792130192751e-01f, 2.619663903219e-01f, 3.894052053399e-03f};
		static const float inferno_B_coeffs[5] = {-1.619801618701e-01f, 1.076632861093e+00f, -2.660968833288e+00f, 2.007936604717e+00f, -2.597743648680e-02f};
		// viridis
		static const float viridis_R_coeffs[5] = {-4.17491528e+00f, 1.12224507e+01f, -7.48458801e+00f, 1.22910649e+00f, 2.35627699e-01f};
		static const float viridis_G_coeffs[5] = {-1.21256605e+00f, 2.25725461e+00f, -1.72980184e+00f, 1.58312791e+00f, 1.64777731e-03f};
		static const float viridis_B_coeffs[5] = {9.01900199e-01f, -1.62765043e+00f, -5.49789252e-01f, 9.66747911e-01f, 3.57845291e-01f};
		const float *cr, *cg, *cb;
		if (whichLUT == LUT_PLASMA) {
			cr = plasma_R_coeffs;
			cg = plasma_G_coeffs;
			cb = plasma_B_coeffs;
		} else if (whichLUT == LUT_MAGMA) {
			cr = magma_R_coeffs;
			cg = magma_G_coeffs;
			cb = magma_B_coeffs;
		} else if (whichLUT == LUT_INFERNO) {
			cr = inferno_R_coeffs;
			cg = inferno_G_coeffs;
			cb = inferno_B_coeffs;
		} else {
			cr = viridis_R_coeffs;
			cg = viridis_G_coeffs;
			cb = viridis_B_coeffs;
		}

		for (int u = 0; u < 256; ++u) {
			float t = (float)u / 255.0f;

			/* Horner's method for degree-4: (((c0*t + c1)*t + c2)*t + c3)*t + c4 */
			float R = (((cr[0] * t + cr[1]) * t + cr[2]) * t + cr[3]) * t + cr[4];
			float G = (((cg[0] * t + cg[1]) * t + cg[2]) * t + cg[3]) * t + cg[4];
			float B = (((cb[0] * t + cb[1]) * t + cb[2]) * t + cb[3]) * t + cb[4];

			/* clamp to [0,1] */
			if (R < 0.0f)
				R = 0.0f;
			else if (R > 1.0f)
				R = 1.0f;
			if (G < 0.0f)
				G = 0.0f;
			else if (G > 1.0f)
				G = 1.0f;
			if (B < 0.0f)
				B = 0.0f;
			else if (B > 1.0f)
				B = 1.0f;

			uint8_t rr = (uint8_t)(R * 255.0f + 0.5f);
			uint8_t gg = (uint8_t)(G * 255.0f + 0.5f);
			uint8_t bb = (uint8_t)(B * 255.0f + 0.5f);

			lut[u] = LUT_PACK(rr, gg, bb, a);
		}
		return;
	}

	// Define some RGBA constants (alpha ignored here)
	static const float BLACK[4] = {0.0f, 0.0f, 0.0f, 1.0f};
	static const float WHITE[4] = {1.0f, 1.0f, 1.0f, 1.0f};
	static const float RED[4] = {1.0f, 0.0f, 0.0f, 1.0f};
	static const float GREEN[4] = {0.0f, 1.0f, 0.0f, 1.0f};
	static const float BLUE[4] = {0.0f, 0.0f, 1.0f, 1.0f};
	static const float CYAN[4] = {0.0f, 1.0f, 1.0f, 1.0f};
	static const float YELLOW[4] = {1.0f, 1.0f, 0.0f, 1.0f};

	const float *RGBA0 = RGBA0_in;
	const float *RGBA255 = RGBA255_in;

	switch (whichLUT) {
	case LUT_GRAY:
		RGBA0 = BLACK;
		RGBA255 = WHITE;
		break;
	case LUT_RED:
		RGBA0 = BLACK;
		RGBA255 = RED;
		break;
	case LUT_GREEN:
		RGBA0 = BLACK;
		RGBA255 = GREEN;
		break;
	case LUT_BLUE:
		RGBA0 = BLACK;
		RGBA255 = BLUE;
		break;
	case LUT_CYAN:
		RGBA0 = BLACK;
		RGBA255 = CYAN;
		break;
	case LUT_YELLOW:
		RGBA0 = BLACK;
		RGBA255 = YELLOW;
		break;
	case LUT_BLUEGREEN:
		RGBA0 = BLUE;
		RGBA255 = GREEN;
		break;
	case LUT_REDYELLOW:
		RGBA0 = RED;
		RGBA255 = YELLOW;
		break;
	case LUT_NONE:
	default:
		/* LUT_NONE or unrecognized, use inputs as-is */
		break;
	}

	/* Interpolate linearly from RGBA0 → RGBA255, constant alpha 'a' */
	float r0 = fminf(fmaxf(RGBA0[0], 0.0f), 1.0f);
	float g0 = fminf(fmaxf(RGBA0[1], 0.0f), 1.0f);
	float b0 = fminf(fmaxf(RGBA0[2], 0.0f), 1.0f);

	float r1 = fminf(fmaxf(RGBA255[0], 0.0f), 1.0f);
	float g1 = fminf(fmaxf(RGBA255[1], 0.0f), 1.0f);
	float b1 = fminf(fmaxf(RGBA255[2], 0.0f), 1.0f);

	for (int u = 0; u < 256; ++u) {
		float t = (float)u / 255.0f;
		float R = r0 + (r1 - r0) * t;
		float G = g0 + (g1 - g0) * t;
		float B = b0 + (b1 - b0) * t;

		uint8_t rr = (uint8_t)(R * 255.0f + 0.5f);
		uint8_t gg = (uint8_t)(G * 255.0f + 0.5f);
		uint8_t bb = (uint8_t)(B * 255.0f + 0.5f);

		lut[u] = LUT_PACK(rr, gg, bb, a);
	}
}

// Fast getters using the LUTs
static inline void get_base_rgb_from_lut(float v, float cal_min, float cal_max,
										 const lut_t base_lut,
										 uint8_t *r, uint8_t *g, uint8_t *b) {
	uint8_t u = float_to_u8((double)v, (double)cal_min, (double)cal_max);
	uint32_t e = base_lut[u];
	*r = LUT_R(e);
	*g = LUT_G(e);
	*b = LUT_B(e);
}

// Fast blend: base from base_lut, overlay from overlay_lut. overlay_lut alpha channel used.
// If overlay intensity maps to 0 => return base directly.
static inline void get_blended_rgb_from_luts(float v_base, float cal_min, float cal_max,
											 float v_overlay, float cal_min2, float cal_max2,
											 const lut_t base_lut,
											 const lut_t negative_overlay_lut,
											 const lut_t overlay_lut,
											 uint8_t *r, uint8_t *g, uint8_t *b,
											 const int isNegativeColorMap) {
	uint8_t ub = float_to_u8((double)v_base, (double)cal_min, (double)cal_max);
	uint32_t be = base_lut[ub];
	uint8_t br = LUT_R(be), bg = LUT_G(be), bb = LUT_B(be);

	uint8_t uov = float_to_u8((double)v_overlay, (double)cal_min2, (double)cal_max2);
	float v_overlay2 = v_overlay;
	if ((v_overlay < 0.0) && (isNegativeColorMap) && (cal_min2 > 0)) {
		uov = float_to_u8((double)-v_overlay, (double)cal_min2, (double)cal_max2);
		v_overlay2 = -v_overlay;
	}
	if (v_overlay2 < cal_min2) {
		// if (uov == 0) {
		*r = br;
		*g = bg;
		*b = bb;
		return;
	}

	uint32_t oe = overlay_lut[uov];
	if ((v_overlay < 0.0) && (isNegativeColorMap) && (cal_min2 > 0)) {
		oe = negative_overlay_lut[uov];
	}
	uint8_t or_ = LUT_R(oe), og_ = LUT_G(oe), ob_ = LUT_B(oe);
	uint8_t alpha = LUT_A(oe); // 0..255

	if (alpha == 0) {
		*r = br;
		*g = bg;
		*b = bb;
		return;
	}

	// integer composite: out = ((255-alpha)*base + alpha*overlay)/255
	int out_r = ((255 - (int)alpha) * (int)br + (int)alpha * (int)or_ + 127) / 255;
	int out_g = ((255 - (int)alpha) * (int)bg + (int)alpha * (int)og_ + 127) / 255;
	int out_b = ((255 - (int)alpha) * (int)bb + (int)alpha * (int)ob_ + 127) / 255;

	*r = (uint8_t)clamp_u8_int(out_r);
	*g = (uint8_t)clamp_u8_int(out_g);
	*b = (uint8_t)clamp_u8_int(out_b);
}

// ---------- Replacement render_and_write_png() that builds LUTs and uses them ----------
static int render_and_write_png(const nifti_image *nim,
								nifti_image *nim2,
								const tile_t *tiles,
								size_t n_tiles,
								const layout_t *layout,
								const char *out_png_path,
								const image_options_t *opts) {
	if (!nim || !tiles || n_tiles == 0 || !layout || !out_png_path || !opts)
		return -1;

	const int channels = 3;
	size_t W = layout->total_width;
	size_t H = layout->total_height;
	if (W == 0 || H == 0) {
		fprintf(stderr, "render: zero output dimension\n");
		return -1;
	}

	// check calibration
	if (!(isfinite(opts->cal_min) && isfinite(opts->cal_max) && opts->cal_min < opts->cal_max)) {
		fprintf(stderr, "render: nim cal_min and cal_max must be finite and cal_min < cal_max (%g %g)\n",
				opts->cal_min, opts->cal_max);
		return -1;
	}
	float cal_min = opts->cal_min;
	float cal_max = opts->cal_max;

	float cal_min2 = opts->cal_min2;
	float cal_max2 = opts->cal_max2;

	int have_nim2 = 0;
	float *data2 = NULL;
	if (nim2) {
		if (!(isfinite(opts->cal_min2) && isfinite(opts->cal_min2) && opts->cal_min2 < opts->cal_max2)) {
			fprintf(stderr, "render: overlay (nim2) cal_min and cal_max must be finite and cal_min < cal_max (%g %g)\n",
					opts->cal_min2, opts->cal_max2);
			return -1;
		}
		if ((nim2->nx != nim->nx) || (nim2->ny != nim->ny) || (nim2->nz != nim->nz)) {
			fprintf(stderr, "render: overlay dimensions do not match\n");
			return -1;
		}
		if (!nim2->data) {
			fprintf(stderr, "render: overlay (nim2) data pointer is NULL\n");
			return -1;
		}
		data2 = (float *)nim2->data;
		have_nim2 = 1;
	}

	// allocate RGB buffer and fill with background color
	size_t row_bytes = W * channels;
	size_t buf_size = row_bytes * H;
	uint8_t *rgb = (uint8_t *)calloc(1, buf_size);
	if (!rgb) {
		fprintf(stderr, "render: out of memory\n");
		return -1;
	}

	uint8_t bg_r = (uint8_t)(opts->RGBA2[0] * 255.0f + 0.5f);
	uint8_t bg_g = (uint8_t)(opts->RGBA2[1] * 255.0f + 0.5f);
	uint8_t bg_b = (uint8_t)(opts->RGBA2[2] * 255.0f + 0.5f);
	for (size_t y = 0; y < H; ++y) {
		uint8_t *row = rgb + y * row_bytes;
		for (size_t x = 0; x < W; ++x) {
			size_t idx = x * 3;
			row[idx + 0] = bg_r;
			row[idx + 1] = bg_g;
			row[idx + 2] = bg_b;
		}
	}

	// data pointer for nim
	float *data = (float *)nim->data;
	if (!data) {
		free(rgb);
		fprintf(stderr, "render: nim->data is NULL\n");
		return -1;
	}

	// Build LUTs from opts->RGBA and opts->overlayRGBA (simple builder)
	lut_t base_lut;
	lut_t overlay_lut;
	lut_t negative_overlay_lut;
	const float RGBA0[4] = {0.0f, 0.0f, 0.0f, 0.0f};
	build_lut(base_lut, RGBA0, opts->RGBA, 255, opts->LUT);
	build_lut(overlay_lut, RGBA0, opts->overlayRGBA, (uint8_t)(opts->overlayRGBA[3] * 255.0), opts->LUT2);
	build_lut(negative_overlay_lut, RGBA0, opts->overlayRGBA, (uint8_t)(opts->overlayRGBA[3] * 255.0), LUT_BLUEGREEN);

	// iterate tiles and rasterize each tile's slice into the appropriate rectangle
	long long nx = nim->nx, ny = nim->ny, nz = nim->nz;
	size_t slice_size = (size_t)nx * (size_t)ny;

	for (size_t ti = 0; ti < n_tiles; ++ti) {
		const tile_t *t = &tiles[ti];
		if (t->axis == 'N')
			continue;

		size_t row = t->row;
		if (row == (size_t)-1 || row >= layout->rows)
			continue;
		size_t base_x = t->x_off;
		size_t base_y = layout->row_y_offsets[row];

		long long slice_idx = clamp_ll(t->value, 0, (t->axis == 'x' ? nx - 1 : (t->axis == 'y' ? ny - 1 : nz - 1)));

		for (size_t yy = 0; yy < t->h; ++yy) {
			for (size_t xx = 0; xx < t->w; ++xx) {
				long long xi, yi, zi;
				map_tile_coords(t->axis, slice_idx, xx, yy, nx, ny, nz, &xi, &yi, &zi);

				if (xi < 0 || xi >= nx || yi < 0 || yi >= ny || zi < 0 || zi >= nz)
					continue;

				size_t data_idx = (size_t)xi + (size_t)yi * (size_t)nx + (size_t)zi * slice_size;
				size_t dst_x;
				if (t->axis == 'x') {
					// Flip sagittal slices: lowest x-indices on the left
					dst_x = base_x + xx;
				} else {
					// Keep current radiological vs neurological behavior for other planes
					dst_x = opts->isRadiological ? (base_x + (t->w - 1 - xx)) : (base_x + xx);
				}
				// size_t dst_x = opts->isRadiological ? (base_x + (t->w - 1 - xx)) : (base_x + xx);
				size_t dst_y = base_y + (t->h - 1 - yy);

				uint8_t rr, gg, bb;
				if (have_nim2) {
					float v = data[data_idx];
					float v2 = data2[data_idx];
					get_blended_rgb_from_luts(v, cal_min, cal_max, v2, cal_min2, cal_max2,
											  base_lut, negative_overlay_lut, overlay_lut, &rr, &gg, &bb, opts->isNegativeColorMap);
				} else {
					float v = data[data_idx];
					get_base_rgb_from_lut(v, cal_min, cal_max, base_lut, &rr, &gg, &bb);
				}
				write_pixel_rgb(rgb, row_bytes, W, H, dst_x, dst_y, rr, gg, bb);
			}
		}

		// draw cross-lines for this tile if requested
		if (t->isCrossLines) {
			for (size_t ot_i = 0; ot_i < n_tiles; ++ot_i) {
				const tile_t *ot = &tiles[ot_i];
				if (ot == t)
					continue;
				if (ot->axis == 'N')
					continue;
				if (ot->axis == t->axis)
					continue; // only other axes
				if (ot->isCrossLines)
					continue; // only if other tile marked

				// clamp other tile's slice into valid range for its axis
				long long ov = clamp_ll(ot->value,
										0,
										(ot->axis == 'x' ? nx - 1 : (ot->axis == 'y' ? ny - 1 : nz - 1)));

				// determine whether this other tile fixes an xx (vertical) or yy (horizontal) coord
				int draw_vertical = 0, draw_horizontal = 0;
				size_t local_x = 0, local_y = 0;

				if (t->axis == 'x') {
					// local xx -> yi, local yy -> zi
					if (ot->axis == 'y') {
						draw_vertical = 1;
						local_x = (size_t)ov;
					} else if (ot->axis == 'z') {
						draw_horizontal = 1;
						local_y = (size_t)ov;
					}
				} else if (t->axis == 'y') {
					// local xx -> xi, local yy -> zi
					if (ot->axis == 'x') {
						draw_vertical = 1;
						local_x = (size_t)ov;
					} else if (ot->axis == 'z') {
						draw_horizontal = 1;
						local_y = (size_t)ov;
					}
				} else { // t->axis == 'z'
					// local xx -> xi, local yy -> yi
					if (ot->axis == 'x') {
						draw_vertical = 1;
						local_x = (size_t)ov;
					} else if (ot->axis == 'y') {
						draw_horizontal = 1;
						local_y = (size_t)ov;
					}
				}

				// draw vertical line at local_x
				if (draw_vertical && local_x < t->w) {
					// compute dst_x using same radiological rule as pixel rasterization
					size_t dst_x;
					if (t->axis == 'x') {
						dst_x = base_x + local_x;
					} else {
						dst_x = opts->isRadiological ? (base_x + (t->w - 1 - local_x)) : (base_x + local_x);
					}
					// draw down the tile height
					for (size_t yy2 = 0; yy2 < t->h; ++yy2) {
						size_t dst_y = base_y + (t->h - 1 - yy2);
						// bounds check (defensive)
						if (dst_x >= W || dst_y >= H)
							continue;
						write_pixel_rgb(rgb, row_bytes, W, H, dst_x, dst_y, 255, 255, 255);
					}
				}

				// draw horizontal line at local_y
				if (draw_horizontal && local_y < t->h) {
					size_t dst_y = base_y + (t->h - 1 - local_y);
					for (size_t xx2 = 0; xx2 < t->w; ++xx2) {
						size_t dst_x;
						if (t->axis == 'x') {
							dst_x = base_x + xx2;
						} else {
							dst_x = opts->isRadiological ? (base_x + (t->w - 1 - xx2)) : (base_x + xx2);
						}
						if (dst_x >= W || dst_y >= H)
							continue;
						write_pixel_rgb(rgb, row_bytes, W, H, dst_x, dst_y, 255, 255, 255);
					}
				}
			} // end other-tile loop
		}

		// optionally draw labels
		if ((opts->isLabel) && (t->axis != 'x')) {
			size_t inset = 2;
			size_t lab_y = base_y + (t->h - 2 - 5 + 1); // base_y + t->h - 6
			if (lab_y + 5 <= H) {
				char left_label = opts->isRadiological ? 'R' : 'L';
				burn_label_rgb(rgb, W, H, base_x + inset, lab_y, left_label);
				char right_label = opts->isRadiological ? 'L' : 'R';
				if (t->w >= 7) {
					size_t px = base_x + t->w - inset - 5;
					burn_label_rgb(rgb, W, H, px, lab_y, right_label);
				}
			}
		}
	}

	// write PNG (with optional scaling)
	int ret = write_png_scale(out_png_path, (int)W, (int)H, channels, rgb, (int)row_bytes, opts->scale, opts->interpolationOrder);
	free(rgb);
	return ret;
}

// ------------------- Public entry: nim2png ---------------------------
// modernized signature: const char *argv[]
// This function ties the three stages together and manages memory.
int nim2png(nifti_image *nim, nifti_image *nim2, int argc, const char *argv[], const char *out_png_path) {
	if (!nim || argc < 0 || !argv || !out_png_path)
		return -1;

	// Stage 1: parse arguments -> tiles
	tile_t *tiles = NULL;
	size_t n_tiles = 0;

	// defaults
	image_options_t opts;
	opts.isRadiological = 1;
	opts.isLabel = 1;
	opts.isNegativeColorMap = 0;
	opts.interpolationOrder = 1;
	opts.scale = 1.0;
	opts.LUT2 = LUT_NONE;
	opts.LUT = LUT_NONE;
	opts.overlayRGBA[0] = 1.0f;
	opts.overlayRGBA[1] = 0.0f;
	opts.overlayRGBA[2] = 0.0f;
	opts.overlayRGBA[3] = 0.5f;
	opts.RGBA[0] = 1.0f;
	opts.RGBA[1] = 1.0f;
	opts.RGBA[2] = 1.0f;
	opts.RGBA[3] = 1.0f;
	opts.RGBA2[0] = 0.0f;
	opts.RGBA2[1] = 0.0f;
	opts.RGBA2[2] = 0.0f;
	opts.RGBA2[3] = 1.0f;
	opts.cal_min = nim->cal_min;
	opts.cal_max = nim->cal_max;
	if (nim2) {
		opts.cal_min2 = nim2->cal_min;
		opts.cal_max2 = nim2->cal_max;
	}
	if (parse_args_build_tiles(nim, argc, argv, &tiles, &n_tiles, &opts) != 0) {
		fprintf(stderr, "-bitmap failed to parse args (e.g. niimath nifti.nii -bitmap -a all.png\n");
		return -1;
	}

	// Stage 2: layout tiles -> layout_t
	layout_t layout;
	memset(&layout, 0, sizeof(layout));
	if (layout_tiles(tiles, n_tiles, &layout) != 0) {
		fprintf(stderr, "-bitmap failed to layout tiles\n");
		free(tiles);
		return -1;
	}

	// Stage 3: render and write PNG
	int rc = render_and_write_png(nim, nim2, tiles, n_tiles, &layout, out_png_path, &opts);

	// cleanup
	free(tiles);
	free(layout.row_widths);
	free(layout.row_heights);
	free(layout.row_y_offsets);
	return rc;
}
