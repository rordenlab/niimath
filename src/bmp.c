// nim2png_refactored.c
// Refactored nim2png implementation: single mapping for axes, smaller inner loop,
// normalized tile->value to signed index, and extracted helpers for maintainability.
//
// Based on your file. Reference version: :contentReference[oaicite:1]{index=1}

#include <math.h>
#include <nifti2_io.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// stb_image_write for PNG output (single-file header)
// If you already include STB_IMAGE_WRITE_IMPLEMENTATION once elsewhere, remove the define here.
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

// Dale Schumacher, "General Filtered Image Rescaling", in Graphics Gems III, Academic Press, 1992
// https://github.com/erich666/GraphicsGems/tree/master/gemsiii
// write_png_scale.c
// Requires filter.c (Graphics Gems zoom/filter functions) to be compiled and linked.
// See the uploaded filter.c (Graphics Gems) for zoom() and filters. :contentReference[oaicite:1]{index=1}
/* --- Begin minimal Graphics Gems resampling code (transplanted) --- */
/* Source: Dale Schumacher, "General Filtered Image Rescaling", Graphics Gems III.
   Public domain. See filter.c for full original code. :contentReference[oaicite:1]{index=1} */

#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Pixel and Image types used by zoom() */
typedef unsigned char Pixel;

typedef struct {
    int xsize;
    int ysize;
    Pixel *data;
    int span;
} Image;

/* small macros */
#ifndef CLAMP
#define CLAMP(v, lo, hi) ( ((v) < (lo)) ? (lo) : (((v) > (hi)) ? (hi) : (v)) )
#endif

/* create/destroy Image */
Image * new_image(int xsize, int ysize) {
    Image *image;
    if ((image = (Image *)malloc(sizeof(Image))) &&
        (image->data = (Pixel *)calloc((size_t)ysize, (size_t)xsize))) {
        image->xsize = xsize;
        image->ysize = ysize;
        image->span  = xsize;
    } else {
        free(image);
        image = NULL;
    }
    return image;
}

void free_image(Image *image) {
    if (!image) return;
    free(image->data);
    free(image);
}

/* basic pixel access helpers used by zoom() */
Pixel get_pixel(Image *image, int x, int y) {
    if (!image) return 0;
    if ((x < 0) || (x >= image->xsize) || (y < 0) || (y >= image->ysize)) return 0;
    return image->data[y * image->span + x];
}

void get_row(Pixel *row, Image *image, int y) {
    if (!image) return;
    if ((y < 0) || (y >= image->ysize)) return;
    memcpy(row, image->data + (y * image->span), (size_t)image->xsize);
}

void get_column(Pixel *col, Image *image, int x) {
    if (!image) return;
    if ((x < 0) || (x >= image->xsize)) return;
    int d = image->span;
    Pixel *p = image->data + x;
    for (int i = image->ysize; i-- > 0; p += d) *col++ = *p;
}

Pixel put_pixel(Image *image, int x, int y, Pixel data) {
    if (!image) return 0;
    if ((x < 0) || (x >= image->xsize) || (y < 0) || (y >= image->ysize)) return 0;
    image->data[y * image->span + x] = data;
    return data;
}

/* filter support constants and filter functions (we expose box and triangle) */
#define box_support      (0.5)
#define triangle_support (1.0)

double box_filter(double t) {
    if ((t > -0.5) && (t <= 0.5)) return 1.0;
    return 0.0;
}

double triangle_filter(double t) {
    if (t < 0.0) t = -t;
    if (t < 1.0) return 1.0 - t;
    return 0.0;
}

/* Contribution lists used by zoom() */
typedef struct {
    int pixel;
    double weight;
} CONTRIB;

typedef struct {
    int n;
    CONTRIB *p;
} CLIST;

static CLIST *contrib = NULL; /* allocated by zoom() */

/* zoom() implementation (two-pass separable resampler) */
void zoom(Image *dst, Image *src, double (*filterf)(double), double fwidth) {
    if (!dst || !src || !filterf) return;
    Image *tmp;
    double xscale = (double)dst->xsize / (double)src->xsize;
    double yscale = (double)dst->ysize / (double)src->ysize;

    /* intermediate image */
    tmp = new_image(dst->xsize, src->ysize);
    if (!tmp) return;

    /* horizontal contribution lists */
    contrib = (CLIST *)calloc((size_t)dst->xsize, sizeof(CLIST));
    if (!contrib) { free_image(tmp); return; }

    if (xscale < 1.0) {
        double width = fwidth / xscale;
        double fscale = 1.0 / xscale;
        for (int i = 0; i < dst->xsize; ++i) {
            contrib[i].n = 0;
            int pcount = (int)(width * 2 + 1);
            contrib[i].p = (CONTRIB *)calloc((size_t)pcount, sizeof(CONTRIB));
            double center = (double)i / xscale;
            int left = (int)ceil(center - width);
            int right = (int)floor(center + width);
            for (int j = left; j <= right; ++j) {
                double weight = center - (double)j;
                weight = (*filterf)(weight / fscale) / fscale;
                int n;
                if (j < 0) n = -j;
                else if (j >= src->xsize) n = (src->xsize - j) + src->xsize - 1;
                else n = j;
                int k = contrib[i].n++;
                contrib[i].p[k].pixel = n;
                contrib[i].p[k].weight = weight;
            }
        }
    } else {
        for (int i = 0; i < dst->xsize; ++i) {
            contrib[i].n = 0;
            int pcount = (int)(fwidth * 2 + 1);
            contrib[i].p = (CONTRIB *)calloc((size_t)pcount, sizeof(CONTRIB));
            double center = (double)i / xscale;
            int left = (int)ceil(center - fwidth);
            int right = (int)floor(center + fwidth);
            for (int j = left; j <= right; ++j) {
                double weight = center - (double)j;
                weight = (*filterf)(weight);
                int n;
                if (j < 0) n = -j;
                else if (j >= src->xsize) n = (src->xsize - j) + src->xsize - 1;
                else n = j;
                int k = contrib[i].n++;
                contrib[i].p[k].pixel = n;
                contrib[i].p[k].weight = weight;
            }
        }
    }

    /* apply horizontal filter (src -> tmp) */
    Pixel *raster = (Pixel *)calloc((size_t)src->xsize, sizeof(Pixel));
    if (!raster) { /* cleanup */ for (int i=0;i<dst->xsize;++i) free(contrib[i].p); free(contrib); free_image(tmp); return; }
    for (int k = 0; k < tmp->ysize; ++k) {
        get_row(raster, src, k);
        for (int i = 0; i < tmp->xsize; ++i) {
            double weight = 0.0;
            for (int j = 0; j < contrib[i].n; ++j) {
                weight += raster[contrib[i].p[j].pixel] * contrib[i].p[j].weight;
            }
            put_pixel(tmp, i, k, (Pixel)CLAMP(weight, 0, 255));
        }
    }
    free(raster);

    /* free horizontal contrib arrays */
    for (int i = 0; i < tmp->xsize; ++i) free(contrib[i].p);
    free(contrib);

    /* vertical contribution lists */
    contrib = (CLIST *)calloc((size_t)dst->ysize, sizeof(CLIST));
    if (!contrib) { free_image(tmp); return; }

    if (yscale < 1.0) {
        double width = fwidth / yscale;
        double fscale = 1.0 / yscale;
        for (int i = 0; i < dst->ysize; ++i) {
            contrib[i].n = 0;
            int pcount = (int)(width * 2 + 1);
            contrib[i].p = (CONTRIB *)calloc((size_t)pcount, sizeof(CONTRIB));
            double center = (double)i / yscale;
            int left = (int)ceil(center - width);
            int right = (int)floor(center + width);
            for (int j = left; j <= right; ++j) {
                double weight = center - (double)j;
                weight = (*filterf)(weight / fscale) / fscale;
                int n;
                if (j < 0) n = -j;
                else if (j >= tmp->ysize) n = (tmp->ysize - j) + tmp->ysize - 1;
                else n = j;
                int k = contrib[i].n++;
                contrib[i].p[k].pixel = n;
                contrib[i].p[k].weight = weight;
            }
        }
    } else {
        for (int i = 0; i < dst->ysize; ++i) {
            contrib[i].n = 0;
            int pcount = (int)(fwidth * 2 + 1);
            contrib[i].p = (CONTRIB *)calloc((size_t)pcount, sizeof(CONTRIB));
            double center = (double)i / yscale;
            int left = (int)ceil(center - fwidth);
            int right = (int)floor(center + fwidth);
            for (int j = left; j <= right; ++j) {
                double weight = center - (double)j;
                weight = (*filterf)(weight);
                int n;
                if (j < 0) n = -j;
                else if (j >= tmp->ysize) n = (tmp->ysize - j) + tmp->ysize - 1;
                else n = j;
                int k = contrib[i].n++;
                contrib[i].p[k].pixel = n;
                contrib[i].p[k].weight = weight;
            }
        }
    }

    /* apply vertical filter (tmp -> dst) */
    raster = (Pixel *)calloc((size_t)tmp->ysize, sizeof(Pixel));
    if (!raster) { for (int i=0;i<dst->ysize;++i) free(contrib[i].p); free(contrib); free_image(tmp); return; }
    for (int k = 0; k < dst->xsize; ++k) {
        get_column(raster, tmp, k);
        for (int i = 0; i < dst->ysize; ++i) {
            double weight = 0.0;
            for (int j = 0; j < contrib[i].n; ++j) {
                weight += raster[contrib[i].p[j].pixel] * contrib[i].p[j].weight;
            }
            put_pixel(dst, k, i, (Pixel)CLAMP(weight, 0, 255));
        }
    }
    free(raster);

    /* free vertical contrib arrays */
    for (int i = 0; i < dst->ysize; ++i) free(contrib[i].p);
    free(contrib);

    free_image(tmp);
}
/* --- End minimal Graphics Gems code --- */

// Helper to write scaled PNG
int write_png_scale(char const *filename, int x, int y, int comp, const void *data, int stride_bytes, float scale, int isLinear) {
    if (!filename || !data || x <= 0 || y <= 0 || comp <= 0) return -1;
    if (!(scale > 0.0f)) {
        fprintf(stderr, "write_png_scale: invalid scale %g\n", (double)scale);
        return -1;
    }

    // If scale == 1.0, call stbi_write_png directly using provided stride
    if (fabsf(scale - 1.0f) < 1e-9f) {
        // Note: stbi_write_png expects row stride in bytes; pass stride_bytes
        int ok = stbi_write_png(filename, x, y, comp, data, stride_bytes);
        return ok ? 0 : -1;
    }

    // Compute destination size: round to nearest integer
    int dst_w = (int)floor((double)x * (double)scale + 0.5);
    int dst_h = (int)floor((double)y * (double)scale + 0.5);

    if (dst_w <= 0 || dst_h <= 0) {
        fprintf(stderr, "write_png_scale: resulting dimensions invalid %d x %d\n", dst_w, dst_h);
        return -1;
    }

    // Allocate interleaved destination buffer: row stride = dst_w * comp
    size_t dst_row_bytes = (size_t)dst_w * (size_t)comp;
    size_t dst_buf_size = dst_row_bytes * (size_t)dst_h;
    uint8_t *dst_buf = (uint8_t *)malloc(dst_buf_size);
    if (!dst_buf) {
        fprintf(stderr, "write_png_scale: out of memory\n");
        return -1;
    }
    memset(dst_buf, 0, dst_buf_size);

    // Filter selection
    double (*filterf)() = NULL;
    double fwidth = 0.0;
    if (isLinear) {
        // tent / triangle (bilinear-like)
        filterf = triangle_filter;
        fwidth = 1.0; // triangle_support
    } else {
        // box (nearest / pixel)
        filterf = box_filter;
        fwidth = 0.5; // box_support
    }

    // For each component, extract channel, resample with zoom(), then interleave
    for (int c = 0; c < comp; ++c) {
        // build source single-channel buffer (x * y), copying from interleaved 'data' with stride
        Pixel *src_pixels = (Pixel *)malloc((size_t)x * (size_t)y);
        if (!src_pixels) {
            fprintf(stderr, "write_png_scale: out of memory (src_pixels)\n");
            free(dst_buf);
            return -1;
        }

        // Copy channel c into src_pixels row by row, using stride_bytes
        const uint8_t *in_row_ptr = (const uint8_t *)data;
        for (int row = 0; row < y; ++row) {
            const uint8_t *p = in_row_ptr + (size_t)row * (size_t)stride_bytes;
            for (int col = 0; col < x; ++col) {
                src_pixels[row * x + col] = p[col * comp + c];
            }
        }

        // Create Image wrappers
        Image src_img;
        src_img.xsize = x;
        src_img.ysize = y;
        src_img.span  = x;
        src_img.data  = src_pixels;

        Image *dst_img = new_image(dst_w, dst_h);
        if (!dst_img) {
            fprintf(stderr, "write_png_scale: new_image failed\n");
            free(src_pixels);
            free(dst_buf);
            return -1;
        }

        // Call Graphics Gems zoom() to resample
        zoom(dst_img, &src_img, filterf, fwidth);

        // Copy dst_img->data into interleaved dst_buf for component c
        for (int row = 0; row < dst_h; ++row) {
            uint8_t *out_row = dst_buf + (size_t)row * dst_row_bytes;
            Pixel *src_row = dst_img->data + (size_t)row * dst_img->span;
            for (int col = 0; col < dst_w; ++col) {
                out_row[col * comp + c] = src_row[col];
            }
        }

        // cleanup
        free_image(dst_img); // frees dst_img->data and the struct inside new_image/free_image
        free(src_pixels);
    }

    // Finally, write PNG
    int ok = stbi_write_png(filename, dst_w, dst_h, comp, dst_buf, (int)dst_row_bytes);
    free(dst_buf);
    if (!ok) {
        fprintf(stderr, "write_png_scale: stbi_write_png failed\n");
        return -1;
    }
    return ok;
}


// --------------------------- Types -------------------------------------
// Image-wide options (defaults applied unless per-tile flags override)
typedef struct {
	int isRadiological;		 // 1 = radiological (default), 0 = neurological
	int isLabel;			 // 1 = draw labels by default
	int interpolationOrder; //0=nearest, 1= linear, [future 2=Mitchell]
	float scale; // zoom factor
	float overlayRGBA[4];	 // overlay color + opacity (0..1) default {1,0,0,0.5}
	float baseRGBA[4];		 // base tint + opacity (alpha currently unused for base) default {1,1,1,1}
	float backgroundRGBA[4]; // background color + alpha default {0,0,0,0.5}
} image_options_t;

// tile description produced by parsing the argv
typedef struct {
	char axis;	   // 'x','y','z' or 'N' for newline/separator
	long long value; // numeric slice index (signed integer index)
	size_t w;	   // tile pixel width (computed from nim dims)
	size_t h;	   // tile pixel height (computed from nim dims)
	size_t x_off;  // x offset in row (filled during layout)
	size_t row;	   // row index (filled during layout)
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
																		const nifti_image *nim, char axis, double frac) {
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
	if (axis == 'x') value = clamp_ll(value, 0, nx - 1);
	if (axis == 'y') value = clamp_ll(value, 0, ny - 1);
	if (axis == 'z') value = clamp_ll(value, 0, nz - 1);

	tiles[n].axis = axis;
	tiles[n].value = value;
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
		if (!strcmp(tok, "-o")) { // overlayRGBA: -o R G B A
			if (i + 4 >= argc) { fprintf(stderr, "parse_args: -o requires 4 numeric args\n"); free(tiles); return -1; }
			if (!is_number_token(argv[i + 1]) || !is_number_token(argv[i + 2]) || !is_number_token(argv[i + 3]) || !is_number_token(argv[i + 4])) {
				fprintf(stderr, "parse_args: -o expects numeric args\n"); free(tiles); return -1;
			}
			opts->overlayRGBA[0] = fminf(fmaxf((float)strtod(argv[i + 1], NULL), 0.0f), 1.0f);
			opts->overlayRGBA[1] = fminf(fmaxf((float)strtod(argv[i + 2], NULL), 0.0f), 1.0f);
			opts->overlayRGBA[2] = fminf(fmaxf((float)strtod(argv[i + 3], NULL), 0.0f), 1.0f);
			opts->overlayRGBA[3] = fminf(fmaxf((float)strtod(argv[i + 4], NULL), 0.0f), 1.0f);
			i += 5;
			continue;
		}
		if (!strcmp(tok, "-c")) { // baseRGBA: -c R G B A
			if (i + 4 >= argc) { fprintf(stderr, "parse_args: -c requires 4 numeric args\n"); free(tiles); return -1; }
			if (!is_number_token(argv[i + 1]) || !is_number_token(argv[i + 2]) || !is_number_token(argv[i + 3]) || !is_number_token(argv[i + 4])) {
				fprintf(stderr, "parse_args: -c expects numeric args\n"); free(tiles); return -1;
			}
			opts->baseRGBA[0] = fminf(fmaxf((float)strtod(argv[i + 1], NULL), 0.0f), 1.0f);
			opts->baseRGBA[1] = fminf(fmaxf((float)strtod(argv[i + 2], NULL), 0.0f), 1.0f);
			opts->baseRGBA[2] = fminf(fmaxf((float)strtod(argv[i + 3], NULL), 0.0f), 1.0f);
			opts->baseRGBA[3] = fminf(fmaxf((float)strtod(argv[i + 4], NULL), 0.0f), 1.0f);
			i += 5;
			continue;
		}
		if (!strcmp(tok, "-b")) { // backgroundRGBA: -b R G B A
			if (i + 4 >= argc) { fprintf(stderr, "parse_args: -b requires 4 numeric args\n"); free(tiles); return -1; }
			if (!is_number_token(argv[i + 1]) || !is_number_token(argv[i + 2]) || !is_number_token(argv[i + 3]) || !is_number_token(argv[i + 4])) {
				fprintf(stderr, "parse_args: -b expects numeric args\n"); free(tiles); return -1;
			}
			opts->backgroundRGBA[0] = fminf(fmaxf((float)strtod(argv[i + 1], NULL), 0.0f), 1.0f);
			opts->backgroundRGBA[1] = fminf(fmaxf((float)strtod(argv[i + 2], NULL), 0.0f), 1.0f);
			opts->backgroundRGBA[2] = fminf(fmaxf((float)strtod(argv[i + 3], NULL), 0.0f), 1.0f);
			opts->backgroundRGBA[3] = fminf(fmaxf((float)strtod(argv[i + 4], NULL), 0.0f), 1.0f);
			i += 5;
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
			if (i + 1 < argc && is_number_token(argv[i+1])) {
				double v = strtod(argv[i+1], NULL);
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

		// ---------- Tile-related tokens ----------
		if (tok[0] == '-' && tok[1] != '\0') {
			char flag = tok[1];

			if (flag == 'r') {
				// row separator
				if (append_separator(&tiles, &n, &cap) != 0) { free(tiles); return -1; }
				++i;
				continue;
			}

			// ALIASES: -a and -m expand into tile lists
			if (flag == 'a') {
				if (append_tile_for_axis_value(&tiles, &n, &cap, nim, 'x', 0.5) != 0) { free(tiles); return -1; }
				if (append_tile_for_axis_value(&tiles, &n, &cap, nim, 'y', 0.5) != 0) { free(tiles); return -1; }
				if (append_tile_for_axis_value(&tiles, &n, &cap, nim, 'z', 0.5) != 0) { free(tiles); return -1; }
				++i;
				continue;
			}
			if (flag == 'm') {
				double vals[3] = {0.25, 0.5, 0.75};
				for (int k = 0; k < 3; ++k) {
					if (append_tile_for_axis_value(&tiles, &n, &cap, nim, 'x', vals[k]) != 0) { free(tiles); return -1; }
				}
				if (append_separator(&tiles, &n, &cap) != 0) { free(tiles); return -1; }
				for (int k = 0; k < 3; ++k) {
					if (append_tile_for_axis_value(&tiles, &n, &cap, nim, 'y', vals[k]) != 0) { free(tiles); return -1; }
				}
				if (append_separator(&tiles, &n, &cap) != 0) { free(tiles); return -1; }
				for (int k = 0; k < 3; ++k) {
					if (append_tile_for_axis_value(&tiles, &n, &cap, nim, 'z', vals[k]) != 0) { free(tiles); return -1; }
				}
				++i;
				continue;
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
				if (append_tile_for_axis_value(&tiles, &n, &cap, nim, flag, val) != 0) { free(tiles); return -1; }
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

// Returns tinted RGB for base nim using baseRGBA tint (RGB components 0..1)
static inline void get_base_rgb(float v, float cal_min, float cal_max,
								const float baseRGBA[4],
								uint8_t *r, uint8_t *g, uint8_t *b) {
	uint8_t u = float_to_u8(v, cal_min, cal_max); // 0..255
	float rr = baseRGBA[0] * u;
	float gg = baseRGBA[1] * u;
	float bb = baseRGBA[2] * u;
	if (rr < 0.0f)
		rr = 0.0f;
	if (rr > 255.0f)
		rr = 255.0f;
	if (gg < 0.0f)
		gg = 0.0f;
	if (gg > 255.0f)
		gg = 255.0f;
	if (bb < 0.0f)
		bb = 0.0f;
	if (bb > 255.0f)
		bb = 255.0f;
	*r = (uint8_t)(rr + 0.5f);
	*g = (uint8_t)(gg + 0.5f);
	*b = (uint8_t)(bb + 0.5f);
}

// Returns blended RGB (base tinted + overlay) using overlayRGBA (overlay color + alpha).
// If overlay intensity is zero, returns base directly.
static inline void get_blended_rgb(float v_base, float cal_min, float cal_max,
								   float v_overlay, float cal_min2, float cal_max2,
								   const float baseRGBA[4],
								   const float overlayRGBA[4],
								   uint8_t *r, uint8_t *g, uint8_t *b) {
	// compute base tinted color first
	uint8_t ub = float_to_u8(v_base, cal_min, cal_max);

	float base_r = baseRGBA[0] * ub;
	float base_g = baseRGBA[1] * ub;
	float base_b = baseRGBA[2] * ub;

	// overlay intensity (0..255)
	uint8_t uov = float_to_u8(v_overlay, cal_min2, cal_max2);
	if (uov == 0) {
		// no overlay contribution — return base tinted color
		if (base_r < 0.0f)
			base_r = 0.0f;
		if (base_r > 255.0f)
			base_r = 255.0f;
		if (base_g < 0.0f)
			base_g = 0.0f;
		if (base_g > 255.0f)
			base_g = 255.0f;
		if (base_b < 0.0f)
			base_b = 0.0f;
		if (base_b > 255.0f)
			base_b = 255.0f;
		*r = (uint8_t)(base_r + 0.5f);
		*g = (uint8_t)(base_g + 0.5f);
		*b = (uint8_t)(base_b + 0.5f);
		return;
	}

	float t2 = (float)uov / 255.0f; // overlay intensity 0..1
	float alpha = overlayRGBA[3];	// overlay opacity 0..1

	// overlay color scaled by its intensity
	float overlay_r = overlayRGBA[0] * t2 * 255.0f;
	float overlay_g = overlayRGBA[1] * t2 * 255.0f;
	float overlay_b = overlayRGBA[2] * t2 * 255.0f;

	// composite: out = (1 - alpha) * base + alpha * overlay
	float out_r = (1.0f - alpha) * base_r + alpha * overlay_r;
	float out_g = (1.0f - alpha) * base_g + alpha * overlay_g;
	float out_b = (1.0f - alpha) * base_b + alpha * overlay_b;

	if (out_r < 0.0f)
		out_r = 0.0f;
	if (out_r > 255.0f)
		out_r = 255.0f;
	if (out_g < 0.0f)
		out_g = 0.0f;
	if (out_g > 255.0f)
		out_g = 255.0f;
	if (out_b < 0.0f)
		out_b = 0.0f;
	if (out_b > 255.0f)
		out_b = 255.0f;

	*r = (uint8_t)(out_r + 0.5f);
	*g = (uint8_t)(out_g + 0.5f);
	*b = (uint8_t)(out_b + 0.5f);
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

// Process one voxel (data_idx) and write the pixel at dst_x,dst_y
static inline void process_and_write_pixel(uint8_t *rgb, size_t row_bytes, size_t W, size_t H,
										  size_t dst_x, size_t dst_y,
										  const float *data, const float *data2,
										  size_t data_idx,
										  int have_nim2,
										  float cal_min, float cal_max,
										  float cal_min2, float cal_max2,
										  const float baseRGBA[4], const float overlayRGBA[4]) {
	uint8_t r, g, b;
	if (have_nim2) {
		float v = data[data_idx];
		float v2 = data2[data_idx];
		get_blended_rgb(v, cal_min, cal_max, v2, cal_min2, cal_max2, baseRGBA, overlayRGBA, &r, &g, &b);
	} else {
		float v = data[data_idx];
		get_base_rgb(v, cal_min, cal_max, baseRGBA, &r, &g, &b);
	}
	write_pixel_rgb(rgb, row_bytes, W, H, dst_x, dst_y, r, g, b);
}

// ------------------- Stage 3: render_and_write_png ----------------------
// Render all tiles into an RGB buffer and write PNG using stbi_write_png.
// - nim: nifti image (assumed float data). Uses nim->cal_min and nim->cal_max if present, otherwise computes min/max.
// - tiles: tiles array with tile->row and tile->x_off filled by layout_tiles
// - n_tiles: number of tiles
// - layout: layout_t produced by layout_tiles
// - out_png_path: output filename
//
// Returns 0 on success, non-zero on failure.
static int render_and_write_png(const nifti_image *nim,
								nifti_image *nim2,
								const tile_t *tiles,
								size_t n_tiles,
								const layout_t *layout,
								const char *out_png_path,
								const image_options_t *opts) {
	if (!nim || !tiles || n_tiles == 0 || !layout || !out_png_path || !opts)
		return -1;

	// channels: RGB
	const int channels = 3;
	size_t W = layout->total_width;
	size_t H = layout->total_height;
	if (W == 0 || H == 0) { fprintf(stderr, "render: zero output dimension\n"); return -1; }

	// determine calibration min/max (use nifti cal_min/cal_max if present; otherwise compute from data)
	if (!(isfinite(nim->cal_min) && isfinite(nim->cal_max) && nim->cal_min < nim->cal_max)) {
		fprintf(stderr, "render: nim cal_min and cal_max must be finite and cal_min < cal_max (%g %g)\n", nim->cal_min, nim->cal_max);
		return -1;
	}
	float cal_min = nim->cal_min;
	float cal_max = nim->cal_max;

	/* initialize overlay cal range to main image's range by default */
	float cal_min2 = nim->cal_min;
	float cal_max2 = nim->cal_max;

	int have_nim2 = 0;
	float *data2 = NULL;
	if (nim2) {
		if (!(isfinite(nim2->cal_min) && isfinite(nim2->cal_max) && nim2->cal_min < nim2->cal_max)) {
			fprintf(stderr, "render: overlay (nim2) cal_min and cal_max must be finite and cal_min < cal_max (%g %g)\n", nim2->cal_min, nim2->cal_max);
			return -1;
		}
		if ((nim2->nx != nim->nx) || (nim2->ny != nim->ny) || (nim2->nz != nim->nz)) {
			fprintf(stderr, "render: overlay dimensions do not match\n");
			return -1;
		}
		if (!nim2->data) { fprintf(stderr, "render: overlay (nim2) data pointer is NULL\n"); return -1; }
		cal_min2 = nim2->cal_min;
		cal_max2 = nim2->cal_max;
		data2 = (float *)nim2->data;
		have_nim2 = 1;
	}

	// allocate RGB buffer: row-major, H rows of W pixels, each pixel channels bytes
	size_t row_bytes = W * channels;
	size_t buf_size = row_bytes * H;
	uint8_t *rgb = (uint8_t *)calloc(1, buf_size);
	if (!rgb) { fprintf(stderr, "render: out of memory\n"); return -1; }

	// Fill the image buffer with the background color from opts->backgroundRGBA
	uint8_t bg_r = (uint8_t)(opts->backgroundRGBA[0] * 255.0f + 0.5f);
	uint8_t bg_g = (uint8_t)(opts->backgroundRGBA[1] * 255.0f + 0.5f);
	uint8_t bg_b = (uint8_t)(opts->backgroundRGBA[2] * 255.0f + 0.5f);
	for (size_t y = 0; y < H; ++y) {
		uint8_t *row = rgb + y * row_bytes;
		for (size_t x = 0; x < W; ++x) {
			size_t idx = x * 3;
			row[idx + 0] = bg_r;
			row[idx + 1] = bg_g;
			row[idx + 2] = bg_b;
		}
	}

	// Inspect data as float* for now. If your images are another datatype, dispatch appropriately.
	float *data = (float *)nim->data;
	if (!data) { free(rgb); fprintf(stderr, "render: nim->data is NULL\n"); return -1; }

	// axis mapping
	long long nx = nim->nx, ny = nim->ny, nz = nim->nz;
	size_t slice_size = (size_t)nx * (size_t)ny;

	// iterate tiles and rasterize each tile's slice into the appropriate rectangle
	for (size_t ti = 0; ti < n_tiles; ++ti) {
		const tile_t *t = &tiles[ti];
		if (t->axis == 'N') continue; // separator tile

		// where to place this tile
		size_t row = t->row;
		if (row == (size_t)-1 || row >= layout->rows) continue;
		size_t base_x = t->x_off;
		size_t base_y = layout->row_y_offsets[row];

		// slice index (already stored as integer)
		long long slice_idx = clamp_ll(t->value, 0, (t->axis == 'x' ? nx - 1 : (t->axis == 'y' ? ny - 1 : nz - 1)));

		// iterate over tile pixels, mapping to image voxel indices
		for (size_t yy = 0; yy < t->h; ++yy) {
			for (size_t xx = 0; xx < t->w; ++xx) {
				long long xi, yi, zi;
				map_tile_coords(t->axis, slice_idx, xx, yy, nx, ny, nz, &xi, &yi, &zi);

				// bounds check (only upper bounds necessary because xi/yi/zi are unsigned-constructed)
				if (xi < 0 || xi >= nx || yi < 0 || yi >= ny || zi < 0 || zi >= nz)
					continue;

				size_t data_idx = (size_t)xi + (size_t)yi * (size_t)nx + (size_t)zi * (size_t)slice_size;

				// compute destination pixel coords (radiological flips horizontally)
				size_t dst_x = opts->isRadiological ? (base_x + (t->w - 1 - xx)) : (base_x + xx);
				size_t dst_y = base_y + (t->h - 1 - yy);

				// process pixel and write
				process_and_write_pixel(rgb, row_bytes, W, H, dst_x, dst_y, data, data2, data_idx, have_nim2,
										cal_min, cal_max, cal_min2, cal_max2, opts->baseRGBA, opts->overlayRGBA);
			}
		}

		// optionally burn labels L/R at left/right margin if requested
		if ((opts->isLabel) && (t->axis != 'x')) {
			size_t inset = 2;
			size_t lab_y = base_y + (t->h - 2 - 5 + 1); // base_y + t->h - 6
			if (lab_y + 5 <= H) {
				char left_label = opts->isRadiological ? 'R' : 'L';
				burn_label_rgb(rgb, W, H, base_x + inset, lab_y, left_label);
				char right_label = opts->isRadiological ? 'L' : 'R';
				if (t->w >= 7) {
					size_t px = base_x + t->w - inset - 5; // 5 px label width
					burn_label_rgb(rgb, W, H, px, lab_y, right_label);
				}
			}
		}
	}

	// write PNG
	int ok = 0;

	if ((opts->scale > 0.0) && (opts->scale != 1.0))
		ok = write_png_scale(out_png_path, (int)W, (int)H, channels, rgb, (int)row_bytes, opts->scale, opts->interpolationOrder);
	else
		ok = stbi_write_png(out_png_path, (int)W, (int)H, channels, rgb, (int)row_bytes);
	free(rgb);
	if (!ok) {
		fprintf(stderr, "render: stbi_write_png failed\n");
		return -1;
	}
	return 0;
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
	opts.interpolationOrder = 1;
	opts.scale = 1.0;
	opts.overlayRGBA[0] = 1.0f; opts.overlayRGBA[1] = 0.0f; opts.overlayRGBA[2] = 0.0f; opts.overlayRGBA[3] = 0.5f;
	opts.baseRGBA[0] = 1.0f; opts.baseRGBA[1] = 1.0f; opts.baseRGBA[2] = 1.0f; opts.baseRGBA[3] = 1.0f;
	opts.backgroundRGBA[0] = 0.0f; opts.backgroundRGBA[1] = 0.0f; opts.backgroundRGBA[2] = 0.0f; opts.backgroundRGBA[3] = 1.0f;

	if (parse_args_build_tiles(nim, argc, argv, &tiles, &n_tiles, &opts) != 0) {
		fprintf(stderr, "nim2png_refactored: failed to parse args\n");
		return -1;
	}

	// Stage 2: layout tiles -> layout_t
	layout_t layout;
	memset(&layout, 0, sizeof(layout));
	if (layout_tiles(tiles, n_tiles, &layout) != 0) {
		fprintf(stderr, "nim2png_refactored: failed to layout tiles\n");
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
