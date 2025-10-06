#ifndef NIIMATH_FILTER_H
#define NIIMATH_FILTER_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef unsigned char Pixel;

/* Image: interleaved row-major storage with 'components' per pixel.
   - xsize, ysize: dimensions in pixels
   - components: 1..4 (or more, but this code assumes small component counts)
   - data: pointer to allocated buffer of size xsize * ysize * components
   - span: bytes per row = xsize * components
*/
typedef struct {
	int xsize;
	int ysize;
	int components;
	Pixel *data;
	int span;
} Image;

/* create/destroy Image */
Image *new_image(int xsize, int ysize, int components);
void free_image(Image *image);

/* basic pixel access helpers (operate on interleaved images) */
Pixel get_pixel(const Image *image, int x, int y, int c);
void get_row(Pixel *row, const Image *image, int y, int c);
void get_column(Pixel *col, const Image *image, int x, int c);
void put_pixel(Image *image, int x, int y, int c, Pixel val);

/* filter functions */
double box_filter(double t);
double triangle_filter(double t);

/* zoom resampler (multi-component aware)
 * - dst, src: Image pointers (may have components >= 1)
 * - filterf: filter function (box_filter or triangle_filter)
 * - fwidth: filter support (0.5 for box, 1.0 for triangle)
 *
 * Returns 0 on success, -1 on allocation error or invalid args.
 */
int zoom(Image *dst, const Image *src, double (*filterf)(double), double fwidth);

#ifdef __cplusplus
}
#endif

#endif // NIIMATH_FILTER_H
