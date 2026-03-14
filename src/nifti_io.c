/* nifti_io.c - Minimal NIfTI I/O for niimath
 *
 * Replaces niftilib (nifti2_io.c) and znzlib (znzlib.c).
 * Public domain code based on work by Robert W Cox, Mark Jenkinson,
 * Rick Reynolds, and Chris Rorden.
 */

#include "nifti_io.h"
#include <sys/types.h>
#include <sys/stat.h>

#ifdef _WIN32
  #include <windows.h>
  #include <io.h>
  #include <fcntl.h>
#else
  #include <unistd.h>
  #include <poll.h>
#endif

#ifdef HAVE_ZLIB
  #include <zlib.h>
#endif
#ifdef HAVE_ZSTD
  #include <zstd.h>
#endif

/* compile-time flags from Makefile */
/* FSLSTYLE, PIGZ, REJECT_COMPLEX are set externally */

/*========== Internal file I/O abstraction (replaces znzlib) ==========*/

/* Compression types for NIIFILE */
#define NII_COMPRESS_NONE 0
#define NII_COMPRESS_GZ   1
#define NII_COMPRESS_ZST  2

/* A minimal file handle that wraps FILE* and optionally gzFile.
   For zstd: decompresses to tmpfile on open, compresses from tmpfile on close.
   All intermediate I/O goes through the standard FILE* path. */
typedef struct {
   FILE  *nzfptr;
#ifdef HAVE_ZLIB
   gzFile zfptr;
#endif
   int    withz;        /* NII_COMPRESS_NONE, NII_COMPRESS_GZ, or NII_COMPRESS_ZST */
#ifdef HAVE_ZSTD
   char  *zst_fname;    /* target filename for zstd write (deferred compression) */
#endif
} *NIIFILE;

#define ZNZ_MAX_BLOCK_SIZE (1<<30)

/* Detect compression type from filename extension */
static int nii_compression_type(const char *fname)
{
   if (!fname) return NII_COMPRESS_NONE;
   int len = (int)strlen(fname);
#ifdef HAVE_ZSTD
   if (len >= 4 && (strcmp(fname + len - 4, ".zst") == 0 ||
                    strcmp(fname + len - 4, ".ZST") == 0))
      return NII_COMPRESS_ZST;
#endif
#ifdef HAVE_ZLIB
   if (len >= 3 && (strcmp(fname + len - 3, ".gz") == 0 ||
                    strcmp(fname + len - 3, ".GZ") == 0))
      return NII_COMPRESS_GZ;
#endif
   return NII_COMPRESS_NONE;
}

#ifdef HAVE_ZSTD
/* Decompress a .zst file into a tmpfile, return the tmpfile handle */
static FILE *zst_decompress_to_tmpfile(const char *path)
{
   FILE *fin = fopen(path, "rb");
   if (!fin) return NULL;
   fseek(fin, 0, SEEK_END);
   long csize = ftell(fin);
   fseek(fin, 0, SEEK_SET);
   if (csize <= 0) { fclose(fin); return NULL; }

   char *cbuf = (char *)malloc(csize);
   if (!cbuf) { fclose(fin); return NULL; }
   if (fread(cbuf, 1, csize, fin) != (size_t)csize) { free(cbuf); fclose(fin); return NULL; }
   fclose(fin);

   unsigned long long dsize = ZSTD_getFrameContentSize(cbuf, csize);
   if (dsize == ZSTD_CONTENTSIZE_ERROR) {
      fprintf(stderr, "** zstd: corrupt frame in '%s'\n", path);
      free(cbuf); return NULL;
   }
   if (dsize == ZSTD_CONTENTSIZE_UNKNOWN)
      dsize = (unsigned long long)csize * 8; /* estimate; will retry if too small */

   char *dbuf = (char *)malloc(dsize);
   if (!dbuf) { free(cbuf); return NULL; }

   size_t result = ZSTD_decompress(dbuf, dsize, cbuf, csize);
   /* If buffer was too small (unknown size case), retry with larger buffer */
   while (ZSTD_isError(result) && ZSTD_getErrorCode(result) == ZSTD_error_dstSize_tooSmall && dsize < (unsigned long long)csize * 256) {
      dsize *= 2;
      char *newbuf = (char *)realloc(dbuf, dsize);
      if (!newbuf) { free(dbuf); free(cbuf); return NULL; }
      dbuf = newbuf;
      result = ZSTD_decompress(dbuf, dsize, cbuf, csize);
   }
   free(cbuf);
   if (ZSTD_isError(result)) {
      fprintf(stderr, "** zstd decompression error: %s\n", ZSTD_getErrorName(result));
      free(dbuf); return NULL;
   }

   FILE *tmpf = tmpfile();
   if (!tmpf) { free(dbuf); return NULL; }
   fwrite(dbuf, 1, result, tmpf);
   free(dbuf);
   rewind(tmpf);
   return tmpf;
}

/* Compress tmpfile contents to a .zst file */
static int zst_compress_from_tmpfile(FILE *tmpf, const char *path)
{
   fseek(tmpf, 0, SEEK_END);
   long dsize = ftell(tmpf);
   fseek(tmpf, 0, SEEK_SET);
   if (dsize <= 0) return -1;

   char *dbuf = (char *)malloc(dsize);
   if (!dbuf) return -1;
   if (fread(dbuf, 1, dsize, tmpf) != (size_t)dsize) { free(dbuf); return -1; }

   size_t cbound = ZSTD_compressBound(dsize);
   char *cbuf = (char *)malloc(cbound);
   if (!cbuf) { free(dbuf); return -1; }

   size_t csize = ZSTD_compress(cbuf, cbound, dbuf, dsize, 3); /* level 3: fast */
   free(dbuf);
   if (ZSTD_isError(csize)) {
      fprintf(stderr, "** zstd compression error: %s\n", ZSTD_getErrorName(csize));
      free(cbuf); return -1;
   }

   FILE *fout = fopen(path, "wb");
   if (!fout) { free(cbuf); return -1; }
   fwrite(cbuf, 1, csize, fout);
   fclose(fout);
   free(cbuf);
   return 0;
}
#endif /* HAVE_ZSTD */

static NIIFILE nii_open(const char *path, const char *mode, int use_gz)
{
   NIIFILE f = (NIIFILE)calloc(1, sizeof(*f));
   if (!f) return NULL;

   /* stdout/stdin special case */
   if (strcmp(path, "-") == 0) {
      if (mode[0] == 'r') f->nzfptr = stdin;
      else                 f->nzfptr = stdout;
      return f;
   }

#ifdef HAVE_ZSTD
   if (nii_compression_type(path) == NII_COMPRESS_ZST) {
      f->withz = NII_COMPRESS_ZST;
      if (mode[0] == 'r') {
         f->nzfptr = zst_decompress_to_tmpfile(path);
         if (!f->nzfptr) { free(f); return NULL; }
      } else {
         /* Deferred write: buffer to tmpfile, compress on close */
         f->zst_fname = nifti_strdup(path);
         f->nzfptr = tmpfile();
         if (!f->nzfptr || !f->zst_fname) { free(f->zst_fname); free(f); return NULL; }
      }
      return f;
   }
#endif

#ifdef HAVE_ZLIB
   if (use_gz == 1) {
      f->withz = NII_COMPRESS_GZ;
      f->zfptr = gzopen(path, mode);
      if (!f->zfptr) { free(f); return NULL; }
      return f;
   }
#else
   (void)use_gz;
#endif

   f->nzfptr = fopen(path, mode);
   if (!f->nzfptr) { free(f); return NULL; }
   return f;
}

static int nii_close(NIIFILE *fp)
{
   int ret = 0;
   if (!fp || !*fp) return 0;
#ifdef HAVE_ZSTD
   if ((*fp)->withz == NII_COMPRESS_ZST && (*fp)->zst_fname) {
      /* Deferred zstd write: compress tmpfile contents to target */
      ret = zst_compress_from_tmpfile((*fp)->nzfptr, (*fp)->zst_fname);
      fclose((*fp)->nzfptr);
      free((*fp)->zst_fname);
      free(*fp); *fp = NULL;
      return ret;
   }
   if ((*fp)->withz == NII_COMPRESS_ZST) {
      /* zstd read: just close the tmpfile */
      if ((*fp)->nzfptr) fclose((*fp)->nzfptr);
      free(*fp); *fp = NULL;
      return 0;
   }
#endif
#ifdef HAVE_ZLIB
   if ((*fp)->zfptr) { ret = gzclose((*fp)->zfptr); }
   else
#endif
   if ((*fp)->nzfptr) {
      if ((*fp)->nzfptr != stdout && (*fp)->nzfptr != stdin)
         ret = fclose((*fp)->nzfptr);
   }
   free(*fp);
   *fp = NULL;
   return ret;
}

static size_t nii_read(void *buf, size_t size, size_t nmemb, NIIFILE f)
{
   if (!f) return 0;
#ifdef HAVE_ZLIB
   if (f->zfptr) {
      size_t remain = size * nmemb;
      char *cbuf = (char *)buf;
      while (remain > 0) {
         unsigned n2r = (remain < ZNZ_MAX_BLOCK_SIZE) ? (unsigned)remain : ZNZ_MAX_BLOCK_SIZE;
         int nr = gzread(f->zfptr, cbuf, n2r);
         if (nr < 0) return 0;
         remain -= nr;
         cbuf += nr;
         if (nr < (int)n2r) break;
      }
      return nmemb - remain / size;
   }
#endif
   return fread(buf, size, nmemb, f->nzfptr);
}

static size_t nii_write(const void *buf, size_t size, size_t nmemb, NIIFILE f)
{
   if (!f) return 0;
#ifdef HAVE_ZLIB
   if (f->zfptr) {
      size_t remain = size * nmemb;
      const char *cbuf = (const char *)buf;
      while (remain > 0) {
         unsigned n2w = (remain < ZNZ_MAX_BLOCK_SIZE) ? (unsigned)remain : ZNZ_MAX_BLOCK_SIZE;
         int nw = gzwrite(f->zfptr, cbuf, n2w);
         if (nw < 0) return (size_t)nw;
         remain -= nw;
         cbuf += nw;
         if (nw < (int)n2w) break;
      }
      return nmemb - remain / size;
   }
#endif
   return fwrite(buf, size, nmemb, f->nzfptr);
}

static long nii_seek(NIIFILE f, long offset, int whence)
{
   if (!f) return 0;
#ifdef HAVE_ZLIB
   if (f->zfptr) return (long)gzseek(f->zfptr, offset, whence);
#endif
   if (f->nzfptr == stdin) {
      if (whence == SEEK_SET && offset > 0) {
         long skipped = 0;
         char buf[512];
         while (skipped < offset) {
            size_t to_read = (offset - skipped < (long)sizeof(buf)) ? (size_t)(offset - skipped) : sizeof(buf);
            size_t bytes = fread(buf, 1, to_read, f->nzfptr);
            if (bytes == 0) return -1;
            skipped += bytes;
         }
         return 0;
      }
      return -1;
   }
   return fseek(f->nzfptr, offset, whence);
}

static long nii_tell(NIIFILE f)
{
   if (!f) return 0;
#ifdef HAVE_ZLIB
   if (f->zfptr) return (long)gztell(f->zfptr);
#endif
   return ftell(f->nzfptr);
}

static int nii_rewind(NIIFILE f)
{
   if (!f) return 0;
#ifdef HAVE_ZLIB
   if (f->zfptr) return (int)gzseek(f->zfptr, 0L, SEEK_SET);
#endif
   rewind(f->nzfptr);
   return 0;
}

/*========== Magic bytes ==========*/

static char nifti1_magic[4] = { 'n', '+', '1', '\0' };
static char nifti2_magic[8] = { 'n', '+', '2', '\0', '\r', '\n', '\032', '\n' };

/*========== Utility functions ==========*/

char *nifti_strdup(const char *str)
{
   char *dup;
   if (!str) return NULL;
   dup = (char *)malloc(strlen(str) + 1);
   if (dup) strcpy(dup, str);
   return dup;
}

int nifti_short_order(void)
{
   union { unsigned char bb[2]; short ss; } test;
   test.bb[0] = 1; test.bb[1] = 0;
   return (test.ss == 1) ? LSB_FIRST : MSB_FIRST;
}

int nifti_compiled_with_zlib(void)
{
#ifdef HAVE_ZLIB
   return 1;
#else
   return 0;
#endif
}

static int stdin_has_data(void)
{
#ifdef _WIN32
   DWORD bytesAvailable = 0;
   HANDLE hStdin = GetStdHandle(STD_INPUT_HANDLE);
   if (hStdin == INVALID_HANDLE_VALUE) return 0;
   if (!PeekNamedPipe(hStdin, NULL, 0, NULL, &bytesAvailable, NULL)) return 0;
   return bytesAvailable > 0;
#else
   struct pollfd pfd = { .fd = STDIN_FILENO, .events = POLLIN };
   int ret = poll(&pfd, 1, 0);
   return (ret > 0 && (pfd.revents & POLLIN));
#endif
}

static int64_t nifti_get_filesize(const char *pathname)
{
   struct stat buf;
   if (!pathname || !*pathname) return -1;
   if (stat(pathname, &buf) != 0) return -1;
   return buf.st_size;
}

void nifti_datatype_sizes(int datatype, int *nbyper, int *swapsize)
{
   int nb = 0, ss = 0;
   switch (datatype) {
      case DT_INT8: case DT_UINT8:    nb = 1; ss = 0; break;
      case DT_INT16: case DT_UINT16:  nb = 2; ss = 2; break;
      case DT_RGB24:                   nb = 3; ss = 0; break;
      case DT_RGBA32:                  nb = 4; ss = 0; break;
      case DT_INT32: case DT_UINT32:
      case DT_FLOAT32:                 nb = 4; ss = 4; break;
      case DT_COMPLEX64:              nb = 8; ss = 4; break;
      case DT_FLOAT64: case DT_INT64:
      case DT_UINT64:                  nb = 8; ss = 8; break;
      case DT_FLOAT128:               nb = 16; ss = 16; break;
      case DT_COMPLEX128:             nb = 16; ss = 8; break;
      case DT_COMPLEX256:             nb = 32; ss = 16; break;
      default: break;
   }
   ASSIF(nbyper, nb);
   ASSIF(swapsize, ss);
}

/*========== Byte swapping ==========*/

void nifti_swap_2bytes(int64_t n, void *ar)
{
   unsigned char *cp1 = (unsigned char *)ar, *cp2, tval;
   for (int64_t ii = 0; ii < n; ii++) {
      cp2 = cp1 + 1;
      tval = *cp1; *cp1 = *cp2; *cp2 = tval;
      cp1 += 2;
   }
}

void nifti_swap_4bytes(int64_t n, void *ar)
{
   unsigned char *cp0 = (unsigned char *)ar, *cp1, *cp2, tval;
   for (int64_t ii = 0; ii < n; ii++) {
      cp1 = cp0; cp2 = cp0 + 3;
      tval = *cp1; *cp1 = *cp2; *cp2 = tval;
      cp1++; cp2--;
      tval = *cp1; *cp1 = *cp2; *cp2 = tval;
      cp0 += 4;
   }
}

void nifti_swap_8bytes(int64_t n, void *ar)
{
   unsigned char *cp0 = (unsigned char *)ar, *cp1, *cp2, tval;
   for (int64_t ii = 0; ii < n; ii++) {
      cp1 = cp0; cp2 = cp0 + 7;
      while (cp2 > cp1) {
         tval = *cp1; *cp1 = *cp2; *cp2 = tval;
         cp1++; cp2--;
      }
      cp0 += 8;
   }
}

void nifti_swap_16bytes(int64_t n, void *ar)
{
   unsigned char *cp0 = (unsigned char *)ar, *cp1, *cp2, tval;
   for (int64_t ii = 0; ii < n; ii++) {
      cp1 = cp0; cp2 = cp0 + 15;
      while (cp2 > cp1) {
         tval = *cp1; *cp1 = *cp2; *cp2 = tval;
         cp1++; cp2--;
      }
      cp0 += 16;
   }
}

void nifti_swap_Nbytes(int64_t n, int siz, void *ar)
{
   switch (siz) {
      case 2:  nifti_swap_2bytes(n, ar); break;
      case 4:  nifti_swap_4bytes(n, ar); break;
      case 8:  nifti_swap_8bytes(n, ar); break;
      case 16: nifti_swap_16bytes(n, ar); break;
      default: fprintf(stderr, "** NIfTI: cannot swap in %d byte blocks\n", siz); break;
   }
}

static void swap_nifti1(nifti_1_header *h)
{
   if (!h) return;
   nifti_swap_4bytes(1, &h->sizeof_hdr);
   nifti_swap_4bytes(1, &h->extents);
   nifti_swap_2bytes(1, &h->session_error);
   nifti_swap_2bytes(8, h->dim);
   nifti_swap_4bytes(1, &h->intent_p1);
   nifti_swap_4bytes(1, &h->intent_p2);
   nifti_swap_4bytes(1, &h->intent_p3);
   nifti_swap_2bytes(1, &h->intent_code);
   nifti_swap_2bytes(1, &h->datatype);
   nifti_swap_2bytes(1, &h->bitpix);
   nifti_swap_2bytes(1, &h->slice_start);
   nifti_swap_4bytes(8, h->pixdim);
   nifti_swap_4bytes(1, &h->vox_offset);
   nifti_swap_4bytes(1, &h->scl_slope);
   nifti_swap_4bytes(1, &h->scl_inter);
   nifti_swap_2bytes(1, &h->slice_end);
   nifti_swap_4bytes(1, &h->cal_max);
   nifti_swap_4bytes(1, &h->cal_min);
   nifti_swap_4bytes(1, &h->slice_duration);
   nifti_swap_4bytes(1, &h->toffset);
   nifti_swap_4bytes(1, &h->glmax);
   nifti_swap_4bytes(1, &h->glmin);
   nifti_swap_2bytes(1, &h->qform_code);
   nifti_swap_2bytes(1, &h->sform_code);
   nifti_swap_4bytes(1, &h->quatern_b);
   nifti_swap_4bytes(1, &h->quatern_c);
   nifti_swap_4bytes(1, &h->quatern_d);
   nifti_swap_4bytes(1, &h->qoffset_x);
   nifti_swap_4bytes(1, &h->qoffset_y);
   nifti_swap_4bytes(1, &h->qoffset_z);
   nifti_swap_4bytes(4, h->srow_x);
   nifti_swap_4bytes(4, h->srow_y);
   nifti_swap_4bytes(4, h->srow_z);
}

static void swap_nifti2(nifti_2_header *h)
{
   if (!h) return;
   nifti_swap_4bytes(1, &h->sizeof_hdr);
   nifti_swap_2bytes(1, &h->datatype);
   nifti_swap_2bytes(1, &h->bitpix);
   nifti_swap_8bytes(8, h->dim);
   nifti_swap_8bytes(1, &h->intent_p1);
   nifti_swap_8bytes(1, &h->intent_p2);
   nifti_swap_8bytes(1, &h->intent_p3);
   nifti_swap_8bytes(8, h->pixdim);
   nifti_swap_8bytes(1, &h->vox_offset);
   nifti_swap_8bytes(1, &h->scl_slope);
   nifti_swap_8bytes(1, &h->scl_inter);
   nifti_swap_8bytes(1, &h->cal_max);
   nifti_swap_8bytes(1, &h->cal_min);
   nifti_swap_8bytes(1, &h->slice_duration);
   nifti_swap_8bytes(1, &h->toffset);
   nifti_swap_8bytes(1, &h->slice_start);
   nifti_swap_8bytes(1, &h->slice_end);
   nifti_swap_4bytes(1, &h->qform_code);
   nifti_swap_4bytes(1, &h->sform_code);
   nifti_swap_8bytes(1, &h->quatern_b);
   nifti_swap_8bytes(1, &h->quatern_c);
   nifti_swap_8bytes(1, &h->quatern_d);
   nifti_swap_8bytes(1, &h->qoffset_x);
   nifti_swap_8bytes(1, &h->qoffset_y);
   nifti_swap_8bytes(1, &h->qoffset_z);
   nifti_swap_8bytes(4, h->srow_x);
   nifti_swap_8bytes(4, h->srow_y);
   nifti_swap_8bytes(4, h->srow_z);
   nifti_swap_4bytes(1, &h->slice_code);
   nifti_swap_4bytes(1, &h->xyzt_units);
   nifti_swap_4bytes(1, &h->intent_code);
}

static void swap_as_analyze(nifti_analyze75 *h)
{
   if (!h) return;
   nifti_swap_4bytes(1, &h->sizeof_hdr);
   nifti_swap_4bytes(1, &h->extents);
   nifti_swap_2bytes(1, &h->session_error);
   nifti_swap_2bytes(8, h->dim);
   nifti_swap_2bytes(1, &h->unused8);  nifti_swap_2bytes(1, &h->unused9);
   nifti_swap_2bytes(1, &h->unused10); nifti_swap_2bytes(1, &h->unused11);
   nifti_swap_2bytes(1, &h->unused12); nifti_swap_2bytes(1, &h->unused13);
   nifti_swap_2bytes(1, &h->unused14);
   nifti_swap_2bytes(1, &h->datatype);
   nifti_swap_2bytes(1, &h->bitpix);
   nifti_swap_2bytes(1, &h->dim_un0);
   nifti_swap_4bytes(8, h->pixdim);
   nifti_swap_4bytes(1, &h->vox_offset);
   nifti_swap_4bytes(1, &h->funused1);
   nifti_swap_4bytes(1, &h->funused2);
   nifti_swap_4bytes(1, &h->funused3);
   nifti_swap_4bytes(1, &h->cal_max);
   nifti_swap_4bytes(1, &h->cal_min);
   nifti_swap_4bytes(1, &h->compressed);
   nifti_swap_4bytes(1, &h->verified);
   nifti_swap_4bytes(1, &h->glmax);
   nifti_swap_4bytes(1, &h->glmin);
   nifti_swap_4bytes(1, &h->views);
   nifti_swap_4bytes(1, &h->vols_added);
   nifti_swap_4bytes(1, &h->start_field);
   nifti_swap_4bytes(1, &h->field_skip);
   nifti_swap_4bytes(1, &h->omax);
   nifti_swap_4bytes(1, &h->omin);
   nifti_swap_4bytes(1, &h->smax);
   nifti_swap_4bytes(1, &h->smin);
}

static void swap_nifti_header(void *hdr, int ni_ver)
{
   if      (ni_ver == 0) swap_as_analyze((nifti_analyze75 *)hdr);
   else if (ni_ver == 1) swap_nifti1((nifti_1_header *)hdr);
   else if (ni_ver == 2) swap_nifti2((nifti_2_header *)hdr);
}

/*========== String helpers (case-insensitive extension matching) ==========*/

static int compare_strlist(const char *str, char **strlist, int len)
{
   for (int c = 0; c < len; c++)
      if (strlist[c] && !strcmp(str, strlist[c])) return c;
   return -1;
}

static int fileext_compare(const char *test_ext, const char *known_ext)
{
   char caps[8] = "";
   const int cmp = strcmp(test_ext, known_ext);
   if (cmp == 0) return 0;
   if (!test_ext || !known_ext) return cmp;
   size_t len = strlen(known_ext);
   if (len > 7) return cmp;
   for (size_t c = 0; c < len; c++) caps[c] = toupper((int)known_ext[c]);
   caps[len] = '\0';
   return strcmp(test_ext, caps);
}

static int fileext_n_compare(const char *test_ext, const char *known_ext, size_t maxlen)
{
   char caps[8] = "";
   const int cmp = strncmp(test_ext, known_ext, maxlen);
   if (cmp == 0) return 0;
   if (!test_ext || !known_ext) return cmp;
   size_t len = strlen(known_ext);
   if (len > maxlen) len = maxlen;
   if (len > 7) return cmp;
   for (size_t c = 0; c < len; c++) caps[c] = toupper((int)known_ext[c]);
   caps[len] = '\0';
   return strncmp(test_ext, caps, maxlen);
}

static int is_uppercase(const char *str)
{
   int hasupper = 0;
   if (!str || !*str) return 0;
   for (size_t c = 0; c < strlen(str); c++) {
      if (islower((int)str[c])) return 0;
      if (!hasupper && isupper((int)str[c])) hasupper = 1;
   }
   return hasupper;
}

static int is_mixedcase(const char *str)
{
   int hasupper = 0, haslower = 0;
   if (!str || !*str) return 0;
   for (size_t c = 0; c < strlen(str); c++) {
      if (!haslower && islower((int)str[c])) haslower = 1;
      if (!hasupper && isupper((int)str[c])) hasupper = 1;
      if (haslower && hasupper) return 1;
   }
   return 0;
}

static int make_lowercase(char *str)
{
   if (!str || !*str) return 0;
   for (size_t c = 0; c < strlen(str); c++)
      if (isupper((int)str[c])) str[c] = tolower((int)str[c]);
   return 0;
}

static int make_uppercase(char *str)
{
   if (!str || !*str) return 0;
   for (size_t c = 0; c < strlen(str); c++)
      if (islower((int)str[c])) str[c] = toupper((int)str[c]);
   return 0;
}

/*========== File name utilities ==========*/

/* Returns compression type: 0=none, 1=gzip, 2=zstd */
int nifti_is_gzfile(const char *fname)
{
   if (!fname) return 0;
   int len = (int)strlen(fname);
#ifdef HAVE_ZSTD
   if (len >= 4 && fileext_compare(fname + len - 4, ".zst") == 0) return 2;
#endif
   if (len < 3) return 0;
   if (fileext_compare(fname + len - 3, ".gz") == 0) {
#ifdef HAVE_ZLIB
      return 1;
#else
      fprintf(stderr, "** nifti_is_gzfile: recompile for compressed data '%s'\n", fname);
#endif
   }
   return 0;
}

char *nifti_find_file_extension(const char *name)
{
   const char *ext;
   char extcopy[10];
   int len;
   char extnii[10] = ".nii";
   char exthdr[10] = ".hdr";
   char extimg[10] = ".img";
   char extnia[10] = ".nia";
   char extgz[5]  = ".gz";
   char extzst[5] = ".zst";
   char *elist[4] = { NULL, NULL, NULL, NULL };

   elist[0] = extnii; elist[1] = exthdr; elist[2] = extimg; elist[3] = extnia;

   if (!name) return NULL;
   len = (int)strlen(name);
   if (len < 4) return NULL;

   ext = name + len - 4;
   strncpy(extcopy, ext, 9); extcopy[9] = '\0';
   make_lowercase(extcopy);  /* always allow uppercase extensions */

   if (compare_strlist(extcopy, elist, 4) >= 0) {
      if (is_mixedcase(ext)) return NULL;
      return (char *)ext;
   }

#ifdef HAVE_ZLIB
   if (len >= 7) {
      char gzlist[3][10];
      char *gzp[3];
      strcpy(gzlist[0], ".nii"); strcat(gzlist[0], extgz);
      strcpy(gzlist[1], ".hdr"); strcat(gzlist[1], extgz);
      strcpy(gzlist[2], ".img"); strcat(gzlist[2], extgz);
      gzp[0] = gzlist[0]; gzp[1] = gzlist[1]; gzp[2] = gzlist[2];
      ext = name + len - 7;
      strncpy(extcopy, ext, 9); extcopy[9] = '\0';
      make_lowercase(extcopy);
      if (compare_strlist(extcopy, gzp, 3) >= 0) {
         if (is_mixedcase(ext)) return NULL;
         return (char *)ext;
      }
   }
#endif

#ifdef HAVE_ZSTD
   if (len >= 8) {
      char zstlist[3][10];
      char *zstp[3];
      strcpy(zstlist[0], ".nii"); strcat(zstlist[0], extzst);
      strcpy(zstlist[1], ".hdr"); strcat(zstlist[1], extzst);
      strcpy(zstlist[2], ".img"); strcat(zstlist[2], extzst);
      zstp[0] = zstlist[0]; zstp[1] = zstlist[1]; zstp[2] = zstlist[2];
      ext = name + len - 8;
      strncpy(extcopy, ext, 9); extcopy[9] = '\0';
      make_lowercase(extcopy);
      if (compare_strlist(extcopy, zstp, 3) >= 0) {
         if (is_mixedcase(ext)) return NULL;
         return (char *)ext;
      }
   }
#endif

   return NULL;
}

static int nifti_fileexists(const char *fname)
{
   NIIFILE fp = nii_open(fname, "rb", nifti_is_gzfile(fname));
   if (fp) { nii_close(&fp); return 1; }
   return 0;
}

static int nifti_validfilename(const char *fname)
{
   if (!fname || !*fname) return 0;
   const char *ext = nifti_find_file_extension(fname);
   if (ext && ext == fname) return 0;
   return 1;
}

static char *nifti_makebasename(const char *fname)
{
   char *basename = nifti_strdup(fname);
   const char *ext = nifti_find_file_extension(basename);
   if (ext) basename[strlen(basename) - strlen(ext)] = '\0';
   return basename;
}

static char *nifti_findhdrname(const char *fname)
{
   char *basename, *hdrname;
   const char *ext;
   char elist[2][5] = { ".hdr", ".nii" };
   char extzip[4] = ".gz";
   char extzst[5] = ".zst";
   int efirst = 1;
   int eisupper = 0;

   if (!nifti_validfilename(fname)) return NULL;
   basename = nifti_makebasename(fname);
   if (!basename) return NULL;

   ext = nifti_find_file_extension(fname);
   if (ext) eisupper = is_uppercase(ext);

   if (ext && nifti_fileexists(fname)) {
      if (fileext_n_compare(ext, ".img", 4) != 0) {
         hdrname = nifti_strdup(fname);
         free(basename);
         return hdrname;
      } else
         efirst = 0;
   }

   if (eisupper) {
      make_uppercase(elist[0]); make_uppercase(elist[1]);
      make_uppercase(extzip); make_uppercase(extzst);
   }

   hdrname = (char *)calloc(1, strlen(basename) + 12);
   if (!hdrname) { free(basename); return NULL; }

   strcpy(hdrname, basename);
   strcat(hdrname, elist[efirst]);
#ifdef FSLSTYLE
   if (nifti_fileexists(hdrname)) {
      free(basename);
      char *gzname = (char *)calloc(1, strlen(hdrname) + 8);
      strcpy(gzname, hdrname);
      strcat(gzname, extzip);
      if (nifti_fileexists(gzname)) {
         fprintf(stderr, "Image Exception : Multiple possible filenames detected for basename (*.nii, *.nii.gz): %s\n", basename);
         free(gzname);
         exit(134);
      }
      free(gzname);
      return hdrname;
   }
#else
   if (nifti_fileexists(hdrname)) { free(basename); return hdrname; }
#endif
#ifdef HAVE_ZLIB
   strcat(hdrname, extzip);
   if (nifti_fileexists(hdrname)) { free(basename); return hdrname; }
#endif
#ifdef HAVE_ZSTD
   strcpy(hdrname, basename); strcat(hdrname, elist[efirst]); strcat(hdrname, extzst);
   if (nifti_fileexists(hdrname)) { free(basename); return hdrname; }
#endif

   efirst = 1 - efirst;
   strcpy(hdrname, basename);
   strcat(hdrname, elist[efirst]);
   if (nifti_fileexists(hdrname)) { free(basename); return hdrname; }
#ifdef HAVE_ZLIB
   strcat(hdrname, extzip);
   if (nifti_fileexists(hdrname)) { free(basename); return hdrname; }
#endif
#ifdef HAVE_ZSTD
   strcpy(hdrname, basename); strcat(hdrname, elist[efirst]); strcat(hdrname, extzst);
   if (nifti_fileexists(hdrname)) { free(basename); return hdrname; }
#endif

   free(basename);
   free(hdrname);
   return NULL;
}

static char *nifti_findimgname(const char *fname, int nifti_type)
{
   char *basename, *imgname;
   char elist[2][5] = { ".nii", ".img" };
   char extzip[4] = ".gz";
   char extzst[5] = ".zst";
   const char *ext;
   int first;

   if (!nifti_validfilename(fname)) return NULL;
   basename = nifti_makebasename(fname);
   imgname = (char *)calloc(1, strlen(basename) + 12);
   if (!imgname) { free(basename); return NULL; }

   ext = nifti_find_file_extension(fname);
   if (ext && is_uppercase(ext)) {
      make_uppercase(elist[0]); make_uppercase(elist[1]);
      make_uppercase(extzip); make_uppercase(extzst);
   }

   if (nifti_type == NIFTI_FTYPE_NIFTI1_1 || nifti_type == NIFTI_FTYPE_NIFTI2_1)
      first = 0;
   else
      first = 1;

   strcpy(imgname, basename);
   strcat(imgname, elist[first]);
   if (nifti_fileexists(imgname)) { free(basename); return imgname; }
#ifdef HAVE_ZLIB
   strcat(imgname, extzip);
   if (nifti_fileexists(imgname)) { free(basename); return imgname; }
#endif
#ifdef HAVE_ZSTD
   strcpy(imgname, basename); strcat(imgname, elist[first]); strcat(imgname, extzst);
   if (nifti_fileexists(imgname)) { free(basename); return imgname; }
#endif

   strcpy(imgname, basename);
   strcat(imgname, elist[1 - first]);
   if (nifti_fileexists(imgname)) { free(basename); return imgname; }
#ifdef HAVE_ZLIB
   strcat(imgname, extzip);
   if (nifti_fileexists(imgname)) { free(basename); return imgname; }
#endif
#ifdef HAVE_ZSTD
   strcpy(imgname, basename); strcat(imgname, elist[1 - first]); strcat(imgname, extzst);
   if (nifti_fileexists(imgname)) { free(basename); return imgname; }
#endif

   free(basename);
   free(imgname);
   return NULL;
}

/* comp: 0=none, 1=gzip(.gz), 2=zstd(.zst) */
static char *nifti_makehdrname(const char *prefix, int nifti_type, int check, int comp)
{
   char *iname;
   const char *ext;
   char extnii[5] = ".nii", exthdr[5] = ".hdr", extimg[5] = ".img";
   char extgz[5] = ".gz";
   char extzst[5] = ".zst";

   if (!nifti_validfilename(prefix)) return NULL;
   iname = (char *)calloc(1, strlen(prefix) + 12);
   if (!iname) return NULL;
   strcpy(iname, prefix);

   if ((ext = nifti_find_file_extension(iname)) != NULL) {
      if (is_uppercase(ext)) {
         make_uppercase(extnii); make_uppercase(exthdr);
         make_uppercase(extimg); make_uppercase(extgz); make_uppercase(extzst);
      }
      if (strncmp(ext, extimg, 4) == 0)
         memcpy(&(iname[strlen(iname) - strlen(ext)]), exthdr, 4);
   } else if (nifti_type == NIFTI_FTYPE_NIFTI1_1 || nifti_type == NIFTI_FTYPE_NIFTI2_1)
      strcat(iname, extnii);
   else
      strcat(iname, exthdr);

#ifdef HAVE_ZLIB
   if (comp == 1 && (!ext || !strstr(iname, extgz))) strcat(iname, extgz);
#endif
#ifdef HAVE_ZSTD
   if (comp == 2 && (!ext || !strstr(iname, extzst))) strcat(iname, extzst);
#endif
   if (comp == 0) { /* no-op */ }
   if (check && nifti_fileexists(iname)) { free(iname); return NULL; }
   return iname;
}

/* comp: 0=none, 1=gzip(.gz), 2=zstd(.zst) */
static char *nifti_makeimgname(const char *prefix, int nifti_type, int check, int comp)
{
   char *iname;
   const char *ext;
   char extnii[5] = ".nii", exthdr[5] = ".hdr", extimg[5] = ".img";
   char extgz[5] = ".gz";
   char extzst[5] = ".zst";

   if (!nifti_validfilename(prefix)) return NULL;
   iname = (char *)calloc(1, strlen(prefix) + 12);
   if (!iname) return NULL;
   strcpy(iname, prefix);

   if ((ext = nifti_find_file_extension(iname)) != NULL) {
      if (is_uppercase(ext)) {
         make_uppercase(extnii); make_uppercase(exthdr);
         make_uppercase(extimg); make_uppercase(extgz); make_uppercase(extzst);
      }
      if (strncmp(ext, exthdr, 4) == 0)
         memcpy(&(iname[strlen(iname) - strlen(ext)]), extimg, 4);
   } else if (nifti_type == NIFTI_FTYPE_NIFTI1_1 || nifti_type == NIFTI_FTYPE_NIFTI2_1)
      strcat(iname, extnii);
   else
      strcat(iname, extimg);

#ifdef HAVE_ZLIB
   if (comp == 1 && (!ext || !strstr(iname, extgz))) strcat(iname, extgz);
#endif
#ifdef HAVE_ZSTD
   if (comp == 2 && (!ext || !strstr(iname, extzst))) strcat(iname, extzst);
#endif
   if (comp == 0) { /* no-op */ }
   if (check && nifti_fileexists(iname)) { free(iname); return NULL; }
   return iname;
}

static int nifti_set_type_from_names(nifti_image *nim)
{
   if (!nim || !nim->fname || !nim->iname) return -1;
   if (!nifti_validfilename(nim->fname) || !nifti_validfilename(nim->iname) ||
       !nifti_find_file_extension(nim->fname) || !nifti_find_file_extension(nim->iname))
      return -1;
   if (strcmp(nim->fname, nim->iname) == 0) {
      nim->nifti_type = (nim->nifti_type >= NIFTI_FTYPE_NIFTI2_1) ?
                         NIFTI_FTYPE_NIFTI2_1 : NIFTI_FTYPE_NIFTI1_1;
   } else if (nim->nifti_type == NIFTI_FTYPE_NIFTI1_1)
      nim->nifti_type = NIFTI_FTYPE_NIFTI1_2;
   else if (nim->nifti_type == NIFTI_FTYPE_NIFTI2_1)
      nim->nifti_type = NIFTI_FTYPE_NIFTI2_2;
   return 0;
}

int nifti_set_filenames(nifti_image *nim, const char *prefix, int check, int set_byte_order)
{
   int comp = nifti_is_gzfile(prefix);
   if (!nim || !prefix) return -1;
   if (nim->fname) free(nim->fname);
   if (nim->iname) free(nim->iname);
   nim->iname = NULL;
   nim->fname = nifti_makehdrname(prefix, nim->nifti_type, check, comp);
   if (nim->fname)
      nim->iname = nifti_makeimgname(prefix, nim->nifti_type, check, comp);
   if (!nim->fname || !nim->iname) return -1;
   if (set_byte_order) nim->byteorder = nifti_short_order();
   if (nifti_set_type_from_names(nim) < 0) return -1;
   return 0;
}

/*========== Header version detection ==========*/

static int nifti_header_version(const char *buf, size_t nbytes)
{
   nifti_1_header *n1p = (nifti_1_header *)buf;
   nifti_2_header *n2p = (nifti_2_header *)buf;
   int sizeof_hdr, sver, nver;

   if (!buf || nbytes < sizeof(nifti_1_header)) return -1;

   sver = -1;
   sizeof_hdr = n1p->sizeof_hdr;
   if      (sizeof_hdr == (int)sizeof(nifti_1_header)) sver = 1;
   else if (sizeof_hdr == (int)sizeof(nifti_2_header)) sver = 2;
   else {
      nifti_swap_4bytes(1, &sizeof_hdr);
      if      (sizeof_hdr == (int)sizeof(nifti_1_header)) sver = 1;
      else if (sizeof_hdr == (int)sizeof(nifti_2_header)) sver = 2;
   }

   if (sver == 1) {
      nver = NIFTI_VERSION(*n1p);
      if (nver == 0) return 0;  /* ANALYZE */
      if (nver == 1) return 1;  /* NIFTI-1 */
      return -1;
   } else if (sver == 2) {
      nver = NIFTI_VERSION(*n2p);
      if (nver == 2) return 2;  /* NIFTI-2 */
      return -1;
   }
   return -1;
}

static int need_nhdr_swap(short dim0, int hdrsize)
{
   short d0 = dim0;
   int hsize = hdrsize;

   if (d0 != 0) {
      if (d0 > 0 && d0 <= 7) return 0;
      nifti_swap_2bytes(1, &d0);
      if (d0 > 0 && d0 <= 7) return 1;
      return -1;
   }
   if (hsize == (int)sizeof(nifti_1_header)) return 0;
   nifti_swap_4bytes(1, &hsize);
   if (hsize == (int)sizeof(nifti_1_header)) return 1;
   return -2;
}

/*========== Math - double precision ==========*/

nifti_dmat44 nifti_quatern_to_dmat44(double qb, double qc, double qd,
                                      double qx, double qy, double qz,
                                      double dx, double dy, double dz, double qfac)
{
   nifti_dmat44 R;
   double a, b = qb, c = qc, d = qd, xd, yd, zd;

   R.m[3][0] = R.m[3][1] = R.m[3][2] = 0.0; R.m[3][3] = 1.0;
   a = 1.0l - (b*b + c*c + d*d);
   if (a < 1.e-7l) {
      a = 1.0l / sqrt(b*b + c*c + d*d);
      b *= a; c *= a; d *= a; a = 0.0l;
   } else {
      a = sqrt(a);
   }
   xd = (dx > 0.0) ? dx : 1.0l;
   yd = (dy > 0.0) ? dy : 1.0l;
   zd = (dz > 0.0) ? dz : 1.0l;
   if (qfac < 0.0) zd = -zd;

   R.m[0][0] = (a*a+b*b-c*c-d*d) * xd;
   R.m[0][1] = 2.0l*(b*c-a*d) * yd;
   R.m[0][2] = 2.0l*(b*d+a*c) * zd;
   R.m[1][0] = 2.0l*(b*c+a*d) * xd;
   R.m[1][1] = (a*a+c*c-b*b-d*d) * yd;
   R.m[1][2] = 2.0l*(c*d-a*b) * zd;
   R.m[2][0] = 2.0l*(b*d-a*c) * xd;
   R.m[2][1] = 2.0l*(c*d+a*b) * yd;
   R.m[2][2] = (a*a+d*d-c*c-b*b) * zd;
   R.m[0][3] = qx; R.m[1][3] = qy; R.m[2][3] = qz;
   return R;
}

mat44 nifti_quatern_to_mat44(float qb, float qc, float qd,
                              float qx, float qy, float qz,
                              float dx, float dy, float dz, float qfac)
{
   mat44 R;
   double a, b = qb, c = qc, d = qd, xd, yd, zd;

   R.m[3][0] = R.m[3][1] = R.m[3][2] = 0.0f; R.m[3][3] = 1.0f;
   a = 1.0l - (b*b + c*c + d*d);
   if (a < 1.e-7l) {
      a = 1.0l / sqrt(b*b + c*c + d*d);
      b *= a; c *= a; d *= a; a = 0.0l;
   } else {
      a = sqrt(a);
   }
   xd = (dx > 0.0) ? dx : 1.0l;
   yd = (dy > 0.0) ? dy : 1.0l;
   zd = (dz > 0.0) ? dz : 1.0l;
   if (qfac < 0.0) zd = -zd;

   R.m[0][0] = (float)((a*a+b*b-c*c-d*d) * xd);
   R.m[0][1] = 2.0l*(b*c-a*d) * yd;
   R.m[0][2] = 2.0l*(b*d+a*c) * zd;
   R.m[1][0] = 2.0l*(b*c+a*d) * xd;
   R.m[1][1] = (float)((a*a+c*c-b*b-d*d) * yd);
   R.m[1][2] = 2.0l*(c*d-a*b) * zd;
   R.m[2][0] = 2.0l*(b*d-a*c) * xd;
   R.m[2][1] = 2.0l*(c*d+a*b) * yd;
   R.m[2][2] = (float)((a*a+d*d-c*c-b*b) * zd);
   R.m[0][3] = qx; R.m[1][3] = qy; R.m[2][3] = qz;
   return R;
}

/* Forward declarations for static helpers used by quatern functions */
static nifti_dmat33 nifti_dmat33_polar(nifti_dmat33 A);
static mat33 nifti_mat33_polar(mat33 A);

void nifti_dmat44_to_quatern(nifti_dmat44 R,
                              double *qb, double *qc, double *qd,
                              double *qx, double *qy, double *qz,
                              double *dx, double *dy, double *dz, double *qfac)
{
   double r11,r12,r13, r21,r22,r23, r31,r32,r33, xd,yd,zd, a,b,c,d;
   nifti_dmat33 P, Q;

   ASSIF(qx, R.m[0][3]); ASSIF(qy, R.m[1][3]); ASSIF(qz, R.m[2][3]);
   r11=R.m[0][0]; r12=R.m[0][1]; r13=R.m[0][2];
   r21=R.m[1][0]; r22=R.m[1][1]; r23=R.m[1][2];
   r31=R.m[2][0]; r32=R.m[2][1]; r33=R.m[2][2];

   xd = sqrt(r11*r11 + r21*r21 + r31*r31);
   yd = sqrt(r12*r12 + r22*r22 + r32*r32);
   zd = sqrt(r13*r13 + r23*r23 + r33*r33);

   if (xd == 0.0l) { r11 = 1.0l; r21 = r31 = 0.0l; xd = 1.0l; }
   if (yd == 0.0l) { r22 = 1.0l; r12 = r32 = 0.0l; yd = 1.0l; }
   if (zd == 0.0l) { r33 = 1.0l; r13 = r23 = 0.0l; zd = 1.0l; }

   ASSIF(dx, xd); ASSIF(dy, yd); ASSIF(dz, zd);

   r11 /= xd; r21 /= xd; r31 /= xd;
   r12 /= yd; r22 /= yd; r32 /= yd;
   r13 /= zd; r23 /= zd; r33 /= zd;

   Q.m[0][0]=r11; Q.m[0][1]=r12; Q.m[0][2]=r13;
   Q.m[1][0]=r21; Q.m[1][1]=r22; Q.m[1][2]=r23;
   Q.m[2][0]=r31; Q.m[2][1]=r32; Q.m[2][2]=r33;
   P = nifti_dmat33_polar(Q);
   r11=P.m[0][0]; r12=P.m[0][1]; r13=P.m[0][2];
   r21=P.m[1][0]; r22=P.m[1][1]; r23=P.m[1][2];
   r31=P.m[2][0]; r32=P.m[2][1]; r33=P.m[2][2];

   zd = r11*r22*r33-r11*r32*r23-r21*r12*r33+r21*r32*r13+r31*r12*r23-r31*r22*r13;
   if (zd > 0) { ASSIF(qfac, 1.0); }
   else         { ASSIF(qfac, -1.0); r13 = -r13; r23 = -r23; r33 = -r33; }

   a = r11 + r22 + r33 + 1.0l;
   if (a > 0.5l) {
      a = 0.5l * sqrt(a);
      b = 0.25l*(r32-r23)/a; c = 0.25l*(r13-r31)/a; d = 0.25l*(r21-r12)/a;
   } else {
      xd = 1.0+r11-(r22+r33); yd = 1.0+r22-(r11+r33); zd = 1.0+r33-(r11+r22);
      if (xd > 1.0) {
         b=0.5l*sqrt(xd); c=0.25l*(r12+r21)/b; d=0.25l*(r13+r31)/b; a=0.25l*(r32-r23)/b;
      } else if (yd > 1.0) {
         c=0.5l*sqrt(yd); b=0.25l*(r12+r21)/c; d=0.25l*(r23+r32)/c; a=0.25l*(r13-r31)/c;
      } else {
         d=0.5l*sqrt(zd); b=0.25l*(r13+r31)/d; c=0.25l*(r23+r32)/d; a=0.25l*(r21-r12)/d;
      }
      if (a < 0.0l) { b=-b; c=-c; d=-d; }
   }
   ASSIF(qb, b); ASSIF(qc, c); ASSIF(qd, d);
}

void nifti_mat44_to_quatern(mat44 R,
                             float *qb, float *qc, float *qd,
                             float *qx, float *qy, float *qz,
                             float *dx, float *dy, float *dz, float *qfac)
{
   double r11,r12,r13, r21,r22,r23, r31,r32,r33, xd,yd,zd, a,b,c,d;
   mat33 P, Q;

   ASSIF(qx, R.m[0][3]); ASSIF(qy, R.m[1][3]); ASSIF(qz, R.m[2][3]);
   r11=R.m[0][0]; r12=R.m[0][1]; r13=R.m[0][2];
   r21=R.m[1][0]; r22=R.m[1][1]; r23=R.m[1][2];
   r31=R.m[2][0]; r32=R.m[2][1]; r33=R.m[2][2];

   xd = sqrt(r11*r11 + r21*r21 + r31*r31);
   yd = sqrt(r12*r12 + r22*r22 + r32*r32);
   zd = sqrt(r13*r13 + r23*r23 + r33*r33);

   if (xd == 0.0l) { r11=1.0l; r21=r31=0.0l; xd=1.0l; }
   if (yd == 0.0l) { r22=1.0l; r12=r32=0.0l; yd=1.0l; }
   if (zd == 0.0l) { r33=1.0l; r13=r23=0.0l; zd=1.0l; }

   ASSIF(dx, (float)xd); ASSIF(dy, (float)yd); ASSIF(dz, (float)zd);
   r11/=xd; r21/=xd; r31/=xd; r12/=yd; r22/=yd; r32/=yd; r13/=zd; r23/=zd; r33/=zd;

   Q.m[0][0]=(float)r11; Q.m[0][1]=(float)r12; Q.m[0][2]=(float)r13;
   Q.m[1][0]=(float)r21; Q.m[1][1]=(float)r22; Q.m[1][2]=(float)r23;
   Q.m[2][0]=(float)r31; Q.m[2][1]=(float)r32; Q.m[2][2]=(float)r33;
   P = nifti_mat33_polar(Q);
   r11=P.m[0][0]; r12=P.m[0][1]; r13=P.m[0][2];
   r21=P.m[1][0]; r22=P.m[1][1]; r23=P.m[1][2];
   r31=P.m[2][0]; r32=P.m[2][1]; r33=P.m[2][2];

   zd = r11*r22*r33-r11*r32*r23-r21*r12*r33+r21*r32*r13+r31*r12*r23-r31*r22*r13;
   if (zd > 0) { ASSIF(qfac, 1.0f); }
   else         { ASSIF(qfac, -1.0f); r13=-r13; r23=-r23; r33=-r33; }

   a = r11+r22+r33+1.0l;
   if (a > 0.5l) {
      a=0.5l*sqrt(a); b=0.25l*(r32-r23)/a; c=0.25l*(r13-r31)/a; d=0.25l*(r21-r12)/a;
   } else {
      xd=1.0+r11-(r22+r33); yd=1.0+r22-(r11+r33); zd=1.0+r33-(r11+r22);
      if (xd > 1.0) { b=0.5l*sqrt(xd); c=0.25l*(r12+r21)/b; d=0.25l*(r13+r31)/b; a=0.25l*(r32-r23)/b; }
      else if (yd > 1.0) { c=0.5l*sqrt(yd); b=0.25l*(r12+r21)/c; d=0.25l*(r23+r32)/c; a=0.25l*(r13-r31)/c; }
      else { d=0.5l*sqrt(zd); b=0.25l*(r13+r31)/d; c=0.25l*(r23+r32)/d; a=0.25l*(r21-r12)/d; }
      if (a < 0.0l) { b=-b; c=-c; d=-d; }
   }
   ASSIF(qb, (float)b); ASSIF(qc, (float)c); ASSIF(qd, (float)d);
}

/*========== Matrix inverse, multiply, determinant ==========*/

#define MAT44_INV(TYPE, FNAME, CAST)                                   \
TYPE FNAME(TYPE R) {                                                    \
   double r11,r12,r13,r21,r22,r23,r31,r32,r33,v1,v2,v3,deti;         \
   TYPE Q;                                                              \
   r11=R.m[0][0]; r12=R.m[0][1]; r13=R.m[0][2];                      \
   r21=R.m[1][0]; r22=R.m[1][1]; r23=R.m[1][2];                      \
   r31=R.m[2][0]; r32=R.m[2][1]; r33=R.m[2][2];                      \
   v1=R.m[0][3]; v2=R.m[1][3]; v3=R.m[2][3];                         \
   deti = r11*r22*r33-r11*r32*r23-r21*r12*r33                         \
         +r21*r32*r13+r31*r12*r23-r31*r22*r13;                        \
   if (deti != 0.0l) deti = 1.0l / deti;                              \
   Q.m[0][0]=CAST(deti*(r22*r33-r32*r23));                            \
   Q.m[0][1]=CAST(deti*(-r12*r33+r32*r13));                           \
   Q.m[0][2]=CAST(deti*(r12*r23-r22*r13));                            \
   Q.m[0][3]=CAST(deti*(-r12*r23*v3+r12*v2*r33+r22*r13*v3            \
                        -r22*v1*r33-r32*r13*v2+r32*v1*r23));          \
   Q.m[1][0]=CAST(deti*(-r21*r33+r31*r23));                           \
   Q.m[1][1]=CAST(deti*(r11*r33-r31*r13));                            \
   Q.m[1][2]=CAST(deti*(-r11*r23+r21*r13));                           \
   Q.m[1][3]=CAST(deti*(r11*r23*v3-r11*v2*r33-r21*r13*v3             \
                        +r21*v1*r33+r31*r13*v2-r31*v1*r23));          \
   Q.m[2][0]=CAST(deti*(r21*r32-r31*r22));                            \
   Q.m[2][1]=CAST(deti*(-r11*r32+r31*r12));                           \
   Q.m[2][2]=CAST(deti*(r11*r22-r21*r12));                            \
   Q.m[2][3]=CAST(deti*(-r11*r22*v3+r11*r32*v2+r21*r12*v3            \
                        -r21*r32*v1-r31*r12*v2+r31*r22*v1));          \
   Q.m[3][0]=Q.m[3][1]=Q.m[3][2]=0.0l;                               \
   Q.m[3][3]=(deti==0.0l)?0.0l:1.0l;                                  \
   return Q;                                                            \
}

#define D_CAST(x) (x)
#define F_CAST(x) ((float)(x))

MAT44_INV(nifti_dmat44, nifti_dmat44_inverse, D_CAST)
MAT44_INV(mat44, nifti_mat44_inverse, F_CAST)

#define MAT33_INV(TYPE, FNAME, CAST)                                   \
static TYPE FNAME(TYPE R) {                                             \
   double r11,r12,r13,r21,r22,r23,r31,r32,r33,deti;                   \
   TYPE Q;                                                              \
   r11=R.m[0][0]; r12=R.m[0][1]; r13=R.m[0][2];                      \
   r21=R.m[1][0]; r22=R.m[1][1]; r23=R.m[1][2];                      \
   r31=R.m[2][0]; r32=R.m[2][1]; r33=R.m[2][2];                      \
   deti = r11*r22*r33-r11*r32*r23-r21*r12*r33                         \
         +r21*r32*r13+r31*r12*r23-r31*r22*r13;                        \
   if (deti != 0.0l) deti = 1.0l/deti;                                \
   Q.m[0][0]=CAST(deti*(r22*r33-r32*r23));                            \
   Q.m[0][1]=CAST(deti*(-r12*r33+r32*r13));                           \
   Q.m[0][2]=CAST(deti*(r12*r23-r22*r13));                            \
   Q.m[1][0]=CAST(deti*(-r21*r33+r31*r23));                           \
   Q.m[1][1]=CAST(deti*(r11*r33-r31*r13));                            \
   Q.m[1][2]=CAST(deti*(-r11*r23+r21*r13));                           \
   Q.m[2][0]=CAST(deti*(r21*r32-r31*r22));                            \
   Q.m[2][1]=CAST(deti*(-r11*r32+r31*r12));                           \
   Q.m[2][2]=CAST(deti*(r11*r22-r21*r12));                            \
   return Q;                                                            \
}

MAT33_INV(nifti_dmat33, nifti_dmat33_inverse, D_CAST)
MAT33_INV(mat33, nifti_mat33_inverse, F_CAST)

static double nifti_dmat33_determ(nifti_dmat33 R)
{
   double r11=R.m[0][0],r12=R.m[0][1],r13=R.m[0][2],
          r21=R.m[1][0],r22=R.m[1][1],r23=R.m[1][2],
          r31=R.m[2][0],r32=R.m[2][1],r33=R.m[2][2];
   return r11*r22*r33-r11*r32*r23-r21*r12*r33+r21*r32*r13+r31*r12*r23-r31*r22*r13;
}

float nifti_mat33_determ(mat33 R)
{
   double r11=R.m[0][0],r12=R.m[0][1],r13=R.m[0][2],
          r21=R.m[1][0],r22=R.m[1][1],r23=R.m[1][2],
          r31=R.m[2][0],r32=R.m[2][1],r33=R.m[2][2];
   return (float)(r11*r22*r33-r11*r32*r23-r21*r12*r33+r21*r32*r13+r31*r12*r23-r31*r22*r13);
}

static double nifti_dmat33_rownorm(nifti_dmat33 A)
{
   double r1=fabs(A.m[0][0])+fabs(A.m[0][1])+fabs(A.m[0][2]);
   double r2=fabs(A.m[1][0])+fabs(A.m[1][1])+fabs(A.m[1][2]);
   double r3=fabs(A.m[2][0])+fabs(A.m[2][1])+fabs(A.m[2][2]);
   if (r1 < r2) r1=r2; if (r1 < r3) r1=r3; return r1;
}

static float nifti_mat33_rownorm(mat33 A)
{
   float r1=(float)(fabs(A.m[0][0])+fabs(A.m[0][1])+fabs(A.m[0][2]));
   float r2=(float)(fabs(A.m[1][0])+fabs(A.m[1][1])+fabs(A.m[1][2]));
   float r3=(float)(fabs(A.m[2][0])+fabs(A.m[2][1])+fabs(A.m[2][2]));
   if (r1 < r2) r1=r2; if (r1 < r3) r1=r3; return r1;
}

static double nifti_dmat33_colnorm(nifti_dmat33 A)
{
   double r1=fabs(A.m[0][0])+fabs(A.m[1][0])+fabs(A.m[2][0]);
   double r2=fabs(A.m[0][1])+fabs(A.m[1][1])+fabs(A.m[2][1]);
   double r3=fabs(A.m[0][2])+fabs(A.m[1][2])+fabs(A.m[2][2]);
   if (r1 < r2) r1=r2; if (r1 < r3) r1=r3; return r1;
}

static float nifti_mat33_colnorm(mat33 A)
{
   float r1=(float)(fabs(A.m[0][0])+fabs(A.m[1][0])+fabs(A.m[2][0]));
   float r2=(float)(fabs(A.m[0][1])+fabs(A.m[1][1])+fabs(A.m[2][1]));
   float r3=(float)(fabs(A.m[0][2])+fabs(A.m[1][2])+fabs(A.m[2][2]));
   if (r1 < r2) r1=r2; if (r1 < r3) r1=r3; return r1;
}

mat44 nifti_mat44_mul(mat44 A, mat44 B)
{
   mat44 C;
   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) {
         C.m[i][j] = 0.0;
         for (int k = 0; k < 4; k++) C.m[i][j] += A.m[i][k]*B.m[k][j];
      }
   return C;
}

static nifti_dmat33 nifti_dmat33_polar(nifti_dmat33 A)
{
   nifti_dmat33 X, Y, Z;
   double alp, bet, gam, gmi, dif = 1.0;
   int k = 0;
   X = A;
   gam = nifti_dmat33_determ(X);
   while (gam == 0.0) {
      gam = 0.00001 * (0.001 + nifti_dmat33_rownorm(X));
      X.m[0][0] += gam; X.m[1][1] += gam; X.m[2][2] += gam;
      gam = nifti_dmat33_determ(X);
   }
   while (1) {
      Y = nifti_dmat33_inverse(X);
      if (dif > 0.3) {
         alp = sqrt(nifti_dmat33_rownorm(X)*nifti_dmat33_colnorm(X));
         bet = sqrt(nifti_dmat33_rownorm(Y)*nifti_dmat33_colnorm(Y));
         gam = sqrt(bet/alp); gmi = 1.0/gam;
      } else { gam = gmi = 1.0; }
      for (int i = 0; i < 3; i++)
         for (int j = 0; j < 3; j++)
            Z.m[i][j] = 0.5*(gam*X.m[i][j] + gmi*Y.m[j][i]);
      dif = 0;
      for (int i = 0; i < 3; i++)
         for (int j = 0; j < 3; j++)
            dif += fabs(Z.m[i][j] - X.m[i][j]);
      k++;
      if (k > 100 || dif < 3.e-6) break;
      X = Z;
   }
   return Z;
}

static mat33 nifti_mat33_polar(mat33 A)
{
   mat33 X, Y, Z;
   float alp, bet, gam, gmi, dif = 1.0f;
   int k = 0;
   X = A;
   gam = nifti_mat33_determ(X);
   while (gam == 0.0) {
      gam = (float)(0.00001 * (0.001 + nifti_mat33_rownorm(X)));
      X.m[0][0] += gam; X.m[1][1] += gam; X.m[2][2] += gam;
      gam = nifti_mat33_determ(X);
   }
   while (1) {
      Y = nifti_mat33_inverse(X);
      if (dif > 0.3) {
         alp = (float)sqrt(nifti_mat33_rownorm(X)*nifti_mat33_colnorm(X));
         bet = (float)sqrt(nifti_mat33_rownorm(Y)*nifti_mat33_colnorm(Y));
         gam = (float)sqrt(bet/alp); gmi = (float)(1.0/gam);
      } else { gam = gmi = 1.0f; }
      for (int i = 0; i < 3; i++)
         for (int j = 0; j < 3; j++)
            Z.m[i][j] = (float)(0.5*(gam*X.m[i][j] + gmi*Y.m[j][i]));
      dif = 0;
      for (int i = 0; i < 3; i++)
         for (int j = 0; j < 3; j++)
            dif += (float)fabs(Z.m[i][j] - X.m[i][j]);
      k++;
      if (k > 100 || dif < 3.e-6) break;
      X = Z;
   }
   return Z;
}

/*========== Header ↔ nifti_image conversion ==========*/

static nifti_image *nifti_convert_n1hdr2nim(nifti_1_header nhdr, const char *fname)
{
   int ii, doswap, ioff, ni_ver, is_onefile;
   nifti_image *nim = (nifti_image *)calloc(1, sizeof(nifti_image));
   if (!nim) { fprintf(stderr, "** n1hdr2nim: alloc failed\n"); return NULL; }

   doswap = need_nhdr_swap(nhdr.dim[0], nhdr.sizeof_hdr);
   if (doswap < 0) { free(nim); return NULL; }

   ni_ver = NIFTI_VERSION(nhdr);
   if (ni_ver == 0) {
      unsigned char c = *((char *)(&nhdr.qform_code));
      nim->analyze75_orient = (analyze_75_orient_code)c;
   }
   if (doswap) swap_nifti_header(&nhdr, ni_ver);

   if (nhdr.datatype == DT_BINARY || nhdr.datatype == DT_UNKNOWN) { free(nim); return NULL; }
   if (nhdr.dim[1] <= 0) { free(nim); return NULL; }
   for (ii = 2; ii <= nhdr.dim[0]; ii++) if (nhdr.dim[ii] <= 0) nhdr.dim[ii] = 1;
   for (ii = nhdr.dim[0]+1; ii <= 7; ii++) if (nhdr.dim[ii] != 1 && nhdr.dim[ii] != 0) nhdr.dim[ii] = 1;
   for (ii = 1; ii <= nhdr.dim[0]; ii++)
      if (nhdr.pixdim[ii] == 0.0 || !IS_GOOD_FLOAT(nhdr.pixdim[ii])) nhdr.pixdim[ii] = 1.0f;

   is_onefile = (ni_ver > 0) && NIFTI_ONEFILE(nhdr);
   if (ni_ver) nim->nifti_type = is_onefile ? NIFTI_FTYPE_NIFTI1_1 : NIFTI_FTYPE_NIFTI1_2;
   else        nim->nifti_type = NIFTI_FTYPE_ANALYZE;

   ii = nifti_short_order();
   nim->byteorder = doswap ? REVERSE_ORDER(ii) : ii;

   if (nhdr.dim[0] < 1 || nhdr.dim[0] > 7) { free(nim); return NULL; }
   nim->ndim=nim->dim[0]=nhdr.dim[0];
   nim->nx=nim->dim[1]=nhdr.dim[1]; nim->ny=nim->dim[2]=nhdr.dim[2];
   nim->nz=nim->dim[3]=nhdr.dim[3]; nim->nt=nim->dim[4]=nhdr.dim[4];
   nim->nu=nim->dim[5]=nhdr.dim[5]; nim->nv=nim->dim[6]=nhdr.dim[6];
   nim->nw=nim->dim[7]=nhdr.dim[7];
   for (ii=1, nim->nvox=1; ii <= nhdr.dim[0]; ii++) {
      if (nhdr.dim[ii] < 1) { free(nim); return NULL; }
      if (nim->nvox > INT64_MAX / nhdr.dim[ii]) { free(nim); return NULL; }
      nim->nvox *= nhdr.dim[ii];
   }

   nim->datatype = nhdr.datatype;
   nifti_datatype_sizes(nim->datatype, &nim->nbyper, &nim->swapsize);
   if (nim->nbyper == 0) { free(nim); return NULL; }

   nim->dx=nim->pixdim[1]=nhdr.pixdim[1]; nim->dy=nim->pixdim[2]=nhdr.pixdim[2];
   nim->dz=nim->pixdim[3]=nhdr.pixdim[3]; nim->dt=nim->pixdim[4]=nhdr.pixdim[4];
   nim->du=nim->pixdim[5]=nhdr.pixdim[5]; nim->dv=nim->pixdim[6]=nhdr.pixdim[6];
   nim->dw=nim->pixdim[7]=nhdr.pixdim[7];

   if (!ni_ver || nhdr.qform_code <= 0) {
      nim->qto_xyz.m[0][0]=nim->dx; nim->qto_xyz.m[1][1]=nim->dy; nim->qto_xyz.m[2][2]=nim->dz;
      nim->qto_xyz.m[0][1]=nim->qto_xyz.m[0][2]=nim->qto_xyz.m[0][3]=0.0;
      nim->qto_xyz.m[1][0]=nim->qto_xyz.m[1][2]=nim->qto_xyz.m[1][3]=0.0;
      nim->qto_xyz.m[2][0]=nim->qto_xyz.m[2][1]=nim->qto_xyz.m[2][3]=0.0;
      nim->qto_xyz.m[3][0]=nim->qto_xyz.m[3][1]=nim->qto_xyz.m[3][2]=0.0; nim->qto_xyz.m[3][3]=1.0;
      nim->qform_code = NIFTI_XFORM_UNKNOWN;
   } else {
      nim->quatern_b=FIXED_FLOAT(nhdr.quatern_b); nim->quatern_c=FIXED_FLOAT(nhdr.quatern_c);
      nim->quatern_d=FIXED_FLOAT(nhdr.quatern_d);
      nim->qoffset_x=FIXED_FLOAT(nhdr.qoffset_x); nim->qoffset_y=FIXED_FLOAT(nhdr.qoffset_y);
      nim->qoffset_z=FIXED_FLOAT(nhdr.qoffset_z);
      nim->qfac = (nhdr.pixdim[0] < 0.0) ? -1.0 : 1.0;
      nim->qto_xyz = nifti_quatern_to_dmat44(nim->quatern_b, nim->quatern_c, nim->quatern_d,
                        nim->qoffset_x, nim->qoffset_y, nim->qoffset_z,
                        nim->dx, nim->dy, nim->dz, nim->qfac);
      nim->qform_code = nhdr.qform_code;
   }
   nim->qto_ijk = nifti_dmat44_inverse(nim->qto_xyz);

   if (!ni_ver || nhdr.sform_code <= 0) {
      nim->sform_code = NIFTI_XFORM_UNKNOWN;
   } else {
      for (int j = 0; j < 4; j++) {
         nim->sto_xyz.m[0][j] = nhdr.srow_x[j];
         nim->sto_xyz.m[1][j] = nhdr.srow_y[j];
         nim->sto_xyz.m[2][j] = nhdr.srow_z[j];
      }
      nim->sto_xyz.m[3][0]=nim->sto_xyz.m[3][1]=nim->sto_xyz.m[3][2]=0.0; nim->sto_xyz.m[3][3]=1.0;
      nim->sto_ijk = nifti_dmat44_inverse(nim->sto_xyz);
      nim->sform_code = nhdr.sform_code;
   }

   if (ni_ver) {
      nim->scl_slope=FIXED_FLOAT(nhdr.scl_slope); nim->scl_inter=FIXED_FLOAT(nhdr.scl_inter);
      nim->intent_code=nhdr.intent_code;
      nim->intent_p1=FIXED_FLOAT(nhdr.intent_p1); nim->intent_p2=FIXED_FLOAT(nhdr.intent_p2);
      nim->intent_p3=FIXED_FLOAT(nhdr.intent_p3);
      nim->toffset=FIXED_FLOAT(nhdr.toffset);
      memcpy(nim->intent_name, nhdr.intent_name, 15); nim->intent_name[15]='\0';
      nim->xyz_units=XYZT_TO_SPACE(nhdr.xyzt_units); nim->time_units=XYZT_TO_TIME(nhdr.xyzt_units);
      nim->freq_dim=DIM_INFO_TO_FREQ_DIM(nhdr.dim_info);
      nim->phase_dim=DIM_INFO_TO_PHASE_DIM(nhdr.dim_info);
      nim->slice_dim=DIM_INFO_TO_SLICE_DIM(nhdr.dim_info);
      nim->slice_code=nhdr.slice_code; nim->slice_start=nhdr.slice_start;
      nim->slice_end=nhdr.slice_end; nim->slice_duration=FIXED_FLOAT(nhdr.slice_duration);
   }
   nim->cal_min=FIXED_FLOAT(nhdr.cal_min); nim->cal_max=FIXED_FLOAT(nhdr.cal_max);
   memcpy(nim->descrip, nhdr.descrip, 79); nim->descrip[79]='\0';
   memcpy(nim->aux_file, nhdr.aux_file, 23); nim->aux_file[23]='\0';

   is_onefile = ni_ver && NIFTI_ONEFILE(nhdr);
   if (is_onefile) {
      ioff = (int)nhdr.vox_offset;
      if (ioff < (int)sizeof(nhdr)) ioff = (int)sizeof(nhdr);
   } else {
      ioff = (int)nhdr.vox_offset;
   }
   nim->iname_offset = ioff;

   if (fname) {
      nifti_set_filenames(nim, fname, 0, 0);
      if (!nim->iname) { free(nim); return NULL; }
   }
   nim->num_ext = 0; nim->ext_list = NULL;
   return nim;
}

static nifti_image *nifti_convert_n2hdr2nim(nifti_2_header nhdr, const char *fname)
{
   int ii, doswap, ni_ver, is_onefile;
   nifti_image *nim = (nifti_image *)calloc(1, sizeof(nifti_image));
   if (!nim) { fprintf(stderr, "** n2hdr2nim: alloc failed\n"); return NULL; }

   doswap = NIFTI2_NEEDS_SWAP(nhdr);
   ni_ver = NIFTI_VERSION(nhdr);
   if (ni_ver != 2) { free(nim); return NULL; }
   if (doswap) swap_nifti_header(&nhdr, 2);

   if (nhdr.datatype == DT_BINARY || nhdr.datatype == DT_UNKNOWN) { free(nim); return NULL; }
   if (nhdr.dim[1] <= 0) { free(nim); return NULL; }
   for (ii = 2; ii <= nhdr.dim[0]; ii++) if (nhdr.dim[ii] <= 0) nhdr.dim[ii] = 1;
   for (ii = nhdr.dim[0]+1; ii <= 7; ii++) if (nhdr.dim[ii] != 1 && nhdr.dim[ii] != 0) nhdr.dim[ii] = 1;
   for (ii = 1; ii <= nhdr.dim[0]; ii++)
      if (nhdr.pixdim[ii] == 0.0 || !IS_GOOD_FLOAT(nhdr.pixdim[ii])) nhdr.pixdim[ii] = 1.0;

   is_onefile = NIFTI_ONEFILE(nhdr);
   nim->nifti_type = is_onefile ? NIFTI_FTYPE_NIFTI2_1 : NIFTI_FTYPE_NIFTI2_2;
   ii = nifti_short_order();
   nim->byteorder = doswap ? REVERSE_ORDER(ii) : ii;

   if (nhdr.dim[0] < 1 || nhdr.dim[0] > 7) { free(nim); return NULL; }
   nim->ndim=nim->dim[0]=nhdr.dim[0];
   nim->nx=nim->dim[1]=nhdr.dim[1]; nim->ny=nim->dim[2]=nhdr.dim[2];
   nim->nz=nim->dim[3]=nhdr.dim[3]; nim->nt=nim->dim[4]=nhdr.dim[4];
   nim->nu=nim->dim[5]=nhdr.dim[5]; nim->nv=nim->dim[6]=nhdr.dim[6];
   nim->nw=nim->dim[7]=nhdr.dim[7];
   for (ii=1, nim->nvox=1; ii <= nhdr.dim[0]; ii++) {
      if (nhdr.dim[ii] < 1) { free(nim); return NULL; }
      if (nim->nvox > INT64_MAX / nhdr.dim[ii]) { free(nim); return NULL; }
      nim->nvox *= nhdr.dim[ii];
   }

   nim->datatype = nhdr.datatype;
   nifti_datatype_sizes(nim->datatype, &nim->nbyper, &nim->swapsize);
   if (nim->nbyper == 0) { free(nim); return NULL; }

   nim->dx=nim->pixdim[1]=nhdr.pixdim[1]; nim->dy=nim->pixdim[2]=nhdr.pixdim[2];
   nim->dz=nim->pixdim[3]=nhdr.pixdim[3]; nim->dt=nim->pixdim[4]=nhdr.pixdim[4];
   nim->du=nim->pixdim[5]=nhdr.pixdim[5]; nim->dv=nim->pixdim[6]=nhdr.pixdim[6];
   nim->dw=nim->pixdim[7]=nhdr.pixdim[7];

   if (nhdr.qform_code <= 0) {
      nim->qto_xyz.m[0][0]=nim->dx; nim->qto_xyz.m[1][1]=nim->dy; nim->qto_xyz.m[2][2]=nim->dz;
      nim->qto_xyz.m[0][1]=nim->qto_xyz.m[0][2]=nim->qto_xyz.m[0][3]=0.0;
      nim->qto_xyz.m[1][0]=nim->qto_xyz.m[1][2]=nim->qto_xyz.m[1][3]=0.0;
      nim->qto_xyz.m[2][0]=nim->qto_xyz.m[2][1]=nim->qto_xyz.m[2][3]=0.0;
      nim->qto_xyz.m[3][0]=nim->qto_xyz.m[3][1]=nim->qto_xyz.m[3][2]=0.0; nim->qto_xyz.m[3][3]=1.0;
      nim->qform_code = NIFTI_XFORM_UNKNOWN;
   } else {
      nim->quatern_b=FIXED_FLOAT(nhdr.quatern_b); nim->quatern_c=FIXED_FLOAT(nhdr.quatern_c);
      nim->quatern_d=FIXED_FLOAT(nhdr.quatern_d);
      nim->qoffset_x=FIXED_FLOAT(nhdr.qoffset_x); nim->qoffset_y=FIXED_FLOAT(nhdr.qoffset_y);
      nim->qoffset_z=FIXED_FLOAT(nhdr.qoffset_z);
      nim->qfac = (nhdr.pixdim[0] < 0.0) ? -1.0 : 1.0;
      nim->qto_xyz = nifti_quatern_to_dmat44(nim->quatern_b, nim->quatern_c, nim->quatern_d,
                        nim->qoffset_x, nim->qoffset_y, nim->qoffset_z,
                        nim->dx, nim->dy, nim->dz, nim->qfac);
      nim->qform_code = nhdr.qform_code;
   }
   nim->qto_ijk = nifti_dmat44_inverse(nim->qto_xyz);

   if (nhdr.sform_code <= 0) {
      nim->sform_code = NIFTI_XFORM_UNKNOWN;
   } else {
      for (int j = 0; j < 4; j++) {
         nim->sto_xyz.m[0][j] = nhdr.srow_x[j];
         nim->sto_xyz.m[1][j] = nhdr.srow_y[j];
         nim->sto_xyz.m[2][j] = nhdr.srow_z[j];
      }
      nim->sto_xyz.m[3][0]=nim->sto_xyz.m[3][1]=nim->sto_xyz.m[3][2]=0.0; nim->sto_xyz.m[3][3]=1.0;
      nim->sto_ijk = nifti_dmat44_inverse(nim->sto_xyz);
      nim->sform_code = nhdr.sform_code;
   }

   nim->scl_slope=FIXED_FLOAT(nhdr.scl_slope); nim->scl_inter=FIXED_FLOAT(nhdr.scl_inter);
   nim->intent_code=nhdr.intent_code;
   nim->intent_p1=FIXED_FLOAT(nhdr.intent_p1); nim->intent_p2=FIXED_FLOAT(nhdr.intent_p2);
   nim->intent_p3=FIXED_FLOAT(nhdr.intent_p3);
   nim->toffset=FIXED_FLOAT(nhdr.toffset);
   memcpy(nim->intent_name, nhdr.intent_name, 15); nim->intent_name[15]='\0';
   nim->xyz_units=XYZT_TO_SPACE(nhdr.xyzt_units); nim->time_units=XYZT_TO_TIME(nhdr.xyzt_units);
   nim->freq_dim=DIM_INFO_TO_FREQ_DIM(nhdr.dim_info);
   nim->phase_dim=DIM_INFO_TO_PHASE_DIM(nhdr.dim_info);
   nim->slice_dim=DIM_INFO_TO_SLICE_DIM(nhdr.dim_info);
   nim->slice_code=nhdr.slice_code; nim->slice_start=nhdr.slice_start;
   nim->slice_end=nhdr.slice_end; nim->slice_duration=FIXED_FLOAT(nhdr.slice_duration);
   nim->cal_min=FIXED_FLOAT(nhdr.cal_min); nim->cal_max=FIXED_FLOAT(nhdr.cal_max);
   memcpy(nim->descrip, nhdr.descrip, 79); nim->descrip[79]='\0';
   memcpy(nim->aux_file, nhdr.aux_file, 23); nim->aux_file[23]='\0';

   nim->iname_offset = nhdr.vox_offset;
   if (is_onefile && nhdr.vox_offset < (int64_t)sizeof(nhdr))
      nim->iname_offset = (int64_t)sizeof(nhdr);

   if (fname) {
      nifti_set_filenames(nim, fname, 0, 0);
      if (!nim->iname) { free(nim); return NULL; }
   }
   nim->num_ext = 0; nim->ext_list = NULL;
   return nim;
}

/*========== nim → header conversion ==========*/

#define N_CHECK_2BYTE_VAL(fn) do { if (!NIFTI_IS_16_BIT_INT(nim->fn)) { \
   fprintf(stderr,"** nim->%s = %" PRId64 " does not fit in NIfTI-1\n", #fn, (int64_t)nim->fn); return 1; } } while(0)

static int nifti_convert_nim2n1hdr(const nifti_image *nim, nifti_1_header *hdr)
{
   nifti_1_header nhdr;
   if (!hdr) return 1;
   memset(&nhdr, 0, sizeof(nhdr));
   nhdr.sizeof_hdr = sizeof(nhdr);
   nhdr.regular = 'r';
   N_CHECK_2BYTE_VAL(ndim); N_CHECK_2BYTE_VAL(nx); N_CHECK_2BYTE_VAL(ny);
   N_CHECK_2BYTE_VAL(nz); N_CHECK_2BYTE_VAL(nt); N_CHECK_2BYTE_VAL(nu);
   N_CHECK_2BYTE_VAL(nv); N_CHECK_2BYTE_VAL(nw);
   N_CHECK_2BYTE_VAL(datatype); N_CHECK_2BYTE_VAL(nbyper);
   nhdr.dim[0]=nim->ndim; nhdr.dim[1]=nim->nx; nhdr.dim[2]=nim->ny; nhdr.dim[3]=nim->nz;
   nhdr.dim[4]=nim->nt; nhdr.dim[5]=nim->nu; nhdr.dim[6]=nim->nv; nhdr.dim[7]=nim->nw;
   nhdr.pixdim[0]=0.0f; nhdr.pixdim[1]=nim->dx; nhdr.pixdim[2]=nim->dy;
   nhdr.pixdim[3]=nim->dz; nhdr.pixdim[4]=nim->dt; nhdr.pixdim[5]=nim->du;
   nhdr.pixdim[6]=nim->dv; nhdr.pixdim[7]=nim->dw;
   nhdr.datatype=nim->datatype; nhdr.bitpix=8*nim->nbyper;
   if (nim->cal_max > nim->cal_min) { nhdr.cal_max=nim->cal_max; nhdr.cal_min=nim->cal_min; }
   if (nim->scl_slope != 0.0) { nhdr.scl_slope=nim->scl_slope; nhdr.scl_inter=nim->scl_inter; }
   if (nim->descrip[0]) { memcpy(nhdr.descrip, nim->descrip, 79); nhdr.descrip[79]='\0'; }
   if (nim->aux_file[0]) { memcpy(nhdr.aux_file, nim->aux_file, 23); nhdr.aux_file[23]='\0'; }

   if (nim->nifti_type > NIFTI_FTYPE_ANALYZE) {
      if (nim->nifti_type == NIFTI_FTYPE_NIFTI1_1) strcpy(nhdr.magic, "n+1");
      else strcpy(nhdr.magic, "ni1");
      for (int i = 1; i <= 7; i++) nhdr.pixdim[i] = (float)fabs(nhdr.pixdim[i]);
      N_CHECK_2BYTE_VAL(intent_code); N_CHECK_2BYTE_VAL(qform_code); N_CHECK_2BYTE_VAL(sform_code);
      nhdr.intent_code=nim->intent_code; nhdr.intent_p1=nim->intent_p1;
      nhdr.intent_p2=nim->intent_p2; nhdr.intent_p3=nim->intent_p3;
      if (nim->intent_name[0]) { memcpy(nhdr.intent_name, nim->intent_name, 15); nhdr.intent_name[15]='\0'; }
      nhdr.vox_offset=(float)nim->iname_offset;
      nhdr.xyzt_units=SPACE_TIME_TO_XYZT(nim->xyz_units, nim->time_units);
      nhdr.toffset=nim->toffset;
      if (nim->qform_code > 0) {
         nhdr.qform_code=nim->qform_code; nhdr.quatern_b=nim->quatern_b;
         nhdr.quatern_c=nim->quatern_c; nhdr.quatern_d=nim->quatern_d;
         nhdr.qoffset_x=nim->qoffset_x; nhdr.qoffset_y=nim->qoffset_y; nhdr.qoffset_z=nim->qoffset_z;
         nhdr.pixdim[0] = (nim->qfac >= 0.0) ? 1.0f : -1.0f;
      }
#ifdef FSLSTYLE
      else nhdr.pixdim[0] = 1.0;
#endif
      if (nim->sform_code > 0) {
         nhdr.sform_code=nim->sform_code;
         for (int j = 0; j < 4; j++) {
            nhdr.srow_x[j]=nim->sto_xyz.m[0][j]; nhdr.srow_y[j]=nim->sto_xyz.m[1][j]; nhdr.srow_z[j]=nim->sto_xyz.m[2][j];
         }
      }
      N_CHECK_2BYTE_VAL(slice_start); N_CHECK_2BYTE_VAL(slice_end);
      nhdr.dim_info=FPS_INTO_DIM_INFO(nim->freq_dim, nim->phase_dim, nim->slice_dim);
      nhdr.slice_code=nim->slice_code; nhdr.slice_start=nim->slice_start;
      nhdr.slice_end=nim->slice_end; nhdr.slice_duration=nim->slice_duration;
   }
   memcpy(hdr, &nhdr, sizeof(nhdr));
   return 0;
}

static int nifti_convert_nim2n2hdr(const nifti_image *nim, nifti_2_header *hdr)
{
   nifti_2_header nhdr;
   if (!hdr) return 1;
   memset(&nhdr, 0, sizeof(nhdr));
   nhdr.sizeof_hdr = sizeof(nhdr);
   memcpy(nhdr.magic, nifti2_magic, 8);
   if (nim->nifti_type == NIFTI_FTYPE_NIFTI2_2) nhdr.magic[1] = 'i';
   nhdr.datatype=nim->datatype; nhdr.bitpix=8*nim->nbyper;
   nhdr.dim[0]=nim->ndim; nhdr.dim[1]=nim->nx; nhdr.dim[2]=nim->ny; nhdr.dim[3]=nim->nz;
   nhdr.dim[4]=nim->nt; nhdr.dim[5]=nim->nu; nhdr.dim[6]=nim->nv; nhdr.dim[7]=nim->nw;
   nhdr.intent_p1=nim->intent_p1; nhdr.intent_p2=nim->intent_p2; nhdr.intent_p3=nim->intent_p3;
   nhdr.pixdim[0]=0.0;
   for (int i = 1; i <= 7; i++) nhdr.pixdim[i] = fabs(nim->pixdim[i]);
   nhdr.vox_offset=nim->iname_offset;
   nhdr.scl_slope=nim->scl_slope; nhdr.scl_inter=nim->scl_inter;
   nhdr.cal_max=nim->cal_max; nhdr.cal_min=nim->cal_min;
   nhdr.slice_duration=nim->slice_duration; nhdr.toffset=nim->toffset;
   nhdr.slice_start=nim->slice_start; nhdr.slice_end=nim->slice_end;
   if (nim->descrip[0]) { memcpy(nhdr.descrip, nim->descrip, 79); nhdr.descrip[79]='\0'; }
   if (nim->aux_file[0]) { memcpy(nhdr.aux_file, nim->aux_file, 23); nhdr.aux_file[23]='\0'; }
   if (nim->qform_code > 0) {
      nhdr.qform_code=nim->qform_code; nhdr.quatern_b=nim->quatern_b;
      nhdr.quatern_c=nim->quatern_c; nhdr.quatern_d=nim->quatern_d;
      nhdr.qoffset_x=nim->qoffset_x; nhdr.qoffset_y=nim->qoffset_y; nhdr.qoffset_z=nim->qoffset_z;
      nhdr.pixdim[0] = (nim->qfac >= 0.0) ? 1.0 : -1.0;
   }
#ifdef FSLSTYLE
   else nhdr.pixdim[0] = 1.0;
#endif
   if (nim->sform_code > 0) {
      nhdr.sform_code=nim->sform_code;
      for (int j = 0; j < 4; j++) {
         nhdr.srow_x[j]=nim->sto_xyz.m[0][j]; nhdr.srow_y[j]=nim->sto_xyz.m[1][j]; nhdr.srow_z[j]=nim->sto_xyz.m[2][j];
      }
   }
   nhdr.slice_code=nim->slice_code;
   nhdr.xyzt_units=SPACE_TIME_TO_XYZT(nim->xyz_units, nim->time_units);
   nhdr.intent_code=nim->intent_code;
   if (nim->intent_name[0]) { memcpy(nhdr.intent_name, nim->intent_name, 15); nhdr.intent_name[15]='\0'; }
   nhdr.dim_info=FPS_INTO_DIM_INFO(nim->freq_dim, nim->phase_dim, nim->slice_dim);
   memcpy(hdr, &nhdr, sizeof(nhdr));
   return 0;
}

/*========== Extensions ==========*/

static int nifti_is_valid_ecode(int ecode)
{
   return (ecode >= NIFTI_ECODE_IGNORE && ecode <= NIFTI_MAX_ECODE && !(ecode & 1));
}

static int nifti_check_extension(int size, int code, int64_t rem)
{
   (void)code; /* invalid ecode is not fatal */
   if (size < 16 || size > rem || (size & 0xf)) return 0;
   return 1;
}

static int nifti_read_extensions(nifti_image *nim, NIIFILE fp, int64_t remain)
{
   char extdr[4];
   nifti1_extension extn, *Elist = NULL;
   int64_t count;
   int swap;

   if (!nim || !fp) return -1;
   if (remain < 16) return 0;

   if (nii_read(extdr, 1, 4, fp) < 4) return 0;
   if (extdr[0] != 1) return 0;
   remain -= 4;

   swap = (nim->byteorder != nifti_short_order());
   count = 0;
   while (remain >= 16) {
      int size, code = -1;
      if (nii_read(&size, 4, 1, fp) != 1) break;
      if (nii_read(&code, 4, 1, fp) != 1) break;
      if (swap) { nifti_swap_4bytes(1, &size); nifti_swap_4bytes(1, &code); }
      if (!nifti_check_extension(size, code, remain)) break;

      extn.esize = size;
      extn.ecode = code;
      size -= 8;
      extn.edata = (char *)malloc(size);
      if (!extn.edata) break;
      if ((int)nii_read(extn.edata, 1, size, fp) < size) { free(extn.edata); break; }

      nifti1_extension *tmp = Elist;
      Elist = (nifti1_extension *)malloc((count + 1) * sizeof(nifti1_extension));
      if (!Elist) { free(extn.edata); Elist = tmp; break; }
      if (tmp) { memcpy(Elist, tmp, count * sizeof(nifti1_extension)); free(tmp); }
      Elist[count] = extn;
      remain -= extn.esize;
      count++;
   }
   nim->num_ext = (int)count;
   nim->ext_list = Elist;
   return (int)count;
}

static int nifti_extension_size(nifti_image *nim)
{
   int size = 0;
   if (!nim || nim->num_ext <= 0) return 0;
   for (int c = 0; c < nim->num_ext; c++) size += nim->ext_list[c].esize;
   return size;
}

static int nifti_write_extensions(NIIFILE fp, nifti_image *nim)
{
   char extdr[4] = { 0, 0, 0, 0 };
   if (!fp || !nim) return -1;
   if (nim->num_ext > 0) extdr[0] = 1;
   if (nii_write(extdr, 1, 4, fp) != 4) return -1;
   for (int c = 0; c < nim->num_ext; c++) {
      nifti1_extension *e = &nim->ext_list[c];
      if (e->esize < 16 || !e->edata) continue;
      if (nii_write(&e->esize, 4, 1, fp) != 1) return -1;
      if (nii_write(&e->ecode, 4, 1, fp) != 1) return -1;
      if (nii_write(e->edata, 1, e->esize - 8, fp) != (size_t)(e->esize - 8)) return -1;
   }
   return nim->num_ext;
}

/*========== Image read ==========*/

nifti_image *nifti_image_read(const char *hname, int read_data)
{
   nifti_1_header n1hdr;
   nifti_2_header n2hdr;
   nifti_image *nim;
   NIIFILE fp;
   int ii, ni_ver, onefile = 0;
   int64_t filesize, remain;
   int64_t h1size = sizeof(nifti_1_header), h2size = sizeof(nifti_2_header);
   char *hfile = NULL, *posn;

   int isStdIn = (hname && hname[0] == '-' && hname[1] == '\0');
   if (!isStdIn) {
      hfile = nifti_findhdrname(hname);
      if (!hfile) {
#ifndef HAVE_ZLIB
         nifti_is_gzfile(hname);
#endif
         fprintf(stderr, "** failed to find header file for '%s'\n", hname);
         return NULL;
      }
      filesize = nifti_is_gzfile(hfile) ? -1 : nifti_get_filesize(hfile);
   }

   if (isStdIn) {
      filesize = -1;
      if (!stdin_has_data()) { fprintf(stderr, " unable to read piped input buffer\n"); return NULL; }
      fp = nii_open("-", "rb", 0);
   } else {
      fp = nii_open(hfile, "rb", nifti_is_gzfile(hfile));
   }
   if (!fp) { free(hfile); return NULL; }

   ii = (int)nii_read(&n1hdr, 1, h1size, fp);
   if (ii < (int)h1size) { nii_close(&fp); free(hfile); return NULL; }

   ni_ver = nifti_header_version((char *)&n1hdr, h1size);

   if (ni_ver == 0 || ni_ver == 1) {
      nim = nifti_convert_n1hdr2nim(n1hdr, hfile);
      onefile = NIFTI_ONEFILE(n1hdr);
   } else if (ni_ver == 2) {
      memcpy(&n2hdr, &n1hdr, h1size);
      remain = h2size - h1size;
      posn = (char *)&n2hdr + h1size;
      ii = (int)nii_read(posn, 1, remain, fp);
      if (ii < (int)remain) { nii_close(&fp); free(hfile); return NULL; }
      nim = nifti_convert_n2hdr2nim(n2hdr, hfile);
      onefile = NIFTI_ONEFILE(n2hdr);
   } else {
      nii_close(&fp); free(hfile); return NULL;
   }

#ifdef REJECT_COMPLEX
   if (nim && (nim->datatype == DT_COMPLEX64 || nim->datatype == DT_COMPLEX128 || nim->datatype == DT_COMPLEX256)) {
      fprintf(stderr, "Image Exception Unsupported datatype (COMPLEX): use fslcomplex to manipulate: %s\n", hname);
      exit(13);
   }
#endif

   if (!nim) { nii_close(&fp); free(hfile); return NULL; }

   /* read extensions */
   if (onefile) remain = nim->iname_offset;
   else         remain = filesize;
   if (ni_ver <= 1) remain -= h1size;
   else             remain -= h2size;
   nifti_read_extensions(nim, fp, remain);
   nii_close(&fp);
   free(hfile);

   if (read_data) {
      if (isStdIn) {
         nim->iname = nifti_strdup("-");
         if (ni_ver != 1) { fprintf(stderr, "piped input only supports NIFTI-1\n"); return NULL; }
      }
      /* load image data */
      int64_t ntot = (int64_t)nim->nbyper * nim->nvox;
      NIIFILE dfp;
      if (isStdIn) {
         dfp = nii_open("-", "rb", 0);
      } else {
         char *imgname = nifti_findimgname(nim->iname, nim->nifti_type);
         if (!imgname) { nifti_image_free(nim); return NULL; }
         dfp = nii_open(imgname, "rb", nifti_is_gzfile(imgname));
         free(imgname);
      }
      if (!dfp) { nifti_image_free(nim); return NULL; }

      int64_t ioff = nim->iname_offset;
      if (isStdIn) ioff -= sizeof(nifti_1_header);
      if (nii_seek(dfp, (long)ioff, SEEK_SET) < 0) { nii_close(&dfp); nifti_image_free(nim); return NULL; }

      if (!nim->data) {
         nim->data = calloc(1, ntot);
         if (!nim->data) { nii_close(&dfp); nifti_image_free(nim); return NULL; }
      }
      int64_t nread = (int64_t)nii_read(nim->data, 1, ntot, dfp);
      if (nread < ntot) {
         fprintf(stderr, "++ WARNING: read %" PRId64 " of %" PRId64 " bytes\n", nread, ntot);
         nii_close(&dfp); free(nim->data); nim->data = NULL; nifti_image_free(nim); return NULL;
      }
      nii_close(&dfp);

      /* byte swap if needed */
      if (nim->swapsize > 1 && nim->byteorder != nifti_short_order())
         nifti_swap_Nbytes(ntot / nim->swapsize, nim->swapsize, nim->data);
      nim->byteorder = nifti_short_order();
   } else {
      nim->data = NULL;
   }

   return nim;
}

/*========== Image write ==========*/

static void nifti_set_iname_offset(nifti_image *nim, int nifti_ver)
{
   int64_t hsize = (nifti_ver == 2) ? (int64_t)sizeof(nifti_2_header) : (int64_t)sizeof(nifti_1_header);
   if (nim->nifti_type == NIFTI_FTYPE_NIFTI1_1 || nim->nifti_type == NIFTI_FTYPE_NIFTI2_1) {
      int64_t offset = nifti_extension_size(nim) + hsize + 4;
      if ((offset % 16) != 0) offset = ((offset + 0xf) & ~0xf);
      nim->iname_offset = offset;
   } else {
      nim->iname_offset = 0;
   }
}

static int isStdOutFcn(nifti_image *nim)
{
   if (!nim || !nim->fname) return 0;
   char *basename = nifti_makebasename(nim->fname);
   int isToStdOut = (basename && strcmp(basename, "-") == 0);
   free(basename);
   return isToStdOut;
}

#ifdef PIGZ
#ifdef HAVE_ZLIB
static int doPigz_write(nifti_image *nim, void *hdr, int hdrsize)
{
   FILE *pigzPipe;
   size_t cmdlen = strlen(nim->fname) + 32;
   char *command = (char *)malloc(cmdlen);
   if (!command) return -1;
   snprintf(command, cmdlen, "pigz -n -f > \"%s\"", nim->fname);
#ifdef _MSC_VER
   pigzPipe = _popen(command, "w");
#else
   pigzPipe = popen(command, "w");
#endif
   free(command);
   if (!pigzPipe) return -1;
   /* Write header */
   fwrite(hdr, 1, hdrsize, pigzPipe);
   /* Write extensions via wrapper */
   NIIFILE wfp = (NIIFILE)calloc(1, sizeof(*wfp));
   if (!wfp) { pclose(pigzPipe); return -1; }
   wfp->nzfptr = pigzPipe;
   if (nim->nifti_type != NIFTI_FTYPE_ANALYZE)
      nifti_write_extensions(wfp, nim);
   /* Write data */
   if (nim->data)
      fwrite(nim->data, 1, (size_t)nim->nbyper * nim->nvox, pigzPipe);
   free(wfp);
#ifdef _MSC_VER
   _pclose(pigzPipe);
#else
   pclose(pigzPipe);
#endif
   return 0;
}
#endif
#endif

void nifti_image_write(nifti_image *nim)
{
   nifti_1_header n1hdr;
   nifti_2_header n2hdr;
   NIIFILE fp;
   int nver = 1;
   int hsize = (int)sizeof(nifti_1_header);

   if (!nim || !nifti_validfilename(nim->fname)) return;
   if (!nim->data) return;

   int isStdOut = isStdOutFcn(nim);

   if (nim->nifti_type == NIFTI_FTYPE_NIFTI2_1 || nim->nifti_type == NIFTI_FTYPE_NIFTI2_2) {
      nifti_set_iname_offset(nim, 2);
      if (nifti_convert_nim2n2hdr(nim, &n2hdr)) return;
      nver = 2;
      hsize = (int)sizeof(nifti_2_header);
   } else {
      nifti_set_iname_offset(nim, 1);
      if (nifti_convert_nim2n1hdr(nim, &n1hdr)) return;
   }

   /* ensure iname is set for 2-file formats */
   if (nim->nifti_type != NIFTI_FTYPE_NIFTI1_1 && nim->nifti_type != NIFTI_FTYPE_NIFTI2_1) {
      if (nim->iname && strcmp(nim->iname, nim->fname) == 0) { free(nim->iname); nim->iname = NULL; }
      if (!nim->iname) {
         nim->iname = nifti_makeimgname(nim->fname, nim->nifti_type, 0, 0);
         if (!nim->iname) return;
      }
   }

   /* try PIGZ if applicable */
   if (!isStdOut) {
#ifdef PIGZ
#ifdef HAVE_ZLIB
      if ((nim->nifti_type == NIFTI_FTYPE_NIFTI1_1 || nim->nifti_type == NIFTI_FTYPE_NIFTI2_1)
          && nifti_is_gzfile(nim->fname) == 1) {
         const char *val = getenv("AFNI_COMPRESSOR");
         if (val && strstr(val, "PIGZ")) {
            void *hdr = (nver == 2) ? (void *)&n2hdr : (void *)&n1hdr;
            if (doPigz_write(nim, hdr, hsize) == 0) return;
         }
      }
#endif
#endif
   }

   /* open output file */
   if (isStdOut) {
#ifdef _WIN32
      _setmode(_fileno(stdout), _O_BINARY);
#endif
      fp = nii_open("-", "wb", 0);
   } else {
      fp = nii_open(nim->fname, "wb", nifti_is_gzfile(nim->fname));
   }
   if (!fp) return;

   /* write header */
   size_t ss;
   if (nver == 2) ss = nii_write(&n2hdr, 1, hsize, fp);
   else           ss = nii_write(&n1hdr, 1, hsize, fp);
   if ((int)ss < hsize) { nii_close(&fp); return; }

   /* write extensions */
   if (nim->nifti_type != NIFTI_FTYPE_ANALYZE)
      nifti_write_extensions(fp, nim);

   /* for 2-file format, close header file and open image file */
   if (nim->nifti_type != NIFTI_FTYPE_NIFTI1_1 && nim->nifti_type != NIFTI_FTYPE_NIFTI2_1) {
      nii_close(&fp);
      fp = nii_open(nim->iname, "wb", nifti_is_gzfile(nim->iname));
      if (!fp) return;
   }

   /* seek to data offset and write data */
   nii_seek(fp, (long)nim->iname_offset, SEEK_SET);
   int64_t ntot = (int64_t)nim->nbyper * nim->nvox;
   nii_write(nim->data, 1, ntot, fp);
   nim->byteorder = nifti_short_order();
   nii_close(&fp);
}

/*========== Free / infodump ==========*/

static int nifti_free_extensions(nifti_image *nim)
{
   if (!nim) return -1;
   if (nim->num_ext > 0 && nim->ext_list) {
      for (int c = 0; c < nim->num_ext; c++)
         if (nim->ext_list[c].edata) free(nim->ext_list[c].edata);
      free(nim->ext_list);
   }
   nim->num_ext = 0;
   nim->ext_list = NULL;
   return 0;
}

void nifti_image_free(nifti_image *nim)
{
   if (!nim) return;
   if (nim->fname) free(nim->fname);
   if (nim->iname) free(nim->iname);
   if (nim->data)  free(nim->data);
   nifti_free_extensions(nim);
   free(nim);
}

void nifti_image_infodump(const nifti_image *nim)
{
   if (!nim) return;
   fprintf(stderr, "nifti_image_infodump:\n");
   fprintf(stderr, "  fname = %s\n", nim->fname ? nim->fname : "(null)");
   fprintf(stderr, "  iname = %s\n", nim->iname ? nim->iname : "(null)");
   fprintf(stderr, "  ndim = %" PRId64 ", dim =", nim->ndim);
   for (int i = 0; i <= nim->ndim; i++) fprintf(stderr, " %" PRId64, nim->dim[i]);
   fprintf(stderr, "\n  datatype = %d, nbyper = %d, nvox = %" PRId64 "\n",
           nim->datatype, nim->nbyper, nim->nvox);
   fprintf(stderr, "  pixdim =");
   for (int i = 1; i <= nim->ndim; i++) fprintf(stderr, " %g", nim->pixdim[i]);
   fprintf(stderr, "\n  scl_slope = %g, scl_inter = %g\n", nim->scl_slope, nim->scl_inter);
   fprintf(stderr, "  cal_min = %g, cal_max = %g\n", nim->cal_min, nim->cal_max);
   fprintf(stderr, "  qform_code = %d, sform_code = %d\n", nim->qform_code, nim->sform_code);
   fprintf(stderr, "  descrip = %.80s\n", nim->descrip);
   fprintf(stderr, "  iname_offset = %" PRId64 "\n", nim->iname_offset);
   fprintf(stderr, "  nifti_type = %d\n", nim->nifti_type);
}
