/* nifti_io.h - Minimal NIfTI I/O for niimath
 *
 * Replaces niftilib (nifti2_io.h, nifti1.h, nifti2.h) and znzlib.
 * Public domain code based on work by Robert W Cox, Mark Jenkinson,
 * Rick Reynolds, and Chris Rorden.
 */

#ifndef NIFTI_IO_H
#define NIFTI_IO_H

#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#ifdef __cplusplus
extern "C" {
#endif

/*========== NIfTI-1 header struct (348 bytes) ==========*/

#pragma pack(push)
#pragma pack(1)

struct nifti_1_header {
   int32_t sizeof_hdr;
   char    data_type[10];
   char    db_name[18];
   int32_t extents;
   int16_t session_error;
   char    regular;
   char    dim_info;
   int16_t dim[8];
   float   intent_p1;
   float   intent_p2;
   float   intent_p3;
   int16_t intent_code;
   int16_t datatype;
   int16_t bitpix;
   int16_t slice_start;
   float   pixdim[8];
   float   vox_offset;
   float   scl_slope;
   float   scl_inter;
   int16_t slice_end;
   char    slice_code;
   char    xyzt_units;
   float   cal_max;
   float   cal_min;
   float   slice_duration;
   float   toffset;
   int32_t glmax;
   int32_t glmin;
   char    descrip[80];
   char    aux_file[24];
   int16_t qform_code;
   int16_t sform_code;
   float   quatern_b;
   float   quatern_c;
   float   quatern_d;
   float   qoffset_x;
   float   qoffset_y;
   float   qoffset_z;
   float   srow_x[4];
   float   srow_y[4];
   float   srow_z[4];
   char    intent_name[16];
   char    magic[4];
};
typedef struct nifti_1_header nifti_1_header;

/*========== NIfTI-2 header struct (540 bytes) ==========*/

struct nifti_2_header {
   int32_t sizeof_hdr;
   char    magic[8];
   int16_t datatype;
   int16_t bitpix;
   int64_t dim[8];
   double  intent_p1;
   double  intent_p2;
   double  intent_p3;
   double  pixdim[8];
   int64_t vox_offset;
   double  scl_slope;
   double  scl_inter;
   double  cal_max;
   double  cal_min;
   double  slice_duration;
   double  toffset;
   int64_t slice_start;
   int64_t slice_end;
   char    descrip[80];
   char    aux_file[24];
   int32_t qform_code;
   int32_t sform_code;
   double  quatern_b;
   double  quatern_c;
   double  quatern_d;
   double  qoffset_x;
   double  qoffset_y;
   double  qoffset_z;
   double  srow_x[4];
   double  srow_y[4];
   double  srow_z[4];
   int32_t slice_code;
   int32_t xyzt_units;
   int32_t intent_code;
   char    intent_name[16];
   char    dim_info;
   char    unused_str[15];
};
typedef struct nifti_2_header nifti_2_header;

#pragma pack(pop)

/*========== Matrix types ==========*/

typedef struct { float  m[4][4]; } mat44;
typedef struct { float  m[3][3]; } mat33;
typedef struct { double m[4][4]; } nifti_dmat44;
typedef struct { double m[3][3]; } nifti_dmat33;

/*========== Extension type ==========*/

typedef struct {
   int    esize;
   int    ecode;
   char  *edata;
} nifti1_extension;

typedef struct {
   unsigned char extension[4];
} nifti1_extender;

/*========== Analyze 7.5 orient codes ==========*/

typedef enum _analyze75_orient_code {
   a75_transverse_unflipped = 0,
   a75_coronal_unflipped    = 1,
   a75_sagittal_unflipped   = 2,
   a75_transverse_flipped   = 3,
   a75_coronal_flipped      = 4,
   a75_sagittal_flipped     = 5,
   a75_orient_unknown       = 6
} analyze_75_orient_code;

/*========== Analyze 7.5 header struct (for swap_nifti_header) ==========*/

typedef struct {
   int32_t sizeof_hdr;
   char    data_type[10];
   char    db_name[18];
   int32_t extents;
   int16_t session_error;
   char    regular;
   char    hkey_un0;
   int16_t dim[8];
   int16_t unused8;
   int16_t unused9;
   int16_t unused10;
   int16_t unused11;
   int16_t unused12;
   int16_t unused13;
   int16_t unused14;
   int16_t datatype;
   int16_t bitpix;
   int16_t dim_un0;
   float   pixdim[8];
   float   vox_offset;
   float   funused1;
   float   funused2;
   float   funused3;
   float   cal_max;
   float   cal_min;
   float   compressed;
   float   verified;
   int32_t glmax;
   int32_t glmin;
   char    descrip[80];
   char    aux_file[24];
   char    orient;
   char    originator[10];
   char    generated[10];
   char    scannum[10];
   char    patient_id[10];
   char    exp_date[10];
   char    exp_time[10];
   char    hist_un0[3];
   int32_t views;
   int32_t vols_added;
   int32_t start_field;
   int32_t field_skip;
   int32_t omax;
   int32_t omin;
   int32_t smax;
   int32_t smin;
} nifti_analyze75;

/*========== nifti_image struct ==========*/

typedef struct {
   int64_t ndim;
   int64_t nx, ny, nz, nt, nu, nv, nw;
   int64_t dim[8];
   int64_t nvox;
   int     nbyper;
   int     datatype;
   double  dx, dy, dz, dt, du, dv, dw;
   double  pixdim[8];
   double  scl_slope;
   double  scl_inter;
   double  cal_min;
   double  cal_max;
   int     qform_code;
   int     sform_code;
   double  quatern_b, quatern_c, quatern_d;
   double  qoffset_x, qoffset_y, qoffset_z;
   double  qfac;
   nifti_dmat44 qto_xyz;
   nifti_dmat44 qto_ijk;
   nifti_dmat44 sto_xyz;
   nifti_dmat44 sto_ijk;
   double  toffset;
   int     xyz_units;
   int     time_units;
   int     nifti_type;
   int     intent_code;
   double  intent_p1;
   double  intent_p2;
   double  intent_p3;
   char    intent_name[16];
   char    descrip[80];
   char    aux_file[24];
   char   *fname;
   char   *iname;
   int64_t iname_offset;
   int     swapsize;
   int     byteorder;
   void   *data;
   int     num_ext;
   nifti1_extension *ext_list;
   analyze_75_orient_code analyze75_orient;
   int     freq_dim;
   int     phase_dim;
   int     slice_dim;
   int     slice_code;
   int64_t slice_start;
   int64_t slice_end;
   double  slice_duration;
} nifti_image;

/*========== Constants ==========*/

/* NIfTI file types */
#define NIFTI_FTYPE_ANALYZE   0
#define NIFTI_FTYPE_NIFTI1_1  1
#define NIFTI_FTYPE_NIFTI1_2  2
#define NIFTI_FTYPE_ASCII     3
#define NIFTI_FTYPE_NIFTI2_1  4
#define NIFTI_FTYPE_NIFTI2_2  5
#define NIFTI_MAX_FTYPE       5

/* Byte order */
#define LSB_FIRST 1
#define MSB_FIRST 2
#define REVERSE_ORDER(x) ((x)==LSB_FIRST ? MSB_FIRST : LSB_FIRST)

/* Datatype codes */
#define DT_UNKNOWN     0
#define DT_BINARY      1
#define DT_UINT8       2
#define DT_INT16       4
#define DT_INT32       8
#define DT_FLOAT32    16
#define DT_COMPLEX64  32
#define DT_FLOAT64    64
#define DT_RGB24     128
#define DT_INT8      256
#define DT_UINT16    512
#define DT_UINT32    768
#define DT_INT64    1024
#define DT_UINT64   1280
#define DT_FLOAT128 1536
#define DT_COMPLEX128 1792
#define DT_COMPLEX256 2048
#define DT_RGBA32   2304

/* Aliases */
#define DT_UNSIGNED_CHAR DT_UINT8
#define DT_SIGNED_SHORT  DT_INT16
#define DT_SIGNED_INT    DT_INT32
#define DT_FLOAT         DT_FLOAT32
#define DT_COMPLEX       DT_COMPLEX64
#define DT_DOUBLE        DT_FLOAT64
#define DT_RGB           DT_RGB24
#define DT_ALL           255

#define NIFTI_TYPE_UINT8       2
#define NIFTI_TYPE_INT16       4
#define NIFTI_TYPE_INT32       8
#define NIFTI_TYPE_FLOAT32    16
#define NIFTI_TYPE_COMPLEX64  32
#define NIFTI_TYPE_FLOAT64    64
#define NIFTI_TYPE_RGB24     128
#define NIFTI_TYPE_INT8      256
#define NIFTI_TYPE_UINT16    512
#define NIFTI_TYPE_UINT32    768
#define NIFTI_TYPE_INT64    1024
#define NIFTI_TYPE_UINT64   1280
#define NIFTI_TYPE_FLOAT128 1536
#define NIFTI_TYPE_COMPLEX128 1792
#define NIFTI_TYPE_COMPLEX256 2048
#define NIFTI_TYPE_RGBA32   2304

/* Intent codes */
#define NIFTI_INTENT_NONE          0
#define NIFTI_INTENT_CORREL        2
#define NIFTI_INTENT_TTEST         3
#define NIFTI_INTENT_FTEST         4
#define NIFTI_INTENT_ZSCORE        5
#define NIFTI_INTENT_CHISQ         6
#define NIFTI_INTENT_BETA          7
#define NIFTI_INTENT_BINOM         8
#define NIFTI_INTENT_GAMMA         9
#define NIFTI_INTENT_POISSON      10
#define NIFTI_INTENT_NORMAL       11
#define NIFTI_INTENT_FTEST_NONC   12
#define NIFTI_INTENT_CHISQ_NONC   13
#define NIFTI_INTENT_LOGISTIC     14
#define NIFTI_INTENT_LAPLACE      15
#define NIFTI_INTENT_UNIFORM      16
#define NIFTI_INTENT_TTEST_NONC   17
#define NIFTI_INTENT_WEIBULL      18
#define NIFTI_INTENT_CHI          19
#define NIFTI_INTENT_INVGAUSS     20
#define NIFTI_INTENT_EXTVAL       21
#define NIFTI_INTENT_PVAL         22
#define NIFTI_INTENT_LOGPVAL      23
#define NIFTI_INTENT_LOG10PVAL    24
#define NIFTI_INTENT_ESTIMATE    1001
#define NIFTI_INTENT_LABEL       1002
#define NIFTI_INTENT_NEURONAME   1003
#define NIFTI_INTENT_GENMATRIX   1004
#define NIFTI_INTENT_SYMMATRIX   1005
#define NIFTI_INTENT_DISPVECT    1006
#define NIFTI_INTENT_VECTOR      1007
#define NIFTI_INTENT_POINTSET    1008
#define NIFTI_INTENT_TRIANGLE    1009
#define NIFTI_INTENT_QUATERNION  1010
#define NIFTI_INTENT_DIMLESS     1011

#define NIFTI_FIRST_STATCODE  2
#define NIFTI_LAST_STATCODE  24

/* Transform codes */
#define NIFTI_XFORM_UNKNOWN       0
#define NIFTI_XFORM_SCANNER_ANAT  1
#define NIFTI_XFORM_ALIGNED_ANAT  2
#define NIFTI_XFORM_TALAIRACH     3
#define NIFTI_XFORM_MNI_152       4

/* Units codes */
#define NIFTI_UNITS_UNKNOWN  0
#define NIFTI_UNITS_METER    1
#define NIFTI_UNITS_MM       2
#define NIFTI_UNITS_MICRON   3
#define NIFTI_UNITS_SEC      8
#define NIFTI_UNITS_MSEC    16
#define NIFTI_UNITS_USEC    24
#define NIFTI_UNITS_HZ      32
#define NIFTI_UNITS_PPM     40
#define NIFTI_UNITS_RADS    48

/* Slice order codes */
#define NIFTI_SLICE_UNKNOWN    0
#define NIFTI_SLICE_SEQ_INC    1
#define NIFTI_SLICE_SEQ_DEC    2
#define NIFTI_SLICE_ALT_INC    3
#define NIFTI_SLICE_ALT_DEC    4
#define NIFTI_SLICE_ALT_INC2   5
#define NIFTI_SLICE_ALT_DEC2   6

/* Extension codes */
#define NIFTI_ECODE_IGNORE     0
#define NIFTI_ECODE_DICOM      2
#define NIFTI_ECODE_AFNI       4
#define NIFTI_ECODE_COMMENT    6
#define NIFTI_ECODE_XCEDE      8
#define NIFTI_ECODE_JIMDIMINFO 10
#define NIFTI_ECODE_WORKFLOW_FWDS 12
#define NIFTI_ECODE_FREESURFER 14
#define NIFTI_ECODE_PYPICKLE   16
#define NIFTI_ECODE_MIND_IDENT 18
#define NIFTI_ECODE_B_VALUE    20
#define NIFTI_ECODE_SPHERICAL_DIRECTION 22
#define NIFTI_ECODE_DT_COMPONENT 24
#define NIFTI_ECODE_SHC_DEGREEORDER 26
#define NIFTI_ECODE_VOXBO      28
#define NIFTI_ECODE_CARET      30
#define NIFTI_ECODE_CIFTI      32
#define NIFTI_ECODE_VARIABLE_FRAME_TIMING 34
#define NIFTI_ECODE_EVAL       38
#define NIFTI_ECODE_MATLAB     40
#define NIFTI_ECODE_QUANTIPHYSE 42
#define NIFTI_ECODE_MRS        44
#define NIFTI_MAX_ECODE        44

/*========== Macros ==========*/

/* Check NIfTI version from magic field */
#define NIFTI_VERSION(h) \
   (( (h).magic[0]=='n' && (h).magic[3]=='\0' && \
     ((h).magic[1]=='i' || (h).magic[1]=='+') && \
     ((h).magic[2]>='1' && (h).magic[2]<='9') ) ? (h).magic[2]-'0' : 0)

/* Check if single file (.nii) */
#define NIFTI_ONEFILE(h) ((h).magic[1] == '+')

/* Swap test for NIfTI-2 */
#define NIFTI2_NEEDS_SWAP(h) \
   ((h).sizeof_hdr == 1543569408 || (h).sizeof_hdr == 469893120)

/* xyzt_units packing/unpacking */
#define XYZT_TO_SPACE(xyzt) ((xyzt) & 0x07)
#define XYZT_TO_TIME(xyzt)  ((xyzt) & 0x38)
#define SPACE_TIME_TO_XYZT(ss,tt) ((ss) | (tt))

/* dim_info packing/unpacking */
#define DIM_INFO_TO_FREQ_DIM(di)   ((di)     & 0x03)
#define DIM_INFO_TO_PHASE_DIM(di)  (((di)>>2) & 0x03)
#define DIM_INFO_TO_SLICE_DIM(di)  (((di)>>4) & 0x03)
#define FPS_INTO_DIM_INFO(f,p,s)   ((f & 0x03) | ((p & 0x03)<<2) | ((s & 0x03)<<4))

/* NaN/Inf check */
#define IS_GOOD_FLOAT(x) (isfinite(x))
#define FIXED_FLOAT(x)   (IS_GOOD_FLOAT(x) ? (x) : 0.0)

/* 16-bit int range check */
#define NIFTI_IS_16_BIT_INT(x) ((x) >= -32768 && (x) <= 32767)

/* Assign if pointer not NULL */
#undef ASSIF
#define ASSIF(p,v)  if((p)!=NULL) *(p) = (v)

/*========== Public API ==========*/

/* I/O */
nifti_image *nifti_image_read(const char *hname, int read_data);
void         nifti_image_write(nifti_image *nim);
void         nifti_image_free(nifti_image *nim);
void         nifti_image_infodump(const nifti_image *nim);

/* File utilities */
char        *nifti_find_file_extension(const char *name);
int          nifti_set_filenames(nifti_image *nim, const char *prefix,
                                 int check, int set_byte_order);
int          nifti_is_gzfile(const char *fname); /* 0=none, 1=gzip(.gz), 2=zstd(.zst) */
void         nifti_datatype_sizes(int datatype, int *nbyper, int *swapsize);
int          nifti_short_order(void);
int          nifti_compiled_with_zlib(void);

/* Math - float precision */
mat44        nifti_mat44_inverse(mat44 R);
mat44        nifti_mat44_mul(mat44 A, mat44 B);
mat44        nifti_quatern_to_mat44(float qb, float qc, float qd,
                                     float qx, float qy, float qz,
                                     float dx, float dy, float dz, float qfac);
void         nifti_mat44_to_quatern(mat44 R,
                                     float *qb, float *qc, float *qd,
                                     float *qx, float *qy, float *qz,
                                     float *dx, float *dy, float *dz, float *qfac);
float        nifti_mat33_determ(mat33 R);

/* Math - double precision */
nifti_dmat44 nifti_dmat44_inverse(nifti_dmat44 R);
nifti_dmat44 nifti_quatern_to_dmat44(double qb, double qc, double qd,
                                      double qx, double qy, double qz,
                                      double dx, double dy, double dz, double qfac);
void         nifti_dmat44_to_quatern(nifti_dmat44 R,
                                      double *qb, double *qc, double *qd,
                                      double *qx, double *qy, double *qz,
                                      double *dx, double *dy, double *dz, double *qfac);

/* Byte swapping */
void nifti_swap_2bytes(int64_t n, void *ar);
void nifti_swap_4bytes(int64_t n, void *ar);
void nifti_swap_8bytes(int64_t n, void *ar);
void nifti_swap_16bytes(int64_t n, void *ar);
void nifti_swap_Nbytes(int64_t n, int siz, void *ar);

/* String duplicate (internal use) */
char *nifti_strdup(const char *str);

#ifdef __cplusplus
}
#endif

#endif /* NIFTI_IO_H */
