#ifndef CR_COREX_H
#define CR_COREX_H

#ifdef  __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <nifti2_io.h>
#include <float.h> //FLT_EPSILON
//#include <immintrin.h>
#include <limits.h>

//CORE32 and CORE64 handle Float32 and Float64 operations, CORE handles shared code 

//store inputs regarding input header
typedef struct {
	int datatype;
	float scl_slope, scl_inter; 
} in_hdr;
//filter details for image resizing
//C Code From Graphics Gems III
//Dale Schumacher, "Optimization of Bitmap Scaling Operations"
//Dale Schumacher, "General Filtered Image Rescaling"
// https://github.com/erich666/GraphicsGems/tree/master/gemsiii
//see github for other filters

typedef struct {
	int	pixel;
	double	weight;
} CONTRIB;
typedef struct {
	int	n;		/* number of contributors */
	CONTRIB	*p;		/* pointer to list of contributions */
} CLIST;

typedef struct {                   /** x4 vector struct **/
    float v[4] ;
} vec4 ;

int nifti_save(nifti_image * nim, const char *postfix);
vec4 setVec4(float x, float y, float z);
vec4 nifti_vect44mat44_mul(vec4 v, mat44 m );
mat44 xform(nifti_image * nim);
int neg_determ(nifti_image * nim);
nifti_image *nifti_image_read2( const char *hname , int read_data );
float max_displacement_mm( nifti_image * nim,  nifti_image * nim2);
float vertexDisplacement(float x, float y, float z, mat44 m, mat44 m2);
in_hdr set_input_hdr(nifti_image * nim);
int nifti_image_change_datatype ( nifti_image * nim, int dt , in_hdr * ihdr);
int * make_kernel_file(nifti_image * nim, int * nkernel,  char * fin);
int * make_kernel(nifti_image * nim, int * nkernel, int x, int y, int z);
int * make_kernel_sphere(nifti_image * nim, int * nkernel, double mm);
//CLIST	* createFilter(int srcXsize, int dstXsize, double (*filterf)(), double fwidth);
CLIST	* createFilter(int srcXsize, int dstXsize, int filterMethod);
double qginv( double p );
double qg( double x );

#ifndef MAX //from Christian Gaser's TFCE example
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#ifndef MIN
#define MIN(A,B) ((A) > (B) ? (B) : (A))
#endif

#ifndef sqr
#define sqr(x) ((x)*(x))
#endif

#define LOAD_MAT44(AA,a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,a34)    \
( AA.m[0][0]=a11 , AA.m[0][1]=a12 , AA.m[0][2]=a13 , AA.m[0][3]=a14 ,   \
AA.m[1][0]=a21 , AA.m[1][1]=a22 , AA.m[1][2]=a23 , AA.m[1][3]=a24 ,   \
AA.m[2][0]=a31 , AA.m[2][1]=a32 , AA.m[2][2]=a33 , AA.m[2][3]=a34 ,   \
AA.m[3][0]=AA.m[3][1]=AA.m[3][2]=0.0f , AA.m[3][3]=1.0f            )

#define LOAD_MAT33(AA,a11,a12,a13 ,a21,a22,a23 ,a31,a32,a33)    \
( AA.m[0][0]=a11 , AA.m[0][1]=a12 , AA.m[0][2]=a13 ,   \
AA.m[1][0]=a21 , AA.m[1][1]=a22 , AA.m[1][2]=a23  ,   \
AA.m[2][0]=a31 , AA.m[2][1]=a32 , AA.m[2][2]=a33            )

#ifdef  __cplusplus
}
#endif

#endif // CR_COREX_H
