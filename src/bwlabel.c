#include <stdlib.h>
#include <math.h>
#include <nifti2_io.h>
#include "core.h"
#include "print.h"
#ifdef __aarch64__
  #include "arm_malloc.h"
#else
  #include <immintrin.h>
#endif

//Jesper Andersson has acknowledged that this port of spm_bwlabel.c may be released using the BSD 2-Clause license

//Copyright 2021 Jesper Andersson
//Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

//for usage, see https://en.wikibooks.org/wiki/SPM/How-to#How_to_remove_clusters_under_a_certain_size_in_a_binary_mask.3F

/****************************************************************
 **
 ** Set of routines implementing a 2D or 3D connected component
 ** labelling algorithm. Its interface is modelled on bwlabel
 ** (which is a routine in the image processing toolbox) and
 ** takes as input a binary image and (optionally) a connectednes
 ** criterion (6, 18 or 26, in 2D 6 will correspond to 4 and 18
 ** and 26 to 8). It will output an image/volume where each 
 ** connected component will have a unique label.
 **
 ** The implementation is not recursive (i.e. will no crash for
 ** large connected components) and is loosely based on
 ** Thurfjell et al. 1992, A new three-dimensional connected
 ** components labeling algorithm with simultaneous object
 ** feature extraction capability. CVGIP: Graphical Models 
 ** and Image Processing 54(4):357-364.
 **
 ***************************************************************/

void fill_tratab(uint32_t  *tt,     /* Translation table */
                 /*uint32_t  ttn,*/     /* Size of translation table */
                 uint32_t  *nabo,   /* Set of neighbours */
                 uint32_t  nr_set)  /* Number of neighbours in nabo */
{
   int           i = 0, j = 0, cntr = 0;
   uint32_t  tn[9];
   uint32_t  ltn = UINT_MAX;

   /*
   Find smallest terminal number in neighbourhood
   */

   for (i=0; i<nr_set; i++)
   {
      j = nabo[i];
      cntr=0;
      while (tt[j-1] != j) 
      {
         j = tt[j-1];
         cntr++;
         if (cntr>100) {printf("\nOoh no!!"); break;}
      }
      tn[i] = j;
      ltn = MIN(ltn,j);
   }
   /*
   Replace all terminal numbers in neighbourhood by the smallest one
   */
   for (i=0; i<nr_set; i++)
   {
      tt[tn[i]-1] = ltn;
   }

   return;
}

#define idx(A,B,C,DIM) ((C)*DIM[0]*DIM[1] + (B)*DIM[0] + (A))

uint32_t check_previous_slice(uint32_t  *il,     /* Initial labelling map */
                                  uint32_t  r,       /* row */
                                  uint32_t  c,       /* column */
                                  uint32_t  sl,      /* slice */
                                  size_t        dim[3],  /* dimensions of il */
                                  uint32_t  conn,    /* Connectivity criterion */
                                  uint32_t  *tt)     /* Translation table */
//                                  uint32_t  ttn)     /* Size of translation table */
{
   uint32_t l=0;
   uint32_t nabo[9];
   uint32_t nr_set = 0;

   if (!sl) return(0);
  
   if (conn >= 6)
   {
      if ((l = il[idx(r,c,sl-1,dim)])) {nabo[nr_set++] = l;}
   }
   if (conn >= 18)
   {
      if (r) {if ((l = il[idx(r-1,c,sl-1,dim)])) {nabo[nr_set++] = l;}}
      if (c) {if ((l = il[idx(r,c-1,sl-1,dim)])) {nabo[nr_set++] = l;}}
      if (r < dim[0]-1) {if ((l = il[idx(r+1,c,sl-1,dim)])) {nabo[nr_set++] = l;}}
      if (c < dim[1]-1) {if ((l = il[idx(r,c+1,sl-1,dim)])) {nabo[nr_set++] = l;}}
   }
   if (conn == 26)
   {
      if (r && c) {if ((l = il[idx(r-1,c-1,sl-1,dim)])) {nabo[nr_set++] = l;}}
      if ((r < dim[0]-1) && c) {if ((l = il[idx(r+1,c-1,sl-1,dim)])) {nabo[nr_set++] = l;}}
      if (r && (c < dim[1]-1)) {if ((l = il[idx(r-1,c+1,sl-1,dim)])) {nabo[nr_set++] = l;}}
      if ((r < dim[0]-1) && (c < dim[1]-1)) {if ((l = il[idx(r+1,c+1,sl-1,dim)])) {nabo[nr_set++] = l;}}
   }

   if (nr_set) 
   {
      fill_tratab(tt,/*ttn,*/nabo,nr_set);
      return(nabo[0]);
   }
   else {return(0);}
}

void * mxRealloc(void *oldArray, size_t oldBytes, size_t newBytes) {
   // https://octave.org/doxygen/3.8/df/d4e/mex_8cc_source.html
   //reallocate memory, preserve previous bytes
   if (newBytes <= 0) {
      _mm_free(oldArray);
      return NULL;
   }
   void *newArray = (void *)_mm_malloc(newBytes, 64);
   memset(newArray, 0, newBytes);
   if (oldBytes > 0) {
     //void * memcpy ( void * destination, const void * source, size_t num );
     oldBytes = MIN(oldBytes, newBytes);
     //void * memcpy ( void * destination, const void * source, size_t num );
     memcpy(newArray, oldArray, oldBytes);
     _mm_free(oldArray);
   }
   return newArray;
}

/* do_initial_labelling */

uint32_t do_initial_labelling(uint8_t        *bw,   /* Binary map */
                                  size_t        *dim,  /* Dimensions of bw */
                                  uint32_t  conn,  /* Connectivity criterion */
                                  uint32_t  *il,   /* Initially labelled map */
                                  uint32_t  **tt)  /* Translation table */
{
   uint32_t  i = 0, j = 0;
   uint32_t  nabo[8];
   uint32_t  label = 1;
   uint32_t  nr_set = 0;
   uint32_t  l = 0;
   int32_t       sl, r, c;
   uint32_t  ttn = 1000;
   *tt = (uint32_t *)_mm_malloc(ttn * sizeof(uint32_t), 64);
   memset(*tt, 0, ttn * sizeof(uint32_t));
   for (sl=0; sl<dim[2]; sl++)
   {
      for (c=0; c<dim[1]; c++)
      {
         for (r=0; r<dim[0]; r++)
         {
            nr_set = 0;
            if (bw[idx(r,c,sl,dim)])
            {
               nabo[0] = check_previous_slice(il,r,c,sl,dim,conn,*tt /*,ttn*/);
               if (nabo[0]) {nr_set += 1;}
               /*
                  For six(surface)-connectivity
               */
               if (conn >= 6)
               {
                  if (r)
                  {
                     if ((l = il[idx(r-1,c,sl,dim)])) {nabo[nr_set++] = l;}
                  }
                  if (c)
                  {
                     if ((l = il[idx(r,c-1,sl,dim)])) {nabo[nr_set++] = l;}
                  }
               }
               /*
                  For 18(edge)-connectivity
                  N.B. In current slice no difference to 26.
               */
               if (conn >= 18)
               {
                  if (c && r)
                  {
                     if ((l = il[idx(r-1,c-1,sl,dim)])) {nabo[nr_set++] = l;}
                  }
                  if (c && (r < dim[0]-1))
                  {
                     if ((l = il[idx(r+1,c-1,sl,dim)])) {nabo[nr_set++] = l;}
                  }
               }
               if (nr_set)
               {
                  il[idx(r,c,sl,dim)] = nabo[0];
                  fill_tratab(*tt,/*ttn,*/nabo,nr_set);
               }
               else
               {
                  il[idx(r,c,sl,dim)] = label;
                  if (label >= ttn) {ttn += 1000; *tt = mxRealloc(*tt, (ttn - 1000)*sizeof(uint32_t), ttn*sizeof(uint32_t));}
                  (*tt)[label-1] = label;
                  label++;
               }
            }
         }
      }
   }

   /*
      Finalise translation table
   */

   for (i=0; i<(label-1); i++)
   {
      j = i;
      while ((*tt)[j] != j+1)
      {
         j = (*tt)[j]-1;
      }
      (*tt)[i] = j+1;
   }
 
   
   return(label-1);
}

/* translate_labels */

double translate_labels(uint32_t  *il,     /* Map of initial labels. */
                        size_t        dim[3],  /* Dimensions of il. */
                        uint32_t  *tt,     /* Translation table. */
                        uint32_t  ttn,     /* Size of translation table. */
                        double        *l)      /* Final map of labels. */
{
   int            n=0;
   int            i=0;
   uint32_t   ml=0;
   double         cl = 0.0;
   n = dim[0]*dim[1]*dim[2];
   for (i=0; i<ttn; i++) {ml = MAX(ml,tt[i]);}
   double *fl = (double *)_mm_malloc(ml * sizeof(double), 64); 
   memset(fl, 0, ml * sizeof(double));
   for (i=0; i<n; i++)
   {
      if (il[i])
      {
         if (!fl[tt[il[i]-1]-1]) 
         {
            cl += 1.0; 
            fl[tt[il[i]-1]-1] = cl;
         }
         l[i] = fl[tt[il[i]-1]-1];
      }
   }
   _mm_free(fl);
   return(cl);
}


int bwlabel(nifti_image *nim, int conn) {
  if ((conn!=6) && (conn!=18) && (conn!=26)) {
     printfx("bwlabel: conn must be 6, 18 or 26.\n");
     return 11;
  }
  if ((nim->ndim < 2) || (nim->ndim > 3)) {
    printfx("bwlabel: img must be 2 or 3-dimensional\n");
    return 33;
  }
  size_t dim[3];
  dim[0]=nim->nx;
  dim[1]=nim->ny;
  dim[2]=nim->nz;
  size_t nvox = nim->nx*nim->ny*nim->nz;
  double *l = (double *)_mm_malloc(nvox * sizeof(double), 64); //output image
  memset(l, 0, nvox * sizeof(double));
  uint32_t *il = (uint32_t *)_mm_malloc(nvox * sizeof(uint32_t), 64);
  memset(il, 0, nvox * sizeof(uint32_t));
  uint8_t *bw = (uint8_t *)_mm_malloc(nvox * sizeof(uint8_t), 64);
  memset(bw, 0, nvox * sizeof(uint8_t));
  if (nim->datatype == DT_FLOAT32) {
    float *img = (float *)nim->data;
    for (size_t i = 0; i < nvox; i++)
      if (img[i] != 0.0) bw[i] = 1;
  } else if (nim->datatype == DT_FLOAT64) {
    double *img = (double *)nim->data;
    for (size_t i = 0; i < nvox; i++)
      if (img[i] != 0.0) bw[i] = 1;
  } else {
    printfx("bwlabel: Unsupported datatype %d\n", nim->datatype);
    return 22;
  }
  uint32_t  *tt = NULL;
  uint32_t ttn = do_initial_labelling(bw,dim,conn,il,&tt);
  double nl = translate_labels(il,dim,tt,ttn,l);
  nim->cal_min = 0;
  nim->cal_max = nl;
  nim->scl_inter = 0.0;
  nim->scl_slope = 1.0;
  _mm_free(il);
  _mm_free(tt);
  if (nim->datatype == DT_FLOAT32) {
    float *img = (float *)nim->data;
    for (size_t i = 0; i < nvox; i++)
      img[i] = l[i];
  } else if (nim->datatype == DT_FLOAT32) {
    double *img = (double *)nim->data;
    for (size_t i = 0; i < nvox; i++)
      img[i] = l[i];
  }
  _mm_free(l);
  return(EXIT_SUCCESS);
}