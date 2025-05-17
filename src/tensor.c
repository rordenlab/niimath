#include "math.h"
#include "tensor.h"
#include "print.h"
#include <stdlib.h>
#include <stdio.h>

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#ifndef MIN
#define MIN(A,B) ((A) > (B) ? (B) : (A))
#endif

#define old
  //Old method uses conversion of Fortran code. It is hard to read but faster than C
#ifdef old

#undef  EPS
#undef  EPSQ
#define EPS  1.e-8
#define EPSQ 1.e-4   /* sqrt(EPS) */
#undef  SQR
#define SQR(a)  ((a)*(a))
#undef  PI
#define PI 3.14159265358979323846
#undef  SWAP
#define SWAP(x,y) (th=(x),(x)=(y),(y)=th)
#undef  CSWAP
#define CSWAP(i,j) (SWAP(a[i],a[j]),SWAP(a[i+1],a[j+1]),SWAP(a[i+2],a[j+2]))

#undef  DET3
#define DET3(m) ( m[0]*m[4]*m[8]-m[0]*m[7]*m[5]-m[1]*m[3]*m[8] \
                 +m[1]*m[6]*m[5]+m[2]*m[3]*m[7]-m[2]*m[6]*m[4] )
                 
/* Subroutine */ int tred2_(int *nm, int *n, double *a, 
	double *d__, double *e, double *z__)
{
    /* System generated locals */
    int a_dim1, a_offset, z_dim1, z_offset, i__1, i__2, i__3;
    double d__1;

    /* Builtin functions */
    double d_sign(double *, double *);

    /* Local variables */
    double f, g, h__;
    int i__, j, k, l;
    double scale, hh;
    int ii, jp1;



/*     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED2, */
/*     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971). */

/*     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX TO A */
/*     SYMMETRIC TRIDIAGONAL MATRIX USING AND ACCUMULATING */
/*     ORTHOGONAL SIMILARITY TRANSFORMATIONS. */

/*     ON INPUT */

/*        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL */
/*          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM */
/*          DIMENSION STATEMENT. */

/*        N IS THE ORDER OF THE MATRIX. */

/*        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE */
/*          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED. */

/*     ON OUTPUT */

/*        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX. */

/*        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL */
/*          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO. */

/*        Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX */
/*          PRODUCED IN THE REDUCTION. */

/*        A AND Z MAY COINCIDE.  IF DISTINCT, A IS UNALTERED. */

/*     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW, */
/*     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY 
*/

/*     THIS VERSION DATED AUGUST 1983. */

/*     ------------------------------------------------------------------ 
*/

    /* Parameter adjustments */
    z_dim1 = *nm;
    z_offset = z_dim1 + 1;
    z__ -= z_offset;
    --e;
    --d__;
    a_dim1 = *nm;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {

	i__2 = *n;
	for (j = i__; j <= i__2; ++j) {
/* L80: */
	    z__[j + i__ * z_dim1] = a[j + i__ * a_dim1];
	}

	d__[i__] = a[*n + i__ * a_dim1];
/* L100: */
    }

    if (*n == 1) {
	goto L510;
    }
/*     .......... FOR I=N STEP -1 UNTIL 2 DO -- .......... */
    i__1 = *n;
    for (ii = 2; ii <= i__1; ++ii) {
	i__ = *n + 2 - ii;
	l = i__ - 1;
	h__ = 0.;
	scale = 0.;
	if (l < 2) {
	    goto L130;
	}
/*     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) .......... */
	i__2 = l;
	for (k = 1; k <= i__2; ++k) {
/* L120: */
	    scale += (d__1 = d__[k], fabs(d__1));
	}

	if (scale != 0.) {
	    goto L140;
	}
L130:
	e[i__] = d__[l];

	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    d__[j] = z__[l + j * z_dim1];
	    z__[i__ + j * z_dim1] = 0.;
	    z__[j + i__ * z_dim1] = 0.;
/* L135: */
	}

	goto L290;

L140:
	i__2 = l;
	for (k = 1; k <= i__2; ++k) {
	    d__[k] /= scale;
	    h__ += d__[k] * d__[k];
/* L150: */
	}

	f = d__[l];
	d__1 = sqrt(h__);
	g = -d_sign(&d__1, &f);
	e[i__] = scale * g;
	h__ -= f * g;
	d__[l] = f - g;
/*     .......... FORM A*U .......... */
	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
/* L170: */
	    e[j] = 0.;
	}

	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    f = d__[j];
	    z__[j + i__ * z_dim1] = f;
	    g = e[j] + z__[j + j * z_dim1] * f;
	    jp1 = j + 1;
	    if (l < jp1) {
		goto L220;
	    }

	    i__3 = l;
	    for (k = jp1; k <= i__3; ++k) {
		g += z__[k + j * z_dim1] * d__[k];
		e[k] += z__[k + j * z_dim1] * f;
/* L200: */
	    }

L220:
	    e[j] = g;
/* L240: */
	}
/*     .......... FORM P .......... */
	f = 0.;

	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    e[j] /= h__;
	    f += e[j] * d__[j];
/* L245: */
	}

	hh = f / (h__ + h__);
/*     .......... FORM Q .......... */
	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
/* L250: */
	    e[j] -= hh * d__[j];
	}
/*     .......... FORM REDUCED A .......... */
	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    f = d__[j];
	    g = e[j];

	    i__3 = l;
	    for (k = j; k <= i__3; ++k) {
/* L260: */
		z__[k + j * z_dim1] = z__[k + j * z_dim1] - f * e[k] - g * 
			d__[k];
	    }

	    d__[j] = z__[l + j * z_dim1];
	    z__[i__ + j * z_dim1] = 0.;
/* L280: */
	}

L290:
	d__[i__] = h__;
/* L300: */
    }
/*     .......... ACCUMULATION OF TRANSFORMATION MATRICES .......... */
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	l = i__ - 1;
	z__[*n + l * z_dim1] = z__[l + l * z_dim1];
	z__[l + l * z_dim1] = 1.;
	h__ = d__[i__];
	if (h__ == 0.) {
	    goto L380;
	}

	i__2 = l;
	for (k = 1; k <= i__2; ++k) {
/* L330: */
	    d__[k] = z__[k + i__ * z_dim1] / h__;
	}

	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    g = 0.;

	    i__3 = l;
	    for (k = 1; k <= i__3; ++k) {
/* L340: */
		g += z__[k + i__ * z_dim1] * z__[k + j * z_dim1];
	    }

	    i__3 = l;
	    for (k = 1; k <= i__3; ++k) {
		z__[k + j * z_dim1] -= g * d__[k];
/* L360: */
	    }
	}

L380:
	i__3 = l;
	for (k = 1; k <= i__3; ++k) {
/* L400: */
	    z__[k + i__ * z_dim1] = 0.;
	}

/* L500: */
    }

L510:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__[i__] = z__[*n + i__ * z_dim1];
	z__[*n + i__ * z_dim1] = 0.;
/* L520: */
    }

    z__[*n + *n * z_dim1] = 1.;
    e[1] = 0.;
    return 0;
} /* tred2_ */

int tred1_(int *nm, int *n, double *a, 
	double *d__, double *e, double *e2)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2, i__3;
    double d__1;

    /* Builtin functions */
    double d_sign(double *, double *);

    /* Local variables */
    double f, g, h__;
    int i__, j, k, l;
    double scale;
    int ii, jp1;



/*     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED1, */
/*     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971). */

/*     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX */
/*     TO A SYMMETRIC TRIDIAGONAL MATRIX USING */
/*     ORTHOGONAL SIMILARITY TRANSFORMATIONS. */

/*     ON INPUT */

/*        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL */
/*          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM */
/*          DIMENSION STATEMENT. */

/*        N IS THE ORDER OF THE MATRIX. */

/*        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE */
/*          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED. */

/*     ON OUTPUT */

/*        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS- */
/*          FORMATIONS USED IN THE REDUCTION IN ITS STRICT LOWER */
/*          TRIANGLE.  THE FULL UPPER TRIANGLE OF A IS UNALTERED. */

/*        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX. */

/*        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL */
/*          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO. */

/*        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E. */
/*          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED. */

/*     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW, */
/*     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY 
*/

/*     THIS VERSION DATED AUGUST 1983. */

/*     ------------------------------------------------------------------ 
*/

    /* Parameter adjustments */
    --e2;
    --e;
    --d__;
    a_dim1 = *nm;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__[i__] = a[*n + i__ * a_dim1];
	a[*n + i__ * a_dim1] = a[i__ + i__ * a_dim1];
/* L100: */
    }
/*     .......... FOR I=N STEP -1 UNTIL 1 DO -- .......... */
    i__1 = *n;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = *n + 1 - ii;
	l = i__ - 1;
	h__ = 0.;
	scale = 0.;
	if (l < 1) {
	    goto L130;
	}
/*     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) .......... */
	i__2 = l;
	for (k = 1; k <= i__2; ++k) {
/* L120: */
	    scale += (d__1 = d__[k], fabs(d__1));
	}

	if (scale != 0.) {
	    goto L140;
	}

	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    d__[j] = a[l + j * a_dim1];
	    a[l + j * a_dim1] = a[i__ + j * a_dim1];
	    a[i__ + j * a_dim1] = 0.;
/* L125: */
	}

L130:
	e[i__] = 0.;
	e2[i__] = 0.;
	goto L300;

L140:
	i__2 = l;
	for (k = 1; k <= i__2; ++k) {
	    d__[k] /= scale;
	    h__ += d__[k] * d__[k];
/* L150: */
	}

	e2[i__] = scale * scale * h__;
	f = d__[l];
	d__1 = sqrt(h__);
	g = -d_sign(&d__1, &f);
	e[i__] = scale * g;
	h__ -= f * g;
	d__[l] = f - g;
	if (l == 1) {
	    goto L285;
	}
/*     .......... FORM A*U .......... */
	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
/* L170: */
	    e[j] = 0.;
	}

	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    f = d__[j];
	    g = e[j] + a[j + j * a_dim1] * f;
	    jp1 = j + 1;
	    if (l < jp1) {
		goto L220;
	    }

	    i__3 = l;
	    for (k = jp1; k <= i__3; ++k) {
		g += a[k + j * a_dim1] * d__[k];
		e[k] += a[k + j * a_dim1] * f;
/* L200: */
	    }

L220:
	    e[j] = g;
/* L240: */
	}
/*     .......... FORM P .......... */
	f = 0.;

	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    e[j] /= h__;
	    f += e[j] * d__[j];
/* L245: */
	}

	h__ = f / (h__ + h__);
/*     .......... FORM Q .......... */
	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
/* L250: */
	    e[j] -= h__ * d__[j];
	}
/*     .......... FORM REDUCED A .......... */
	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    f = d__[j];
	    g = e[j];

	    i__3 = l;
	    for (k = j; k <= i__3; ++k) {
/* L260: */
		a[k + j * a_dim1] = a[k + j * a_dim1] - f * e[k] - g * d__[k];
	    }

/* L280: */
	}

L285:
	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    f = d__[j];
	    d__[j] = a[l + j * a_dim1];
	    a[l + j * a_dim1] = a[i__ + j * a_dim1];
	    a[i__ + j * a_dim1] = f * scale;
/* L290: */
	}

L300:
	;
    }

    return 0;
} /* tred1_ */

double epslon_(double *x)
{
    /* System generated locals */
    double ret_val, d__1;

    /* Local variables */
    double a, b, c__, eps;


/*     ESTIMATE UNIT ROUNDOFF IN QUANTITIES OF SIZE X. */


/*     THIS PROGRAM SHOULD FUNCTION PROPERLY ON ALL SYSTEMS */
/*     SATISFYING THE FOLLOWING TWO ASSUMPTIONS, */
/*        1.  THE BASE USED IN REPRESENTING FLOATING POINT */
/*            NUMBERS IS NOT A POWER OF THREE. */
/*        2.  THE QUANTITY  A  IN STATEMENT 10 IS REPRESENTED TO */
/*            THE ACCURACY USED IN FLOATING POINT VARIABLES */
/*            THAT ARE STORED IN MEMORY. */
/*     THE STATEMENT NUMBER 10 AND THE GO TO 10 ARE INTENDED TO */
/*     FORCE OPTIMIZING COMPILERS TO GENERATE CODE SATISFYING */
/*     ASSUMPTION 2. */
/*     UNDER THESE ASSUMPTIONS, IT SHOULD BE TRUE THAT, */
/*            A  IS NOT EXACTLY EQUAL TO FOUR-THIRDS, */
/*            B  HAS A ZERO FOR ITS LAST BIT OR DIGIT, */
/*            C  IS NOT EXACTLY EQUAL TO ONE, */
/*            EPS  MEASURES THE SEPARATION OF 1.0 FROM */
/*                 THE NEXT LARGER FLOATING POINT NUMBER. */
/*     THE DEVELOPERS OF EISPACK WOULD APPRECIATE BEING INFORMED */
/*     ABOUT ANY SYSTEMS WHERE THESE ASSUMPTIONS DO NOT HOLD. */

/*     THIS VERSION DATED 4/6/83. */

    a = 1.3333333333333333;
L10:
    b = a - 1.;
    c__ = b + b + b;
    eps = (d__1 = c__ - 1., fabs(d__1));
    if (eps == 0.) {
	goto L10;
    }
    ret_val = eps * fabs(*x);
    return ret_val;
} /* epslon_ */


static double c_b11 = 1.;

/* Subroutine */ int tqlrat_(int *n, double *d__, double *e2, 
	int *ierr)
{
    /* System generated locals */
    int i__1, i__2;
    double d__1, d__2;

    /* Builtin functions */
    double d_sign(double *, double *);

    /* Local variables */
    double b=0.0, c__=0.0, f, g, h__;
    int i__, j, l, m;
    double p, r__, s, t;
    int l1, ii;
    extern double pythag_(double *, double *), epslon_(double 
	    *);
    int mml;



/*     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT, */
/*     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH. */

/*     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC */
/*     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD. */

/*     ON INPUT */

/*        N IS THE ORDER OF THE MATRIX. */

/*        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX. */

/*        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE */
/*          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY. 
*/

/*      ON OUTPUT */

/*        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN */
/*          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND */
/*          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE */
/*          THE SMALLEST EIGENVALUES. */

/*        E2 HAS BEEN DESTROYED. */

/*        IERR IS SET TO */
/*          ZERO       FOR NORMAL RETURN, */
/*          J          IF THE J-TH EIGENVALUE HAS NOT BEEN */
/*                     DETERMINED AFTER 30 ITERATIONS. */

/*     CALLS PYTHAG FOR  DSQRT(A*A + B*B) . */

/*     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW, */
/*     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY 
*/

/*     THIS VERSION DATED AUGUST 1983. */

/*     ------------------------------------------------------------------ 
*/

    /* Parameter adjustments */
    --e2;
    --d__;

    /* Function Body */
    *ierr = 0;
    if (*n == 1) {
	goto L1001;
    }

    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L100: */
	e2[i__ - 1] = e2[i__];
    }

    f = 0.;
    t = 0.;
    e2[*n] = 0.;

    i__1 = *n;
    for (l = 1; l <= i__1; ++l) {
	j = 0;
	h__ = (d__1 = d__[l], fabs(d__1)) + sqrt(e2[l]);
	if (t > h__) {
	    goto L105;
	}
	t = h__;
	b = epslon_(&t);
	c__ = b * b;
/*     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ........
.. */
L105:
	i__2 = *n;
	for (m = l; m <= i__2; ++m) {
	    if (e2[m] <= c__) {
		goto L120;
	    }
/*     .......... E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT */
/*                THROUGH THE BOTTOM OF THE LOOP .......... */
/* L110: */
	}

L120:
	if (m == l) {
	    goto L210;
	}
L130:
	if (j == 30) {
	    goto L1000;
	}
	++j;
/*     .......... FORM SHIFT .......... */
	l1 = l + 1;
	s = sqrt(e2[l]);
	g = d__[l];
	p = (d__[l1] - g) / (s * 2.);
	r__ = pythag_(&p, &c_b11);
	d__[l] = s / (p + d_sign(&r__, &p));
	h__ = g - d__[l];

	i__2 = *n;
	for (i__ = l1; i__ <= i__2; ++i__) {
/* L140: */
	    d__[i__] -= h__;
	}

	f += h__;
/*     .......... RATIONAL QL TRANSFORMATION .......... */
	g = d__[m];
	if (g == 0.) {
	    g = b;
	}
	h__ = g;
	s = 0.;
	mml = m - l;
/*     .......... FOR I=M-1 STEP -1 UNTIL L DO -- .......... */
	i__2 = mml;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = m - ii;
	    p = g * h__;
	    r__ = p + e2[i__];
	    e2[i__ + 1] = s * r__;
	    s = e2[i__] / r__;
	    d__[i__ + 1] = h__ + s * (h__ + d__[i__]);
	    g = d__[i__] - e2[i__] / g;
	    if (g == 0.) {
		g = b;
	    }
	    h__ = g * p / r__;
/* L200: */
	}

	e2[l] = s * g;
	d__[l] = h__;
/*     .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ........
.. */
	if (h__ == 0.) {
	    goto L210;
	}
	if ((d__1 = e2[l], fabs(d__1)) <= (d__2 = c__ / h__, fabs(d__2))) {
	    goto L210;
	}
	e2[l] = h__ * e2[l];
	if (e2[l] != 0.) {
	    goto L130;
	}
L210:
	p = d__[l] + f;
/*     .......... ORDER EIGENVALUES .......... */
	if (l == 1) {
	    goto L250;
	}
/*     .......... FOR I=L STEP -1 UNTIL 2 DO -- .......... */
	i__2 = l;
	for (ii = 2; ii <= i__2; ++ii) {
	    i__ = l + 2 - ii;
	    if (p >= d__[i__ - 1]) {
		goto L270;
	    }
	    d__[i__] = d__[i__ - 1];
/* L230: */
	}

L250:
	i__ = 1;
L270:
	d__[i__] = p;
/* L290: */
    }

    goto L1001;
/*     .......... SET ERROR -- NO CONVERGENCE TO AN */
/*                EIGENVALUE AFTER 30 ITERATIONS .......... */
L1000:
    *ierr = l;
L1001:
    return 0;
} /* tqlrat_ */

double pythag_(double *a, double *b)
{
    /* System generated locals */
    double ret_val, d__1, d__2, d__3;

    /* Local variables */
    double p, r__, s, t, u;


/*     FINDS DSQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW */

/* Computing MAX */
    d__1 = fabs(*a), d__2 = fabs(*b);
    p = MAX(d__1,d__2);
    if (p == 0.) {
	goto L20;
    }
/* Computing MIN */
    d__2 = fabs(*a), d__3 = fabs(*b);
/* Computing 2nd power */
    d__1 = MIN(d__2,d__3) / p;
    r__ = d__1 * d__1;
L10:
    t = r__ + 4.;
    if (t == 4.) {
	goto L20;
    }
    s = r__ / t;
    u = s * 2. + 1.;
    p = u * p;
/* Computing 2nd power */
    d__1 = s / u;
    r__ = d__1 * d__1 * r__;
    goto L10;
L20:
    ret_val = p;
    return ret_val;
} /* pythag_ */

double d_sign(double *a, double *b) {
	double x = (*a >= 0 ? *a : - *a);
	return( *b >= 0 ? x : -x);
}

static double c_b10 = 1.;

/* Subroutine */ int tql2_(int *nm, int *n, double *d__, 
	double *e, double *z__, int *ierr)
{
    /* System generated locals */
    int z_dim1, z_offset, i__1, i__2, i__3;
    double d__1, d__2;

    /* Builtin functions */
    double d_sign(double *, double *);

    /* Local variables */
    double c__, f, g, h__;
    int i__, j, k, l, m;
    double p=0.0, r__=0.0, s=0.0, c2=0.0, c3=0.0;
    int l1, l2;
    double s2=0.0;
    int ii;
    extern double pythag_(double *, double *);
    double dl1, el1;
    int mml;
    double tst1, tst2;



/*     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2, */
/*     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND */
/*     WILKINSON. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971). */

/*     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS */
/*     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD. */
/*     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO */
/*     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS */
/*     FULL MATRIX TO TRIDIAGONAL FORM. */

/*     ON INPUT */

/*        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL */
/*          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM */
/*          DIMENSION STATEMENT. */

/*        N IS THE ORDER OF THE MATRIX. */

/*        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX. */

/*        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX */
/*          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY. */

/*        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE */
/*          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS */
/*          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN */
/*          THE IDENTITY MATRIX. */

/*      ON OUTPUT */

/*        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN */
/*          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT */
/*          UNORDERED FOR INDICES 1,2,...,IERR-1. */

/*        E HAS BEEN DESTROYED. */

/*        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC */
/*          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE, */
/*          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED */
/*          EIGENVALUES. */

/*        IERR IS SET TO */
/*          ZERO       FOR NORMAL RETURN, */
/*          J          IF THE J-TH EIGENVALUE HAS NOT BEEN */
/*                     DETERMINED AFTER 30 ITERATIONS. */

/*     CALLS PYTHAG FOR  DSQRT(A*A + B*B) . */

/*     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW, */
/*     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY 
*/

/*     THIS VERSION DATED AUGUST 1983. */

/*     ------------------------------------------------------------------ 
*/

    /* Parameter adjustments */
    z_dim1 = *nm;
    z_offset = z_dim1 + 1;
    z__ -= z_offset;
    --e;
    --d__;

    /* Function Body */
    *ierr = 0;
    if (*n == 1) {
	goto L1001;
    }

    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L100: */
	e[i__ - 1] = e[i__];
    }

    f = 0.;
    tst1 = 0.;
    e[*n] = 0.;

    i__1 = *n;
    for (l = 1; l <= i__1; ++l) {
	j = 0;
	h__ = (d__1 = d__[l], fabs(d__1)) + (d__2 = e[l], fabs(d__2));
	if (tst1 < h__) {
	    tst1 = h__;
	}
/*     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT .......... */
	i__2 = *n;
	for (m = l; m <= i__2; ++m) {
	    tst2 = tst1 + (d__1 = e[m], fabs(d__1));
	    if (tst2 == tst1) {
		goto L120;
	    }
/*     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT */
/*                THROUGH THE BOTTOM OF THE LOOP .......... */
/* L110: */
	}

L120:
	if (m == l) {
	    goto L220;
	}
L130:
	if (j == 30) {
	    goto L1000;
	}
	++j;
/*     .......... FORM SHIFT .......... */
	l1 = l + 1;
	l2 = l1 + 1;
	g = d__[l];
	p = (d__[l1] - g) / (e[l] * 2.);
	r__ = pythag_(&p, &c_b10);
	d__[l] = e[l] / (p + d_sign(&r__, &p));
	d__[l1] = e[l] * (p + d_sign(&r__, &p));
	dl1 = d__[l1];
	h__ = g - d__[l];
	if (l2 > *n) {
	    goto L145;
	}

	i__2 = *n;
	for (i__ = l2; i__ <= i__2; ++i__) {
/* L140: */
	    d__[i__] -= h__;
	}

L145:
	f += h__;
/*     .......... QL TRANSFORMATION .......... */
	p = d__[m];
	c__ = 1.;
	c2 = c__;
	el1 = e[l1];
	s = 0.;
	mml = m - l;
/*     .......... FOR I=M-1 STEP -1 UNTIL L DO -- .......... */
	i__2 = mml;
	for (ii = 1; ii <= i__2; ++ii) {
	    c3 = c2;
	    c2 = c__;
	    s2 = s;
	    i__ = m - ii;
	    g = c__ * e[i__];
	    h__ = c__ * p;
	    r__ = pythag_(&p, &e[i__]);
	    e[i__ + 1] = s * r__;
	    s = e[i__] / r__;
	    c__ = p / r__;
	    p = c__ * d__[i__] - s * g;
	    d__[i__ + 1] = h__ + s * (c__ * g + s * d__[i__]);
/*     .......... FORM VECTOR .......... */
	    i__3 = *n;
	    for (k = 1; k <= i__3; ++k) {
		h__ = z__[k + (i__ + 1) * z_dim1];
		z__[k + (i__ + 1) * z_dim1] = s * z__[k + i__ * z_dim1] + c__ 
			* h__;
		z__[k + i__ * z_dim1] = c__ * z__[k + i__ * z_dim1] - s * h__;
/* L180: */
	    }

/* L200: */
	}

	p = -s * s2 * c3 * el1 * e[l] / dl1;
	e[l] = s * p;
	d__[l] = c__ * p;
	tst2 = tst1 + (d__1 = e[l], fabs(d__1));
	if (tst2 > tst1) {
	    goto L130;
	}
L220:
	d__[l] += f;
/* L240: */
    }
/*     .......... ORDER EIGENVALUES AND EIGENVECTORS .......... */
    i__1 = *n;
    for (ii = 2; ii <= i__1; ++ii) {
	i__ = ii - 1;
	k = i__;
	p = d__[i__];

	i__2 = *n;
	for (j = ii; j <= i__2; ++j) {
	    if (d__[j] >= p) {
		goto L260;
	    }
	    k = j;
	    p = d__[j];
L260:
	    ;
	}

	if (k == i__) {
	    goto L300;
	}
	d__[k] = d__[i__];
	d__[i__] = p;

	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    p = z__[j + i__ * z_dim1];
	    z__[j + i__ * z_dim1] = z__[j + k * z_dim1];
	    z__[j + k * z_dim1] = p;
/* L280: */
	}

L300:
	;
    }

    goto L1001;
/*     .......... SET ERROR -- NO CONVERGENCE TO AN */
/*                EIGENVALUE AFTER 30 ITERATIONS .......... */
L1000:
    *ierr = l;
L1001:
    return 0;
} /* tql2_ */

int rs_(int *nm, int *n, double *a, double *
	w, int *matz, double *z__, double *fv1, double *fv2, 
	int *ierr) {
    /* System generated locals */
    int a_dim1, a_offset, z_dim1, z_offset;

    /* Local variables */
    extern /* Subroutine */ int tred1_(int *, int *, double *, 
	    double *, double *, double *), tred2_(int *, 
	    int *, double *, double *, double *, double *)
	    , tqlrat_(int *, double *, double *, int *), 
	    tql2_(int *, int *, double *, double *, 
	    double *, int *);



/*     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF */
/*     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK) */
/*     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED) */
/*     OF A REAL SYMMETRIC MATRIX. */

/*     ON INPUT */

/*        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL */
/*        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM */
/*        DIMENSION STATEMENT. */

/*        N  IS THE ORDER OF THE MATRIX  A. */

/*        A  CONTAINS THE REAL SYMMETRIC MATRIX. */

/*        MATZ  IS AN int VARIABLE SET EQUAL TO ZERO IF */
/*        ONLY EIGENVALUES ARE DESIRED.  OTHERWISE IT IS SET TO */
/*        ANY NON-ZERO int FOR BOTH EIGENVALUES AND EIGENVECTORS. */

/*     ON OUTPUT */

/*        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER. */

/*        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO. */

/*        IERR  IS AN int OUTPUT VARIABLE SET EQUAL TO AN ERROR */
/*           COMPLETION CODE DESCRIBED IN THE DOCUMENTATION FOR TQLRAT */
/*           AND TQL2.  THE NORMAL COMPLETION CODE IS ZERO. */

/*        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS. */

/*     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW, */
/*     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY 
*/

/*     THIS VERSION DATED AUGUST 1983. */

/*     ------------------------------------------------------------------ 
*/

    /* Parameter adjustments */
    --fv2;
    --fv1;
    z_dim1 = *nm;
    z_offset = z_dim1 + 1;
    z__ -= z_offset;
    --w;
    a_dim1 = *nm;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    if (*n <= *nm) {
	goto L10;
    }
    *ierr = *n * 10;
    goto L50;

L10:
    if (*matz != 0) {
	goto L20;
    }
/*     .......... FIND EIGENVALUES ONLY .......... */
    tred1_(nm, n, &a[a_offset], &w[1], &fv1[1], &fv2[1]);
    tqlrat_(n, &w[1], &fv2[1], ierr);
    goto L50;
/*     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS .......... */
L20:
    tred2_(nm, n, &a[a_offset], &w[1], &fv1[1], &z__[z_offset]);
    tql2_(nm, n, &w[1], &fv1[1], &z__[z_offset], ierr);
L50:
    return 0;
} /* rs_ */

void symeig_3( double *a , double *e , int dovec )
{
   double aa,bb,cc,dd,ee,ff ;
   double a1,a2,a3 , qq,rr, qs,th , lam1,lam2,lam3 ;
   double aba,abb,abc,abd,abe,abf , ann ;
   double d12,d13,d23 ;
   double u1,u2,u3 , v1,v2,v3 , w1,w2,w3 , t1,t2,t3 , tn ;
   double anni ;

   if( a == NULL || e == NULL ) return ;

   /*----- unload matrix into local variables -----*/

   aa = a[0] ; bb = a[1] ; cc = a[2] ;  /* matrix is [ aa bb cc ]  */
   dd = a[4] ; ee = a[5] ; ff = a[8] ;  /*           [ bb dd ee ]  */
                                        /*           [ cc ee ff ]  */
   aba = fabs(aa) ; abb = fabs(bb) ; abc = fabs(cc) ;
   abd = fabs(dd) ; abe = fabs(ee) ; abf = fabs(ff) ;
   ann = aba+abb+abc+abd+abe+abf   ;                 /* matrix 'norm' */

   if( ann == 0.0 ){             /* matrix is all zero! */
     e[0] = e[1] = e[2] = 0.0 ;
     if( dovec ){
       a[0] = a[4] = a[8] = 1.0 ;
       a[1] = a[2] = a[3] = a[5] = a[6] = a[7] = 0.0 ;
     }
     return ;
   }

   /*----- check for matrix that is essentially diagonal -----*/

   if( abb+abc+abe == 0.0 ||
       ( EPS*aba > (abb+abc) && EPS*abd > (abb+abe) && EPS*abf > (abc+abe) ) ){

     lam1 = aa ; lam2 = dd ; lam3 = ff ;

     if( dovec ){
       a[0] = a[4] = a[8] = 1.0 ;
       a[1] = a[2] = a[3] = a[5] = a[6] = a[7] = 0.0 ;

       if( lam1 > lam2 ){ SWAP(lam1,lam2) ; CSWAP(0,3) ; }
       if( lam1 > lam3 ){ SWAP(lam1,lam3) ; CSWAP(0,6) ; }
       if( lam2 > lam3 ){ SWAP(lam2,lam3) ; CSWAP(3,6) ; }
       if( DET3(a) < 0.0 ){ a[6] = -a[6]; a[7] = -a[7]; a[8] = -a[8]; }
     } else {
       if( lam1 > lam2 )  SWAP(lam1,lam2) ;
       if( lam1 > lam3 )  SWAP(lam1,lam3) ;
       if( lam2 > lam3 )  SWAP(lam2,lam3) ;
     }
     e[0] = lam1 ; e[1] = lam2 ; e[2] = lam3 ;
     return ;
   }

   /*-- Scale matrix so abs sum is 1; unscale e[i] on output [26 Oct 2005] --*/

   anni = 1.0 / ann ;                      /* ann != 0, from above */
   aa *= anni ; bb *= anni ; cc *= anni ;
   dd *= anni ; ee *= anni ; ff *= anni ;

   /*----- not diagonal ==> must solve cubic polynomial for eigenvalues -----*/
   /*      the cubic polynomial is x**3 + a1*x**2 + a2*x + a3 = 0            */

   a1 = -(aa+dd+ff) ;
   a2 =  (aa*ff+aa*dd+dd*ff - bb*bb-cc*cc-ee*ee) ;
   a3 =  ( aa*(ee*ee-dd*ff) + bb*(bb*ff-cc*ee) + cc*(cc*dd-bb*ee) ) ;

   /*-- Rewrite classical formula for qq as a sum of squares [26 Oct 2005] --*/
#if 0
   qq = (a1*a1 - 3.0*a2) / 9.0 ;
#else
   qq = (  0.5 * ( SQR(dd-aa) + SQR(ff-aa) + SQR(ff-dd) )
         + 3.0 * ( bb*bb      + cc*cc      + ee*ee      ) ) / 9.0 ;
#endif
   rr = (2.0*a1*a1*a1 - 9.0*a1*a2 + 27.0*a3) / 54.0 ;

   if( qq <= 0.0 ){       /*** This should never happen!!! ***/
     static int nerr=0 ;
     {
     if( ++nerr < 4 )
       printfx("** ERROR in symeig_3: discrim=%g numerator=%g\n",qq,rr) ;
     }
     qs = qq = rr = 0.0 ;
   } else {
     qs = sqrt(qq) ; rr = rr / (qs*qq) ;
     if( rr < -1.0 ) rr = -1.0 ; else if( rr > 1.0 ) rr = 1.0 ;
   }
   th = acos(rr) ;

   lam1 = -2.0 * qs * cos(  th        /3.0 ) - a1 / 3.0 ;
   lam2 = -2.0 * qs * cos( (th+2.0*PI)/3.0 ) - a1 / 3.0 ;
   lam3 = -2.0 * qs * cos( (th+4.0*PI)/3.0 ) - a1 / 3.0 ;

   /*-- if not doing eigenvectors, just sort the eigenvalues to be done --*/

   if( !dovec ){
     if( lam1 > lam2 ) SWAP(lam1,lam2) ;
     if( lam1 > lam3 ) SWAP(lam1,lam3) ;
     if( lam2 > lam3 ) SWAP(lam2,lam3) ;
     e[0] = ann*lam1 ; e[1] = ann*lam2 ; e[2] = ann*lam3 ;  /* re-scale */
     return ;
   }

   /*-- are doing eigenvectors; must do double root as a special case --*/

#undef  CROSS  /* cross product (x1,x2,x3) X (y1,y2,y3) -> (z1,z2,z3) */
#define CROSS(x1,x2,x3,y1,y2,y3,z1,z2,z3) \
 ( (z1)=(x2)*(y3)-(x3)*(y2), (z2)=(x3)*(y1)-(x1)*(y3), (z3)=(x1)*(y2)-(x2)*(y1) )

   d12 = fabs(lam1-lam2) ; d13 = fabs(lam1-lam3) ; d23 = fabs(lam2-lam3) ;
   rr  = MIN(d12,d13)    ; rr  = MIN(rr,d23)     ;

   if( rr > EPS*ann ){  /*---- not a double root ----*/

     if( lam1 > lam2 )  SWAP(lam1,lam2) ;  /* start by sorting eigenvalues */
     if( lam1 > lam3 )  SWAP(lam1,lam3) ;
     if( lam2 > lam3 )  SWAP(lam2,lam3) ;
     e[0] = ann*lam1 ; e[1] = ann*lam2 ; e[2] = ann*lam3 ;  /* re-scale */

     /* find eigenvector for lam1 by computing Ay-lam1*y for
        vectors y1=[1,0,0], y2=[0,1,0], and y3=[0,0,1]; the eigenvector
        is orthogonal to all of these, so use the cross product to get it */

     u1 = aa-lam1 ; u2 = bb      ; u3 = cc ;   /* A*y1 - lam1*y1 */
     v1 = bb      ; v2 = dd-lam1 ; v3 = ee ;   /* A*y2 - lam1*y2 */
     CROSS(u1,u2,u3 , v1,v2,v3 , t1,t2,t3 ) ;     tn = sqrt(t1*t1+t2*t2+t3*t3) ;
     if( tn < EPSQ*ann ){                      /* u and v were parallel? */
       w1 = cc ; w2 = ee ; w3 = ff-lam1 ;      /* A*y3 - lam1*y3 */
       CROSS(u1,u2,u3 , w1,w2,w3 , t1,t2,t3 ) ;   tn = sqrt(t1*t1+t2*t2+t3*t3) ;
       if( tn < EPSQ*ann ){                    /* u and w were parallel? */
         CROSS(v1,v2,v3 , w1,w2,w3 , t1,t2,t3 ) ; tn = sqrt(t1*t1+t2*t2+t3*t3) ;
       }
     }
     a[0] = t1/tn ; a[1] = t2/tn ; a[2] = t3/tn ;  /* normalize */

     /* do same for lam2 */

     u1 = aa-lam2 ; u2 = bb      ; u3 = cc ;
     v1 = bb      ; v2 = dd-lam2 ; v3 = ee ;
     CROSS(u1,u2,u3 , v1,v2,v3 , t1,t2,t3 ) ;     tn = sqrt(t1*t1+t2*t2+t3*t3) ;
     if( tn < EPSQ*ann ){
       w1 = cc ; w2 = ee ; w3 = ff-lam2 ;
       CROSS(u1,u2,u3 , w1,w2,w3 , t1,t2,t3 ) ;   tn = sqrt(t1*t1+t2*t2+t3*t3) ;
       if( tn < EPSQ*ann ){
         CROSS(v1,v2,v3 , w1,w2,w3 , t1,t2,t3 ) ; tn = sqrt(t1*t1+t2*t2+t3*t3) ;
       }
     }
     a[3] = t1/tn ; a[4] = t2/tn ; a[5] = t3/tn ;  /* normalize */

     /* orthgonality of eigenvectors ==> can get last one by cross product */

#if 1
     CROSS( a[0],a[1],a[2] , a[3],a[4],a[5] , a[6],a[7],a[8] ) ;
#else
     u1 = aa-lam3 ; u2 = bb      ; u3 = cc ;
     v1 = bb      ; v2 = dd-lam3 ; v3 = ee ;
     CROSS(u1,u2,u3 , v1,v2,v3 , t1,t2,t3 ) ;     tn = sqrt(t1*t1+t2*t2+t3*t3) ;
     if( tn < EPSQ*ann ){
       w1 = cc ; w2 = ee ; w3 = ff-lam3 ;
       CROSS(u1,u2,u3 , w1,w2,w3 , t1,t2,t3 ) ;   tn = sqrt(t1*t1+t2*t2+t3*t3) ;
       if( tn < EPSQ*ann ){
         CROSS(v1,v2,v3 , w1,w2,w3 , t1,t2,t3 ) ; tn = sqrt(t1*t1+t2*t2+t3*t3) ;
       }
     }
     a[6] = t1/tn ; a[7] = t2/tn ; a[8] = t3/tn ;  /* normalize */
#endif

   } else { /*---- if here, we have a double root ----*/

     /* make sure that we have lam1=lam2 and lam3 is the outlier */

          if( d13 < d12 && d13 < d23 ) SWAP(lam2,lam3) ;
     else if( d23 < d12 && d23 < d13 ) SWAP(lam1,lam3) ;
     lam1 = lam2 = 0.5*(lam1+lam2) ;

     /* compute eigenvector for lam3 using method as above */

     u1 = aa-lam3 ; u2 = bb      ; u3 = cc ;
     v1 = bb      ; v2 = dd-lam3 ; v3 = ee ;
     CROSS(u1,u2,u3 , v1,v2,v3 , t1,t2,t3 ) ;     tn = sqrt(t1*t1+t2*t2+t3*t3) ;
     if( tn < EPSQ*ann ){
       w1 = cc ; w2 = ee ; w3 = ff-lam3 ;
       CROSS(u1,u2,u3 , w1,w2,w3 , t1,t2,t3 ) ;   tn = sqrt(t1*t1+t2*t2+t3*t3) ;
       if( tn < EPSQ*ann ){
         CROSS(v1,v2,v3 , w1,w2,w3 , t1,t2,t3 ) ; tn = sqrt(t1*t1+t2*t2+t3*t3) ;
       }
     }
     w1 = a[6] = t1/tn ; w2 = a[7] = t2/tn ; w3 = a[8] = t3/tn ;

     /* find a vector orthogonal to it */

     CROSS(w1,w2,w3 , 1,0,0 , t1,t2,t3 ) ; tn = sqrt(t1*t1+t2*t2+t3*t3) ;
     if( tn < EPSQ ){
       CROSS(w1,w2,w3 , 0,1,0 , t1,t2,t3 ) ; tn = sqrt(t1*t1+t2*t2+t3*t3) ;
       if( tn < EPSQ ){
         CROSS(w1,w2,w3 , 0,0,1 , t1,t2,t3 ) ; tn = sqrt(t1*t1+t2*t2+t3*t3) ;
       }
     }
     a[0] = t1/tn ; a[1] = t2/tn ; a[2] = t3/tn ;

     /* and the final vector is the cross product of these two */

     CROSS( w1,w2,w3 , a[0],a[1],a[2] , t1,t2,t3 ) ; tn = sqrt(t1*t1+t2*t2+t3*t3) ;
     a[3] = t1/tn ; a[4] = t2/tn ; a[5] = t3/tn ;

     /* sort results (we know lam1==lam2) */

     if( lam1 > lam3 ){
       SWAP(lam1,lam3) ; CSWAP(0,6) ;
       if( DET3(a) < 0.0 ){ a[6] = -a[6]; a[7] = -a[7]; a[8] = -a[8]; }
     }

     e[0] = ann*lam1 ; e[1] = ann*lam2 ; e[2] = ann*lam3 ;  /* re-scale */
   }

   return ;
}

/*---------------------------------------------------------------------------*/
/*! 2x2 symmetric eigenvalue/vector problem, like symeig_3() above. */

void symeig_2( double *a , double *e , int dovec )
{
   double sxx,sxy,syy , lam1,lam2 , ss,tt , x,y ;

   if( a == NULL || e == NULL ) return ;

   /*----- unload matrix into local variables -----*/

   sxx = a[0] ; sxy = a[1] ; syy = a[3] ;

   ss = fabs(sxx) ; tt = fabs(syy) ; if( ss > tt ) ss = tt ;

   if( fabs(sxy) < EPS*ss ){   /*--- essentially a diagonal matrix ---*/
     if( sxx <= syy ){
       lam1 = sxx ; lam2 = syy ;
       if( dovec ){ a[0]=a[3]=1.0; a[1]=a[2]=0.0; }
     } else {
       lam1 = syy ; lam2 = sxx ;
       if( dovec ){ a[0]=a[3]=1.0 ; a[1]=a[2]=1.0; }
     }
     e[0] = lam1 ; e[1] = lam2 ;
     return ;
   }

   /*--- non-diagonal matrix ==> solve quadratic equation for eigenvalues ---*/

   ss = sqrt( (sxx-syy)*(sxx-syy) + 4.0*sxy*sxy ) ;   /* positive */
   lam1 = 0.5 * ( sxx + syy - ss ) ;                  /* smaller */
   lam2 = 0.5 * ( sxx + syy + ss ) ;                  /* larger */

   if( dovec ){
     x = 2.0*sxy ; y = syy - sxx - ss ; tt = sqrt(x*x+y*y) ;
     a[0] = x/tt ; a[1] = y/tt ;

     y = syy - sxx + ss ; tt = sqrt(x*x+y*y) ;
     a[2] = x/tt ; a[3] = y/tt ;
   }
   e[0] = lam1 ; e[1] = lam2 ;
   return ;
}

static int forbid_23 = 0 ;
void symeig_double( int n , double *a , double *e )
{
   int nm , matz , ierr ;
   double *fv1 , *fv2 ;

   if( a == NULL || e == NULL || n < 1 ) return ;

   /* special cases of small n (much faster than EISPACK) */

   if( n == 1 ){
     e[0] = a[0] ; a[0] = 1.0 ; return ;  /* degenerate case */
   } else if( !forbid_23 ){
     if( n == 2 ){ symeig_2( a , e , 1 ) ; return ; }
     if( n == 3 ){ symeig_3( a , e , 1 ) ; return ; }
   }

   /*-- default code: n > 3 --*/

   fv1 = (double *) malloc(sizeof(double)*n) ;  /* workspaces */
   fv2 = (double *) malloc(sizeof(double)*n) ;

   nm = n ; matz = 1 ; ierr = 0 ;

   rs_( &nm , &nm , a , e , &matz , a , fv1 , fv2 , &ierr ) ;

   free((void *)fv1) ; free((void *)fv2) ;
   return ;
}

#else //if old method (ported from fortran) else newer C method. Older method faster...

//http://barnesc.blogspot.com/2007/02/eigenvectors-of-3x3-symmetric-matrix.html
// Symmetric Householder reduction to tridiagonal form.
// Matlab
//  A = [1 1/2 1/3; 1/2 1 2/3; 1/3 2/3 1]
//  e = eig(A);
//  e' = [0.3020    0.6855    2.0124]
#define n 3

static double hypot2(double x, double y) {
  return sqrt(x*x+y*y);
}



static void tred2(double V[n][n], double d[n], double e[n]) {

//  This is derived from the Algol procedures tred2 by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.

  for (int j = 0; j < n; j++) {
    d[j] = V[n-1][j];
  }

  // Householder reduction to tridiagonal form.

  for (int i = n-1; i > 0; i--) {

    // Scale to avoid under/overflow.

    double scale = 0.0;
    double h = 0.0;
    for (int k = 0; k < i; k++) {
      scale = scale + fabs(d[k]);
    }
    if (scale == 0.0) {
      e[i] = d[i-1];
      for (int j = 0; j < i; j++) {
        d[j] = V[i-1][j];
        V[i][j] = 0.0;
        V[j][i] = 0.0;
      }
    } else {

      // Generate Householder vector.

      for (int k = 0; k < i; k++) {
        d[k] /= scale;
        h += d[k] * d[k];
      }
      double f = d[i-1];
      double g = sqrt(h);
      if (f > 0) {
        g = -g;
      }
      e[i] = scale * g;
      h = h - f * g;
      d[i-1] = f - g;
      for (int j = 0; j < i; j++) {
        e[j] = 0.0;
      }

      // Apply similarity transformation to remaining columns.

      for (int j = 0; j < i; j++) {
        f = d[j];
        V[j][i] = f;
        g = e[j] + V[j][j] * f;
        for (int k = j+1; k <= i-1; k++) {
          g += V[k][j] * d[k];
          e[k] += V[k][j] * f;
        }
        e[j] = g;
      }
      f = 0.0;
      for (int j = 0; j < i; j++) {
        e[j] /= h;
        f += e[j] * d[j];
      }
      double hh = f / (h + h);
      for (int j = 0; j < i; j++) {
        e[j] -= hh * d[j];
      }
      for (int j = 0; j < i; j++) {
        f = d[j];
        g = e[j];
        for (int k = j; k <= i-1; k++) {
          V[k][j] -= (f * e[k] + g * d[k]);
        }
        d[j] = V[i-1][j];
        V[i][j] = 0.0;
      }
    }
    d[i] = h;
  }

  // Accumulate transformations.

  for (int i = 0; i < n-1; i++) {
    V[n-1][i] = V[i][i];
    V[i][i] = 1.0;
    double h = d[i+1];
    if (h != 0.0) {
      for (int k = 0; k <= i; k++) {
        d[k] = V[k][i+1] / h;
      }
      for (int j = 0; j <= i; j++) {
        double g = 0.0;
        for (int k = 0; k <= i; k++) {
          g += V[k][i+1] * V[k][j];
        }
        for (int k = 0; k <= i; k++) {
          V[k][j] -= g * d[k];
        }
      }
    }
    for (int k = 0; k <= i; k++) {
      V[k][i+1] = 0.0;
    }
  }
  for (int j = 0; j < n; j++) {
    d[j] = V[n-1][j];
    V[n-1][j] = 0.0;
  }
  V[n-1][n-1] = 1.0;
  e[0] = 0.0;
} 

// Symmetric tridiagonal QL algorithm.

static void tql2(double V[n][n], double d[n], double e[n]) {

//  This is derived from the Algol procedures tql2, by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.

  for (int i = 1; i < n; i++) {
    e[i-1] = e[i];
  }
  e[n-1] = 0.0;

  double f = 0.0;
  double tst1 = 0.0;
  double eps = pow(2.0,-52.0);
  for (int l = 0; l < n; l++) {

    // Find small subdiagonal element

    tst1 = MAX(tst1,fabs(d[l]) + fabs(e[l]));
    int m = l;
    while (m < n) {
      if (fabs(e[m]) <= eps*tst1) {
        break;
      }
      m++;
    }

    // If m == l, d[l] is an eigenvalue,
    // otherwise, iterate.

    if (m > l) {
      int iter = 0;
      do {
        iter = iter + 1;  // (Could check iteration count here.)

        // Compute implicit shift

        double g = d[l];
        double p = (d[l+1] - g) / (2.0 * e[l]);
        double r = hypot2(p,1.0);
        if (p < 0) {
          r = -r;
        }
        d[l] = e[l] / (p + r);
        d[l+1] = e[l] * (p + r);
        double dl1 = d[l+1];
        double h = g - d[l];
        for (int i = l+2; i < n; i++) {
          d[i] -= h;
        }
        f = f + h;

        // Implicit QL transformation.

        p = d[m];
        double c = 1.0;
        double c2 = c;
        double c3 = c;
        double el1 = e[l+1];
        double s = 0.0;
        double s2 = 0.0;
        for (int i = m-1; i >= l; i--) {
          c3 = c2;
          c2 = c;
          s2 = s;
          g = c * e[i];
          h = c * p;
          r = hypot2(p,e[i]);
          e[i+1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i+1] = h + s * (c * g + s * d[i]);

          // Accumulate transformation.

          for (int k = 0; k < n; k++) {
            h = V[k][i+1];
            V[k][i+1] = s * V[k][i] + c * h;
            V[k][i] = c * V[k][i] - s * h;
          }
        }
        p = -s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;

        // Check for convergence.

      } while (fabs(e[l]) > eps*tst1);
    }
    d[l] = d[l] + f;
    e[l] = 0.0;
  }
  
  // Sort eigenvalues and corresponding vectors.

  for (int i = 0; i < n-1; i++) {
    int k = i;
    double p = d[i];
    for (int j = i+1; j < n; j++) {
      if (d[j] < p) {
        k = j;
        p = d[j];
      }
    }
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (int j = 0; j < n; j++) {
        p = V[j][i];
        V[j][i] = V[j][k];
        V[j][k] = p;
      }
    }
  }
}
#undef n
void symeig_double( int n , double *a , double *e ) {
	//if n != 3 generate error
	double A[3][3];
	double V[3][3];
	int k = 0;
  	for (int i = 0; i < 3; i++) {
    	for (int j = 0; j < 3; j++) {
      		A[i][j] = a[k];
      		V[i][j] = a[k];
      		k++;
    	}
    }
    double d[3];
    tred2(V, d, e);
  	tql2(V, d, e);
    //eigen_decomposition(A, V, d);
    e[0] = d[0];
    e[1] = d[1];
    e[2] = d[2];	
}
#endif //#ifdef old (fortran->C) else (readable C) 

/*********************** 3dDTeig.c **********************************************/
/* Author: Daniel Glen, 15 Nov 2004 */
//https://github.com/afni/afni/blob/b6a9f7a21c1f3231ff09efbd861f8975ad48e525/src/3dDTeig.c
//  Components developed at the US National Institutes of Health (after 15 Jan 2001) are not copyrighted.
void EIG_tsfunc(
                          int npts, float ts[],
                          float * val, int isUpperTriangle )
{
  #define SMALLNUMBER 1E-4
  int i,j;
  // static int nvox , ncall ;
   int maxindex, minindex, midindex;
   float temp, minvalue, maxvalue;
   int sortvector[3];
   double a[9], e[3];
   int astart, vstart;
   double ssq, dsq;
   double dv0, dv1, dv2;

  /* ts is input vector data of 6 floating point numbers.
     For each point in volume brik convert vector data to
     symmetric matrix */
  /* ts should come from data sub-briks in form of Dxx,Dxy,Dxz,Dyy,Dyz,Dzz */
  /* convert to matrix of form 
     [ Dxx Dxy Dxz]
     [ Dxy Dyy Dyz]
     [ Dxz Dyz Dzz]  */

   /** is this a "notification"? **/
   if( val == NULL ){

      if( npts > 0 ){  /* the "start notification" */

         // nvox  = npts ;                       /* keep track of   */
         // ncall = 0 ;                          /* number of calls */

      } else {  /* the "end notification" */

         /* nothing to do here */
      }
      return ;
   }

   /* load the symmetric matrix vector from the "timeseries" subbrik vector values */
   if(isUpperTriangle) {               /* read in as upper diagonal elements */
      a[0]=ts[0]; a[1]=ts[1]; a[2]=ts[2];  
      a[3]=ts[1]; a[4]=ts[3]; a[5]=ts[4];
      a[6]=ts[2]; a[7]=ts[4]; a[8]=ts[5];
   }
   else {         /* read D tensor in as lower diagonal elements - NIFTI standard */ 
      a[0]=ts[0]; a[1]=ts[1]; a[2]=ts[3];
      a[3]=ts[1]; a[4]=ts[2]; a[5]=ts[4];
      a[6]=ts[3]; a[7]=ts[4]; a[8]=ts[5];
   }
  symeig_double(3, a, e);    /* compute eigenvalues in e, eigenvectors in a */
  
  maxindex=2;                      /* find the lowest, middle and highest eigenvalue */
  maxvalue=e[2];
  minindex=0;
  minvalue=e[0];
  midindex = 1;
  for (i=0;i<3;i++) {        
    temp = e[i];
    if(temp>maxvalue) {            /* find the maximum */
      maxindex = i;
      maxvalue = temp;
    }
    if(temp<minvalue) {            /* find the minimum */
      minindex = i;
      minvalue = temp;
    }
  }

  for(i=0;i<3;i++){                /* find the middle */
    if((i!=maxindex) && (i!=minindex))
      {
	midindex = i;
        break;
      }        
  }

  sortvector[0] = maxindex;
  sortvector[1] = midindex;
  sortvector[2] = minindex;

  /* put the eigenvalues at the beginning of the matrix */
  for(i=0;i<3;i++) {
     val[i] = e[sortvector[i]];    /* copy sorted eigenvalues */
  
                                 /* start filling in eigenvector values */
     astart=sortvector[i]*3;    /* start index of double eigenvector */    
     vstart=(i+1)*3;            /* start index of float val vector to 
                                   copy eigenvector */

     for(j=0;j<3;j++){
        val[vstart+j] = a[astart+j];
      }
  }

  for(i=0;i<3;i++) {
    if(fabs(val[i])<SMALLNUMBER)
       val[i] = 0.0;
  }

  /* calculate the Fractional Anisotropy, FA */
  /*   reference, Pierpaoli C, Basser PJ. Microstructural and
       physiological features of tissues elucidated by
       quantitative-diffusion tensor MRI,J Magn Reson B 1996;
       111:209-19 */
  if((val[0]<=0.0)||(val[1]<0.0)||(val[2]<0.0)) {   /* any negative eigenvalues?*/
    val[12]=0.0;                                      /* set FA to 0 */  
    val[13]=0.0;                                      /* set MD to 0 */
    return;//EXRETURN;
  }

  ssq = (val[0]*val[0])+(val[1]*val[1])+(val[2]*val[2]);        /* sum of squares of eigenvalues */
  /* dsq = pow((val[0]-val[1]),2.0) + pow((val[1]-val[2]),2.0) + pow((val[2]-val[0]),2.0);*/ /* sum of differences squared */

  dv0 = val[0]-val[1];
  dv0 *= dv0;
  dv1 = val[1]-val[2];
  dv1 *= dv1;
  dv2 = val[2]-val[0];
  dv2 *= dv2;
  dsq = dv0+dv1+dv2;                 /* sum of differences squared */

  if(ssq!=0.0)
    val[12] = sqrt(dsq/(2.0*ssq));   /* FA calculated here */
  else
    val[12] = 0.0;
  val[13] = (val[0]+val[1]+val[2]) / 3;  /* MD - mean diffusivity=average of eigenvalues */

  return;//EXRETURN;
}
