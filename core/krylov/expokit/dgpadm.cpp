/* ../dgpadm.f -- translated by f2c (version 20050501).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

/* Table of constant values */

static integer c__2 = 2;
static doublereal c_b7 = 0.;
static doublereal c_b11 = 1.;
static doublereal c_b19 = -1.;
static integer c__1 = 1;
static doublereal c_b23 = 2.;

/* ----------------------------------------------------------------------| */
/* Subroutine */ int dgpadm_(integer *ideg, integer *m, doublereal *t, 
	doublereal *h__, integer *ldh, doublereal *wsp, integer *lwsp, 
	integer *ipiv, integer *iexph, integer *ns, integer *iflag)
{
    /* System generated locals */
    integer h_dim1, h_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);
    double log(doublereal);
    integer pow_ii(integer *, integer *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal cp, cq;
    static integer ip, mm, iq, ih2, iodd, iget, iput, icoef;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal scale;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer ifree;
    extern /* Subroutine */ int dgesv_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *);
    static integer iused;
    static doublereal hnorm;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static doublereal scale2;

/* -----Purpose----------------------------------------------------------| */

/*     Computes exp(t*H), the matrix exponential of a general matrix in */
/*     full, using the irreducible rational Pade approximation to the */
/*     exponential function exp(x) = r(x) = (+/-)( I + 2*(q(x)/p(x)) ), */
/*     combined with scaling-and-squaring. */

/* -----Arguments--------------------------------------------------------| */

/*     ideg      : (input) the degre of the diagonal Pade to be used. */
/*                 a value of 6 is generally satisfactory. */

/*     m         : (input) order of H. */

/*     H(ldh,m)  : (input) argument matrix. */

/*     t         : (input) time-scale (can be < 0). */

/*     wsp(lwsp) : (workspace/output) lwsp .ge. 4*m*m+ideg+1. */

/*     ipiv(m)   : (workspace) */

/* >>>> iexph     : (output) number such that wsp(iexph) points to exp(tH) */
/*                 i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1) */
/*                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */
/*                 NOTE: if the routine was called with wsp(iptr), */
/*                       then exp(tH) will start at wsp(iptr+iexph-1). */

/*     ns        : (output) number of scaling-squaring used. */

/*     iflag     : (output) exit flag. */
/*                      0 - no problem */
/*                     <0 - problem */

/* ----------------------------------------------------------------------| */
/*     Roger B. Sidje (rbs@maths.uq.edu.au) */
/*     EXPOKIT: Software Package for Computing Matrix Exponentials. */
/*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998 */
/* ----------------------------------------------------------------------| */

/* ---  check restrictions on input parameters ... */
    /* Parameter adjustments */
    --ipiv;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --wsp;

    /* Function Body */
    mm = *m * *m;
    *iflag = 0;
    if (*ldh < *m) {
	*iflag = -1;
    }
    if (*lwsp < (mm << 2) + *ideg + 1) {
	*iflag = -2;
    }
    if (*iflag != 0) {
	s_stop("bad sizes (in input of DGPADM)", (ftnlen)30);
    }

/* ---  initialise pointers ... */

    icoef = 1;
    ih2 = icoef + (*ideg + 1);
    ip = ih2 + mm;
    iq = ip + mm;
    ifree = iq + mm;

/* ---  scaling: seek ns such that ||t*H/2^ns|| < 1/2; */
/*     and set scale = t/2^ns ... */

    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wsp[i__] = 0.;
    }
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    wsp[i__] += (d__1 = h__[i__ + j * h_dim1], abs(d__1));
	}
    }
    hnorm = 0.;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	d__1 = hnorm, d__2 = wsp[i__];
	hnorm = max(d__1,d__2);
    }
    hnorm = (d__1 = *t * hnorm, abs(d__1));
    if (hnorm == 0.) {
	s_stop("Error - null H in input of DGPADM.", (ftnlen)34);
    }
/* Computing MAX */
    i__1 = 0, i__2 = (integer) (log(hnorm) / log(2.)) + 2;
    *ns = max(i__1,i__2);
    scale = *t / (doublereal) pow_ii(&c__2, ns);
    scale2 = scale * scale;

/* ---  compute Pade coefficients ... */

    i__ = *ideg + 1;
    j = (*ideg << 1) + 1;
    wsp[icoef] = 1.;
    i__1 = *ideg;
    for (k = 1; k <= i__1; ++k) {
	wsp[icoef + k] = wsp[icoef + k - 1] * (doublereal) (i__ - k) / (
		doublereal) (k * (j - k));
    }

/* ---  H2 = scale2*H*H ... */

    dgemm_("n", "n", m, m, m, &scale2, &h__[h_offset], ldh, &h__[h_offset], 
	    ldh, &c_b7, &wsp[ih2], m, (ftnlen)1, (ftnlen)1);

/* ---  initialize p (numerator) and q (denominator) ... */

    cp = wsp[icoef + *ideg - 1];
    cq = wsp[icoef + *ideg];
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    wsp[ip + (j - 1) * *m + i__ - 1] = 0.;
	    wsp[iq + (j - 1) * *m + i__ - 1] = 0.;
	}
	wsp[ip + (j - 1) * (*m + 1)] = cp;
	wsp[iq + (j - 1) * (*m + 1)] = cq;
    }

/* ---  Apply Horner rule ... */

    iodd = 1;
    k = *ideg - 1;
L100:
    iused = iodd * iq + (1 - iodd) * ip;
    dgemm_("n", "n", m, m, m, &c_b11, &wsp[iused], m, &wsp[ih2], m, &c_b7, &
	    wsp[ifree], m, (ftnlen)1, (ftnlen)1);
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	wsp[ifree + (j - 1) * (*m + 1)] += wsp[icoef + k - 1];
    }
    ip = (1 - iodd) * ifree + iodd * ip;
    iq = iodd * ifree + (1 - iodd) * iq;
    ifree = iused;
    iodd = 1 - iodd;
    --k;
    if (k > 0) {
	goto L100;
    }

/* ---  Obtain (+/-)(I + 2*(p\q)) ... */

    if (iodd == 1) {
	dgemm_("n", "n", m, m, m, &scale, &wsp[iq], m, &h__[h_offset], ldh, &
		c_b7, &wsp[ifree], m, (ftnlen)1, (ftnlen)1);
	iq = ifree;
    } else {
	dgemm_("n", "n", m, m, m, &scale, &wsp[ip], m, &h__[h_offset], ldh, &
		c_b7, &wsp[ifree], m, (ftnlen)1, (ftnlen)1);
	ip = ifree;
    }
    daxpy_(&mm, &c_b19, &wsp[ip], &c__1, &wsp[iq], &c__1);
    dgesv_(m, m, &wsp[iq], m, &ipiv[1], &wsp[ip], m, iflag);
    if (*iflag != 0) {
	s_stop("Problem in DGESV (within DGPADM)", (ftnlen)32);
    }
    dscal_(&mm, &c_b23, &wsp[ip], &c__1);
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	wsp[ip + (j - 1) * (*m + 1)] += 1.;
    }
    iput = ip;
    if (*ns == 0 && iodd == 1) {
	dscal_(&mm, &c_b19, &wsp[ip], &c__1);
	goto L200;
    }

/* --   squaring : exp(t*H) = (exp(t*H))^(2^ns) ... */

    iodd = 1;
    i__1 = *ns;
    for (k = 1; k <= i__1; ++k) {
	iget = iodd * ip + (1 - iodd) * iq;
	iput = (1 - iodd) * ip + iodd * iq;
	dgemm_("n", "n", m, m, m, &c_b11, &wsp[iget], m, &wsp[iget], m, &c_b7,
		 &wsp[iput], m, (ftnlen)1, (ftnlen)1);
	iodd = 1 - iodd;
    }
L200:
    *iexph = iput;
    return 0;
} /* dgpadm_ */

#ifdef __cplusplus
	}
#endif
