/* ../zgpadm.f -- translated by f2c (version 20050501).
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

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__2 = 2;
static integer c__1 = 1;
static doublereal c_b19 = 2.;
static doublereal c_b21 = -1.;

/* ----------------------------------------------------------------------| */
/* Subroutine */ int zgpadm_(integer *ideg, integer *m, doublereal *t, 
	doublecomplex *h__, integer *ldh, doublecomplex *wsp, integer *lwsp, 
	integer *ipiv, integer *iexph, integer *ns, integer *iflag)
{
    /* System generated locals */
    integer h_dim1, h_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);
    double z_abs(doublecomplex *), log(doublereal);
    integer pow_ii(integer *, integer *);

    /* Local variables */
    static integer i__, j, k;
    static doublecomplex cp, cq;
    static integer ip, mm, iq, ih2, iodd, iget, iput, icoef;
    static doublecomplex scale;
    static integer ifree, iused;
    extern /* Subroutine */ int zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static doublereal hnorm;
    extern /* Subroutine */ int zgesv_(integer *, integer *, doublecomplex *, 
	    integer *, integer *, doublecomplex *, integer *, integer *);
    static doublecomplex scale2;
    extern /* Subroutine */ int zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), zdscal_(
	    integer *, doublereal *, doublecomplex *, integer *);

/* -----Purpose----------------------------------------------------------| */

/*     Computes exp(t*H), the matrix exponential of a general complex */
/*     matrix in full, using the irreducible rational Pade approximation */
/*     to the exponential exp(z) = r(z) = (+/-)( I + 2*(q(z)/p(z)) ), */
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
/*                       0 - no problem */
/*                      <0 - problem */

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
	s_stop("bad sizes (in input of ZGPADM)", (ftnlen)30);
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
	i__2 = i__;
	wsp[i__2].r = 0., wsp[i__2].i = 0.;
    }
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    i__4 = i__;
	    d__1 = z_abs(&h__[i__ + j * h_dim1]);
	    z__1.r = wsp[i__4].r + d__1, z__1.i = wsp[i__4].i;
	    wsp[i__3].r = z__1.r, wsp[i__3].i = z__1.i;
	}
    }
    hnorm = 0.;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	i__2 = i__;
	d__1 = hnorm, d__2 = wsp[i__2].r;
	hnorm = max(d__1,d__2);
    }
    hnorm = (d__1 = *t * hnorm, abs(d__1));
    if (hnorm == 0.) {
	s_stop("Error - null H in input of ZGPADM.", (ftnlen)34);
    }
/* Computing MAX */
    i__1 = 0, i__2 = (integer) (log(hnorm) / log(2.)) + 2;
    *ns = max(i__1,i__2);
    d__1 = *t / (doublereal) pow_ii(&c__2, ns);
    z__1.r = d__1, z__1.i = 0.;
    scale.r = z__1.r, scale.i = z__1.i;
    z__1.r = scale.r * scale.r - scale.i * scale.i, z__1.i = scale.r * 
	    scale.i + scale.i * scale.r;
    scale2.r = z__1.r, scale2.i = z__1.i;

/* ---  compute Pade coefficients ... */

    i__ = *ideg + 1;
    j = (*ideg << 1) + 1;
    i__1 = icoef;
    wsp[i__1].r = 1., wsp[i__1].i = 0.;
    i__1 = *ideg;
    for (k = 1; k <= i__1; ++k) {
	i__2 = icoef + k;
	i__3 = icoef + k - 1;
	d__1 = (doublereal) (i__ - k);
	z__2.r = d__1 * wsp[i__3].r, z__2.i = d__1 * wsp[i__3].i;
	d__2 = (doublereal) (k * (j - k));
	z__1.r = z__2.r / d__2, z__1.i = z__2.i / d__2;
	wsp[i__2].r = z__1.r, wsp[i__2].i = z__1.i;
    }

/* ---  H2 = scale2*H*H ... */

    zgemm_("n", "n", m, m, m, &scale2, &h__[h_offset], ldh, &h__[h_offset], 
	    ldh, &c_b1, &wsp[ih2], m, (ftnlen)1, (ftnlen)1);

/* ---  initialise p (numerator) and q (denominator) ... */

    i__1 = icoef + *ideg - 1;
    cp.r = wsp[i__1].r, cp.i = wsp[i__1].i;
    i__1 = icoef + *ideg;
    cq.r = wsp[i__1].r, cq.i = wsp[i__1].i;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = ip + (j - 1) * *m + i__ - 1;
	    wsp[i__3].r = 0., wsp[i__3].i = 0.;
	    i__3 = iq + (j - 1) * *m + i__ - 1;
	    wsp[i__3].r = 0., wsp[i__3].i = 0.;
	}
	i__2 = ip + (j - 1) * (*m + 1);
	wsp[i__2].r = cp.r, wsp[i__2].i = cp.i;
	i__2 = iq + (j - 1) * (*m + 1);
	wsp[i__2].r = cq.r, wsp[i__2].i = cq.i;
    }

/* ---  Apply Horner rule ... */

    iodd = 1;
    k = *ideg - 1;
L100:
    iused = iodd * iq + (1 - iodd) * ip;
    zgemm_("n", "n", m, m, m, &c_b2, &wsp[iused], m, &wsp[ih2], m, &c_b1, &
	    wsp[ifree], m, (ftnlen)1, (ftnlen)1);
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = ifree + (j - 1) * (*m + 1);
	i__3 = ifree + (j - 1) * (*m + 1);
	i__4 = icoef + k - 1;
	z__1.r = wsp[i__3].r + wsp[i__4].r, z__1.i = wsp[i__3].i + wsp[i__4]
		.i;
	wsp[i__2].r = z__1.r, wsp[i__2].i = z__1.i;
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

    if (iodd != 0) {
	zgemm_("n", "n", m, m, m, &scale, &wsp[iq], m, &h__[h_offset], ldh, &
		c_b1, &wsp[ifree], m, (ftnlen)1, (ftnlen)1);
	iq = ifree;
    } else {
	zgemm_("n", "n", m, m, m, &scale, &wsp[ip], m, &h__[h_offset], ldh, &
		c_b1, &wsp[ifree], m, (ftnlen)1, (ftnlen)1);
	ip = ifree;
    }
    z__1.r = -1., z__1.i = -0.;
    zaxpy_(&mm, &z__1, &wsp[ip], &c__1, &wsp[iq], &c__1);
    zgesv_(m, m, &wsp[iq], m, &ipiv[1], &wsp[ip], m, iflag);
    if (*iflag != 0) {
	s_stop("Problem in ZGESV (within ZGPADM)", (ftnlen)32);
    }
    zdscal_(&mm, &c_b19, &wsp[ip], &c__1);
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = ip + (j - 1) * (*m + 1);
	i__3 = ip + (j - 1) * (*m + 1);
	z__1.r = wsp[i__3].r + 1., z__1.i = wsp[i__3].i + 0.;
	wsp[i__2].r = z__1.r, wsp[i__2].i = z__1.i;
    }
    iput = ip;
    if (*ns == 0 && iodd != 0) {
	zdscal_(&mm, &c_b21, &wsp[ip], &c__1);
	goto L200;
    }

/* --   squaring : exp(t*H) = (exp(t*H))^(2^ns) ... */

    iodd = 1;
    i__1 = *ns;
    for (k = 1; k <= i__1; ++k) {
	iget = iodd * ip + (1 - iodd) * iq;
	iput = (1 - iodd) * ip + iodd * iq;
	zgemm_("n", "n", m, m, m, &c_b2, &wsp[iget], m, &wsp[iget], m, &c_b1, 
		&wsp[iput], m, (ftnlen)1, (ftnlen)1);
	iodd = 1 - iodd;
    }
L200:
    *iexph = iput;
    return 0;
} /* zgpadm_ */

#ifdef __cplusplus
	}
#endif
