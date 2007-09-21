/* ../dsexpv.f -- translated by f2c (version 20050501).
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

static doublereal c_b4 = 1.;
static integer c__1 = 1;
static doublereal c_b8 = 10.;
static integer c__9 = 9;
static integer c__3 = 3;
static integer c__5 = 5;
static integer c__6 = 6;
static doublereal c_b67 = 0.;

/* ----------------------------------------------------------------------| */
/* Subroutine */ int dsexpv(integer *n, integer *m, doublereal *t, 
	doublereal *v, doublereal *w, doublereal *tol, doublereal *anorm, 
	doublereal *wsp, integer *lwsp, integer *iwsp, integer *liwsp, S_fp 
	matvec, void *matvecdata, integer *itrace, integer *iflag)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);
    double sqrt(doublereal), d_sign(doublereal *, doublereal *), pow_di(
	    doublereal *, integer *), pow_dd(doublereal *, doublereal *), 
	    d_lg10(doublereal *);
    integer i_dnnt(doublereal *);
    double d_int(doublereal *);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle();

    /* Local variables */
    static integer ibrkflag;
    static doublereal step_min__, step_max__;
    static integer i__, j;
    static doublereal break_tol__;
    static integer k1;
    static doublereal p1, p2, p3;
    static integer ih, mh, iv, ns, mx;
    static doublereal xm;
    static integer j1v;
    static doublereal hjj, sgn, eps, hj1j, sqr1, beta;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal hump;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer ifree, lfree;
    static doublereal t_old__;
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static integer iexph;
    static doublereal t_new__;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer nexph;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static doublereal t_now__;
    static integer nstep;
    static doublereal t_out__;
    static integer nmult;
    static doublereal vnorm;
    extern /* Subroutine */ int dgpadm_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *), dnchbv_(integer *, doublereal *,
	     doublereal *, integer *, doublereal *, doublereal *);
    static integer nscale;
    static doublereal rndoff, t_step__, avnorm;
    static integer ireject;
    static doublereal err_loc__;
    static integer nreject, mbrkdwn;
    static doublereal tbrkdwn, s_error__, x_error__;

    /* Fortran I/O blocks */
    static cilist io___40 = { 0, 6, 0, 0, 0 };
    static cilist io___48 = { 0, 6, 0, 0, 0 };
    static cilist io___49 = { 0, 6, 0, 0, 0 };
    static cilist io___50 = { 0, 6, 0, 0, 0 };
    static cilist io___51 = { 0, 6, 0, 0, 0 };
    static cilist io___52 = { 0, 6, 0, 0, 0 };
    static cilist io___53 = { 0, 6, 0, 0, 0 };
    static cilist io___54 = { 0, 6, 0, 0, 0 };
    static cilist io___55 = { 0, 6, 0, 0, 0 };
    static cilist io___56 = { 0, 6, 0, 0, 0 };
    static cilist io___57 = { 0, 6, 0, 0, 0 };
    static cilist io___58 = { 0, 6, 0, 0, 0 };
    static cilist io___59 = { 0, 6, 0, 0, 0 };


/* -----Purpose----------------------------------------------------------| */

/* ---  DSEXPV computes w = exp(t*A)*v - for a Symmetric matrix A. */

/*     It does not compute the matrix exponential in isolation but */
/*     instead, it computes directly the action of the exponential */
/*     operator on the operand vector. This way of doing so allows */
/*     for addressing large sparse problems. */

/*     The method used is based on Krylov subspace projection */
/*     techniques and the matrix under consideration interacts only */
/*     via the external routine `matvec' performing the matrix-vector */
/*     product (matrix-free method). */

/* -----Arguments--------------------------------------------------------| */

/*     n      : (input) order of the principal matrix A. */

/*     m      : (input) maximum size for the Krylov basis. */

/*     t      : (input) time at wich the solution is needed (can be < 0). */

/*     v(n)   : (input) given operand vector. */

/*     w(n)   : (output) computed approximation of exp(t*A)*v. */

/*     tol    : (input/output) the requested accuracy tolerance on w. */
/*              If on input tol=0.0d0 or tol is too small (tol.le.eps) */
/*              the internal value sqrt(eps) is used, and tol is set to */
/*              sqrt(eps) on output (`eps' denotes the machine epsilon). */
/*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol) */

/*     anorm  : (input) an approximation of some norm of A. */

/*   wsp(lwsp): (workspace) lwsp .ge. n*(m+1)+n+(m+2)^2+4*(m+2)^2+ideg+1 */
/*                                   +---------+-------+---------------+ */
/*              (actually, ideg=6)        V        H      wsp for PADE */

/* iwsp(liwsp): (workspace) liwsp .ge. m+2 */

/*     matvec : external subroutine for matrix-vector multiplication. */
/*              synopsis: matvec( x, y ) */
/*                        double precision x(*), y(*) */
/*              computes: y(1:n) <- A*x(1:n) */
/*                        where A is the principal matrix. */

/*     itrace : (input) running mode. 0=silent, 1=print step-by-step info */

/*     iflag  : (output) exit flag. */
/*              <0 - bad input arguments */
/*               0 - no problem */
/*               1 - maximum number of steps reached without convergence */
/*               2 - requested tolerance was too high */

/* -----Accounts on the computation--------------------------------------| */
/*     Upon exit, an interested user may retrieve accounts on the */
/*     computations. They are located in the workspace arrays wsp and */
/*     iwsp as indicated below: */

/*     location  mnemonic                 description */
/*     -----------------------------------------------------------------| */
/*     iwsp(1) = nmult, number of matrix-vector multiplications used */
/*     iwsp(2) = nexph, nbr of Tridiagonal matrix exponential evaluated */
/*     iwsp(3) = nscale, number of repeated squaring involved in Pade */
/*     iwsp(4) = nstep, nbr of integration steps used up to completion */
/*     iwsp(5) = nreject, number of rejected step-sizes */
/*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise */
/*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured */
/*     -----------------------------------------------------------------| */
/*     wsp(1)  = step_min, minimum step-size used during integration */
/*     wsp(2)  = step_max, maximum step-size used during integration */
/*     wsp(3)  = dummy */
/*     wsp(4)  = dummy */
/*     wsp(5)  = x_error, maximum among all local truncation errors */
/*     wsp(6)  = s_error, global sum of local truncation errors */
/*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured */
/*     wsp(8)  = t_now, integration domain successfully covered */
/*     wsp(9)  = hump, i.e., max||exp(sA)||, s in [0,t] (or [t,0] if t<0) */
/*     wsp(10) = ||w||/||v||, scaled norm of the solution w. */
/*     -----------------------------------------------------------------| */
/*     The `hump' is a measure of the conditioning of the problem. The */
/*     matrix exponential is well-conditioned if hump = 1, whereas it is */
/*     poorly-conditioned if hump >> 1. However the solution can still be */
/*     relatively fairly accurate even when the hump is large (the hump */
/*     is an upper bound), especially when the hump and the scaled norm */
/*     of w [this is also computed and returned in wsp(10)] are of the */
/*     same order of magnitude (further details in reference below). */

/* ----------------------------------------------------------------------| */
/* -----The following parameters may also be adjusted herein-------------| */

/*     mxstep  : maximum allowable number of integration steps. */
/*               The value 0 means an infinite number of steps. */

/*     mxreject: maximum allowable number of rejections at each step. */
/*               The value 0 means an infinite number of rejections. */

/*     ideg    : the Pade approximation of type (ideg,ideg) is used as */
/*               an approximation to exp(H). The value 0 switches to the */
/*               uniform rational Chebyshev approximation of type (14,14) */

/*     delta   : local truncation error `safety factor' */

/*     gamma   : stepsize `shrinking factor' */

/* ----------------------------------------------------------------------| */
/*     Roger B. Sidje (rbs@maths.uq.edu.au) */
/*     EXPOKIT: Software Package for Computing Matrix Exponentials. */
/*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998 */
/* ----------------------------------------------------------------------| */

/* ---  check restrictions on input parameters ... */
    /* Parameter adjustments */
    --w;
    --v;
    --wsp;
    --iwsp;

    /* Function Body */
    *iflag = 0;
/* Computing 2nd power */
    i__1 = *m + 2;
    if (*lwsp < *n * (*m + 2) + i__1 * i__1 * 5 + 7) {
	*iflag = -1;
    }
    if (*liwsp < *m + 2) {
	*iflag = -2;
    }
    if (*m >= *n || *m <= 0) {
	*iflag = -3;
    }
    if (*iflag != 0) {
	s_stop("bad sizes (in input of DSEXPV)", (ftnlen)30);
    }

/* ---  initialisations ... */

    k1 = 2;
    mh = *m + 2;
    iv = 1;
    ih = iv + *n * (*m + 1) + *n;
    ifree = ih + mh * mh;
    lfree = *lwsp - ifree + 1;
    ibrkflag = 0;
    mbrkdwn = *m;
    nmult = 0;
    nreject = 0;
    nexph = 0;
    nscale = 0;
    t_out__ = abs(*t);
    tbrkdwn = 0.;
    step_min__ = t_out__;
    step_max__ = 0.;
    nstep = 0;
    s_error__ = 0.;
    x_error__ = 0.;
    t_now__ = 0.;
    t_new__ = 0.;
    p1 = 1.3333333333333333;
L1:
    p2 = p1 - 1.;
    p3 = p2 + p2 + p2;
    eps = (d__1 = p3 - 1., abs(d__1));
    if (eps == 0.) {
	goto L1;
    }
    if (*tol <= eps) {
	*tol = sqrt(eps);
    }
    rndoff = eps * *anorm;
    break_tol__ = 1e-7;
/* >>>  break_tol = tol */
/* >>>  break_tol = anorm*tol */
    sgn = d_sign(&c_b4, t);
    dcopy_(n, &v[1], &c__1, &w[1], &c__1);
    beta = dnrm2_(n, &w[1], &c__1);
    vnorm = beta;
    hump = beta;

/* ---  obtain the very first stepsize ... */

    sqr1 = sqrt(.1);
    xm = 1. / (doublereal) (*m);
    d__1 = (*m + 1) / 2.72;
    i__1 = *m + 1;
    p1 = *tol * pow_di(&d__1, &i__1) * sqrt((*m + 1) * 6.2800000000000002);
    d__1 = p1 / (beta * 4. * *anorm);
    t_new__ = 1. / *anorm * pow_dd(&d__1, &xm);
    d__1 = d_lg10(&t_new__) - sqr1;
    i__1 = i_dnnt(&d__1) - 1;
    p1 = pow_di(&c_b8, &i__1);
    d__1 = t_new__ / p1 + .55;
    t_new__ = d_int(&d__1) * p1;

/* ---  step-by-step integration ... */

L100:
    if (t_now__ >= t_out__) {
	goto L500;
    }
    ++nstep;
/* Computing MIN */
    d__1 = t_out__ - t_now__;
    t_step__ = min(d__1,t_new__);
    p1 = 1. / beta;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wsp[iv + i__ - 1] = p1 * w[i__];
    }
    i__1 = mh * mh;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wsp[ih + i__ - 1] = 0.;
    }

/* ---  Lanczos loop ... */

    j1v = iv + *n;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	++nmult;
	(*matvec)(matvecdata, &wsp[j1v - *n], &wsp[j1v]);
	if (j > 1) {
	    d__1 = -wsp[ih + (j - 1) * mh + j - 2];
	    daxpy_(n, &d__1, &wsp[j1v - (*n << 1)], &c__1, &wsp[j1v], &c__1);
	}
	hjj = ddot_(n, &wsp[j1v - *n], &c__1, &wsp[j1v], &c__1);
	d__1 = -hjj;
	daxpy_(n, &d__1, &wsp[j1v - *n], &c__1, &wsp[j1v], &c__1);
	hj1j = dnrm2_(n, &wsp[j1v], &c__1);
	wsp[ih + (j - 1) * (mh + 1)] = hjj;
/* ---     if `happy breakdown' go straightforward at the end ... */
	if (hj1j <= break_tol__) {
	    s_wsle(&io___40);
	    do_lio(&c__9, &c__1, "happy breakdown: mbrkdwn =", (ftnlen)26);
	    do_lio(&c__3, &c__1, (char *)&j, (ftnlen)sizeof(integer));
	    do_lio(&c__9, &c__1, " h =", (ftnlen)4);
	    do_lio(&c__5, &c__1, (char *)&hj1j, (ftnlen)sizeof(doublereal));
	    e_wsle();
	    k1 = 0;
	    ibrkflag = 1;
	    mbrkdwn = j;
	    tbrkdwn = t_now__;
	    t_step__ = t_out__ - t_now__;
	    goto L300;
	}
	wsp[ih + (j - 1) * mh + j] = hj1j;
	wsp[ih + j * mh + j - 1] = hj1j;
	d__1 = 1. / hj1j;
	dscal_(n, &d__1, &wsp[j1v], &c__1);
	j1v += *n;
/* L200: */
    }
    ++nmult;
    (*matvec)(matvecdata, &wsp[j1v - *n], &wsp[j1v]);
    avnorm = dnrm2_(n, &wsp[j1v], &c__1);

/* ---  set 1 for the 2-corrected scheme ... */

L300:
    wsp[ih + *m * mh + *m - 1] = 0.;
    wsp[ih + *m * mh + *m + 1] = 1.;

/* ---  loop while ireject<mxreject until the tolerance is reached ... */

    ireject = 0;
L401:

/* ---  compute w = beta*V*exp(t_step*H)*e1 ... */

    ++nexph;
    mx = mbrkdwn + k1;
    if (TRUE_) {
/* ---     irreducible rational Pade approximation ... */
	d__1 = sgn * t_step__;
	dgpadm_(&c__6, &mx, &d__1, &wsp[ih], &mh, &wsp[ifree], &lfree, &iwsp[
		1], &iexph, &ns, iflag);
	iexph = ifree + iexph - 1;
	nscale += ns;
    } else {
/* ---     uniform rational Chebyshev approximation ... */
	iexph = ifree;
	i__1 = mx;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    wsp[iexph + i__ - 1] = 0.;
	}
	wsp[iexph] = 1.;
	d__1 = sgn * t_step__;
	dnchbv_(&mx, &d__1, &wsp[ih], &mh, &wsp[iexph], &wsp[ifree + mx]);
    }
/* L402: */

/* ---  error estimate ... */

    if (k1 == 0) {
	err_loc__ = *tol;
    } else {
	p1 = (d__1 = wsp[iexph + *m], abs(d__1)) * beta;
	p2 = (d__1 = wsp[iexph + *m + 1], abs(d__1)) * beta * avnorm;
	if (p1 > p2 * 10.) {
	    err_loc__ = p2;
	    xm = 1. / (doublereal) (*m);
	} else if (p1 > p2) {
	    err_loc__ = p1 * p2 / (p1 - p2);
	    xm = 1. / (doublereal) (*m);
	} else {
	    err_loc__ = p1;
	    xm = 1. / (doublereal) (*m - 1);
	}
    }

/* ---  reject the step-size if the error is not acceptable ... */

    if (k1 != 0 && err_loc__ > t_step__ * 1.2 * *tol) {
	t_old__ = t_step__;
	d__1 = t_step__ * *tol / err_loc__;
	t_step__ = t_step__ * .9 * pow_dd(&d__1, &xm);
	d__1 = d_lg10(&t_step__) - sqr1;
	i__1 = i_dnnt(&d__1) - 1;
	p1 = pow_di(&c_b8, &i__1);
	d__1 = t_step__ / p1 + .55;
	t_step__ = d_int(&d__1) * p1;
	if (*itrace != 0) {
	    s_wsle(&io___48);
	    do_lio(&c__9, &c__1, "t_step =", (ftnlen)8);
	    do_lio(&c__5, &c__1, (char *)&t_old__, (ftnlen)sizeof(doublereal))
		    ;
	    e_wsle();
	    s_wsle(&io___49);
	    do_lio(&c__9, &c__1, "err_loc =", (ftnlen)9);
	    do_lio(&c__5, &c__1, (char *)&err_loc__, (ftnlen)sizeof(
		    doublereal));
	    e_wsle();
	    s_wsle(&io___50);
	    do_lio(&c__9, &c__1, "err_required =", (ftnlen)14);
	    d__1 = t_old__ * 1.2 * *tol;
	    do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	    e_wsle();
	    s_wsle(&io___51);
	    do_lio(&c__9, &c__1, "stepsize rejected, stepping down to:", (
		    ftnlen)36);
	    do_lio(&c__5, &c__1, (char *)&t_step__, (ftnlen)sizeof(doublereal)
		    );
	    e_wsle();
	}
	++ireject;
	++nreject;
	if (FALSE_) {
	    s_wsle(&io___52);
	    do_lio(&c__9, &c__1, "Failure in DSEXPV: ---", (ftnlen)22);
	    e_wsle();
	    s_wsle(&io___53);
	    do_lio(&c__9, &c__1, "The requested tolerance is too high.", (
		    ftnlen)36);
	    e_wsle();
	    s_wsle(&io___54);
	    do_lio(&c__9, &c__1, "Rerun with a smaller value.", (ftnlen)27);
	    e_wsle();
	    *iflag = 2;
	    return 0;
	}
	goto L401;
    }

/* ---  now update w = beta*V*exp(t_step*H)*e1 and the hump ... */

/* Computing MAX */
    i__1 = 0, i__2 = k1 - 1;
    mx = mbrkdwn + max(i__1,i__2);
    dgemv_("n", n, &mx, &beta, &wsp[iv], n, &wsp[iexph], &c__1, &c_b67, &w[1],
	     &c__1, (ftnlen)1);
    beta = dnrm2_(n, &w[1], &c__1);
    hump = max(hump,beta);

/* ---  suggested value for the next stepsize ... */

    d__1 = t_step__ * *tol / err_loc__;
    t_new__ = t_step__ * .9 * pow_dd(&d__1, &xm);
    d__1 = d_lg10(&t_new__) - sqr1;
    i__1 = i_dnnt(&d__1) - 1;
    p1 = pow_di(&c_b8, &i__1);
    d__1 = t_new__ / p1 + .55;
    t_new__ = d_int(&d__1) * p1;
    err_loc__ = max(err_loc__,rndoff);

/* ---  update the time covered ... */

    t_now__ += t_step__;

/* ---  display and keep some information ... */

    if (*itrace != 0) {
	s_wsle(&io___55);
	do_lio(&c__9, &c__1, "integration", (ftnlen)11);
	do_lio(&c__3, &c__1, (char *)&nstep, (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, "---------------------------------", (ftnlen)33);
	e_wsle();
	s_wsle(&io___56);
	do_lio(&c__9, &c__1, "scale-square =", (ftnlen)14);
	do_lio(&c__3, &c__1, (char *)&ns, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___57);
	do_lio(&c__9, &c__1, "step_size =", (ftnlen)11);
	do_lio(&c__5, &c__1, (char *)&t_step__, (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___58);
	do_lio(&c__9, &c__1, "err_loc   =", (ftnlen)11);
	do_lio(&c__5, &c__1, (char *)&err_loc__, (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___59);
	do_lio(&c__9, &c__1, "next_step =", (ftnlen)11);
	do_lio(&c__5, &c__1, (char *)&t_new__, (ftnlen)sizeof(doublereal));
	e_wsle();
    }
    step_min__ = min(step_min__,t_step__);
    step_max__ = max(step_max__,t_step__);
    s_error__ += err_loc__;
    x_error__ = max(x_error__,err_loc__);
    if (nstep < 500) {
	goto L100;
    }
    *iflag = 1;
L500:
    iwsp[1] = nmult;
    iwsp[2] = nexph;
    iwsp[3] = nscale;
    iwsp[4] = nstep;
    iwsp[5] = nreject;
    iwsp[6] = ibrkflag;
    iwsp[7] = mbrkdwn;
    wsp[1] = step_min__;
    wsp[2] = step_max__;
    wsp[3] = 0.;
    wsp[4] = 0.;
    wsp[5] = x_error__;
    wsp[6] = s_error__;
    wsp[7] = tbrkdwn;
    wsp[8] = sgn * t_now__;
    wsp[9] = hump / vnorm;
    wsp[10] = beta / vnorm;
    return 0;
} /* dsexpv_ */

#ifdef __cplusplus
	}
#endif
