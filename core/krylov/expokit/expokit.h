#ifndef EXPOKIT_H
#define EXPOKIT_H

#include <core/common.h>

namespace expokit
{

typedef void(*matvecfunc_double)(void *data, double *in, double *out);
typedef void(*matvecfunc_cplx)(void *data, cplx *in, cplx *out);

extern "C"
{
#include "f2c.h"

int dgexpv(int *n, int *m, double *t, 
	double *v, double *w, double *tol, double *anorm, 
	double *wsp, int *lwsp, int *iwsp, int *liwsp, matvecfunc_double 
	matvec, void* matvecdata, int *itrace, int *iflag);

int dsexpv(int *n, int *m, double *t, 
	double *v, double *w, double *tol, double *anorm, 
	double *wsp, int *lwsp, int *iwsp, int *liwsp, matvecfunc_double 
	matvec, void *matvecdata, int *itrace, int *iflag);

int zgexpv(int *n, int *m, double *t, 
	cplx *v, cplx *w, double *tol, double *
	anorm, cplx *wsp, int *lwsp, int *iwsp, int *
	liwsp, matvecfunc_cplx matvec, void *matvecdata, int *itrace, int *iflag);

int zhexpv(int *n, int *m, double *t, 
	cplx *v, cplx *w, double *tol, double *
	anorm, cplx *wsp, int *lwsp, int *iwsp, int *
	liwsp, matvecfunc_cplx matvec, void *matvecdata, int *itrace, int *iflag);

}

};

#endif

