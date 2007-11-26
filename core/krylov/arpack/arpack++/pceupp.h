#ifndef CEUPP_H
#define CEUPP_H

#include <mpi.h>
#include <stddef.h>
#include "arpackwrapper.h"
#include "parpackf.h"

inline void pceupp(MPI_Comm comm, bool rvec, char HowMny, arcomplex<double> d[],
                  arcomplex<double> Z[], int ldz, arcomplex<double> sigma,
                  arcomplex<double> workev[], char bmat, int n, char* which,
                  int nev, double tol, arcomplex<double> resid[], int ncv,
                  arcomplex<double> V[], int ldv, int iparam[], int ipntr[],
                  arcomplex<double> workd[], arcomplex<double> workl[],
                  int lworkl, double rwork[], int& info)
{

  int                irvec;
  logical*           iselect;
  arcomplex<double>* iZ;

  irvec   = (int) rvec;
  iselect = new logical[ncv];
  iZ = (Z == NULL) ? &V[0] : Z;

  MPI_Fint commHandle = MPI_Comm_c2f(comm);

  F77NAME(pzneupd)(&commHandle, &irvec, &HowMny, iselect, d, iZ, &ldz, &sigma,
                  &workev[0], &bmat, &n, which, &nev, &tol, resid,
                  &ncv, &V[0], &ldv, &iparam[0], &ipntr[0],
                  &workd[0], &workl[0], &lworkl, &rwork[0], &info);

  delete[] iselect;

} // ceupp (arcomplex<double>).

inline void pceupp(MPI_Comm comm, bool rvec, char HowMny, arcomplex<float> d[],
                  arcomplex<float> Z[], int ldz, arcomplex<float> sigma,
                  arcomplex<float> workev[], char bmat, int n, char* which,
                  int nev, float tol, arcomplex<float> resid[], int ncv,
                  arcomplex<float> V[], int ldv, int iparam[], int ipntr[],
                  arcomplex<float> workd[], arcomplex<float> workl[],
                  int lworkl, float rwork[], int& info)

/*
  c++ version of ARPACK routine cneupd. The only difference between
  cneupd and zneupd is that in the former function all vectors have
  single precision elements and in the latter all vectors have double
  precision elements.
*/

{

  int               irvec;
  logical*          iselect;
  arcomplex<float>* iZ;

  irvec   = (int) rvec;
  iselect = new logical[ncv];
  iZ = (Z == NULL) ? &V[0] : Z;

  MPI_Fint commHandle = MPI_Comm_c2f(comm);

  F77NAME(pcneupd)(&commHandle, &irvec, &HowMny, iselect, d, iZ, &ldz, &sigma,
                  &workev[0], &bmat, &n, which, &nev, &tol, resid,
                  &ncv, &V[0], &ldv, &iparam[0], &ipntr[0],
                  &workd[0], &workl[0], &lworkl, &rwork[0], &info);

  delete[] iselect;

} // ceupp (arcomplex<float>).

#endif // CEUPP_H
