#ifndef PCAUPP_H
#define PCAUPP_H

#include <mpi.h>
#include "arpackwrapper.h"
#include "parpackf.h"

inline void pcaupp(MPI_Comm comm, int& ido, char bmat, int n, char* which, int nev,
                  double& tol, arcomplex<double> resid[], int ncv,
                  arcomplex<double> V[], int ldv, int iparam[], int ipntr[],
                  arcomplex<double> workd[], arcomplex<double> workl[],
                  int lworkl, double rwork[], int& info)
{
  MPI_Fint commHandle = MPI_Comm_c2f(comm);
  F77NAME(pznaupd)(&commHandle, &ido, &bmat, &n, which, &nev, &tol, resid, &ncv,
                  &V[0], &ldv, &iparam[0], &ipntr[0], &workd[0],
                  &workl[0], &lworkl, &rwork[0], &info);

}

inline void pcaupp(MPI_Comm comm, int& ido, char bmat, int n, char* which, int nev,
                  float& tol, arcomplex<float> resid[], int ncv,
                  arcomplex<float> V[], int ldv, int iparam[], int ipntr[],
                  arcomplex<float> workd[], arcomplex<float> workl[],
                  int lworkl, float rwork[], int& info)
{
  MPI_Fint commHandle = MPI_Comm_c2f(comm);
  F77NAME(pcnaupd)(&commHandle, &ido, &bmat, &n, which, &nev, &tol, resid, &ncv,
                  &V[1], &ldv, &iparam[0], &ipntr[0], &workd[0],
                  &workl[0], &lworkl, &rwork[0], &info);

} // caupp (arcomplex<float>).

#endif // PCAUPP_H



