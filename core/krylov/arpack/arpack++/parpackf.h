#ifndef PARPACKF_H
#define PARPACKF_H

#include <mpi.h>
#include "arpackwrapper.h"

extern "C"
{

  void F77NAME(pdsaupd)(MPI_Fint *comm, integer *ido, char *bmat, integer *n, char *which,
                       integer *nev, double *tol, double *resid,
                       integer *ncv, double *V, integer *ldv,
                       integer *iparam, integer *ipntr, double *workd,
                       double *workl, integer *lworkl, integer *info);

  void F77NAME(pdseupd)(MPI_Fint *comm, logical *rvec, char *HowMny, logical *select,
                       double *d, double *Z, integer *ldz,
                       double *sigma, char *bmat, integer *n,
                       char *which, integer *nev, double *tol,
                       double *resid, integer *ncv, double *V,
                       integer *ldv, integer *iparam, integer *ipntr,
                       double *workd, double *workl,
                       integer *lworkl, integer *info);

// double precision nonsymmetric routines.

  void F77NAME(pdnaupd)(MPI_Fint *comm, integer *ido, char *bmat, integer *n, char *which,
                       integer *nev, double *tol, double *resid,
                       integer *ncv, double *V, integer *ldv,
                       integer *iparam, integer *ipntr, double *workd,
                       double *workl, integer *lworkl, integer *info);

  void F77NAME(pdneupd)(MPI_Fint *comm, logical *rvec, char *HowMny, logical *select,
                       double *dr, double *di, double *Z,
                       integer *ldz, double *sigmar,
                       double *sigmai, double *workev,
                       char *bmat, integer *n, char *which,
                       integer *nev, double *tol, double *resid,
                       integer *ncv, double *V, integer *ldv,
                       integer *iparam, integer *ipntr,
                       double *workd, double *workl,
                       integer *lworkl, integer *info);

// single precision symmetric routines.

  void F77NAME(pssaupd)(MPI_Fint *comm, integer *ido, char *bmat, integer *n, char *which,
                       integer *nev, float *tol, float *resid,
                       integer *ncv, float *V, integer *ldv,
                       integer *iparam, integer *ipntr, float *workd,
                       float *workl, integer *lworkl, integer *info);

  void F77NAME(psseupd)(MPI_Fint *comm, logical *rvec, char *HowMny, logical *select,
                       float *d, float *Z, integer *ldz,
                       float *sigma, char *bmat, integer *n,
                       char *which, integer *nev, float *tol,
                       float *resid, integer *ncv, float *V,
                       integer *ldv, integer *iparam, integer *ipntr,
                       float *workd, float *workl,
                       integer *lworkl, integer *info);

// single precision nonsymmetric routines.

  void F77NAME(psnaupd)(MPI_Fint *comm, integer *ido, char *bmat, integer *n, char *which,
                       integer *nev, float *tol, float *resid,
                       integer *ncv, float *V, integer *ldv,
                       integer *iparam, integer *ipntr, float *workd,
                       float *workl, integer *lworkl, integer *info);

  void F77NAME(psneupd)(MPI_Fint *comm, logical *rvec, char *HowMny, logical *select,
                       float *dr, float *di, float *Z,
                       integer *ldz, float *sigmar,
                       float *sigmai, float *workev, char *bmat,
                       integer *n, char *which, integer *nev,
                       float *tol, float *resid, integer *ncv,
                       float *V, integer *ldv, integer *iparam,
                       integer *ipntr, float *workd, float *workl,
                       integer *lworkl, integer *info);

#ifdef ARCOMP_H

// single precision complex routines.

  void F77NAME(pcnaupd)(MPI_Fint *comm, integer *ido, char *bmat, integer *n, char *which,
                       integer *nev, float *tol, arcomplex<float> *resid,
                       integer *ncv, arcomplex<float> *V, integer *ldv,
                       integer *iparam, integer *ipntr, arcomplex<float> *workd,
                       arcomplex<float> *workl, integer *lworkl,
                       float *rwork, integer *info);

  void F77NAME(pcneupd)(MPI_Fint *comm, logical *rvec, char *HowMny, logical *select,
                       arcomplex<float> *d, arcomplex<float> *Z, integer *ldz,
                       arcomplex<float> *sigma, arcomplex<float> *workev,
                       char *bmat, integer *n, char *which, integer *nev,
                       float *tol, arcomplex<float> *resid, integer *ncv,
                       arcomplex<float> *V, integer *ldv, integer *iparam,
                       integer *ipntr, arcomplex<float> *workd,
                       arcomplex<float> *workl, integer *lworkl,
                       float *rwork, integer *info);

// double precision complex routines.

  void F77NAME(pznaupd)(MPI_Fint *comm, integer *ido, char *bmat, integer *n, char *which,
                       integer *nev, double *tol, arcomplex<double> *resid,
                       integer *ncv, arcomplex<double> *V, integer *ldv,
                       integer *iparam, integer *ipntr, arcomplex<double> *workd,
                       arcomplex<double> *workl, integer *lworkl,
                       double *rwork, integer *info);

  void F77NAME(pzneupd)(MPI_Fint *comm, logical *rvec, char *HowMny, logical *select,
                       arcomplex<double> *d, arcomplex<double> *Z, integer *ldz,
                       arcomplex<double> *sigma, arcomplex<double> *workev,
                       char *bmat, integer *n, char *which, integer *nev,
                       double *tol, arcomplex<double> *resid, integer *ncv,
                       arcomplex<double> *V, integer *ldv, integer *iparam,
                       integer *ipntr, arcomplex<double> *workd,
                       arcomplex<double> *workl, integer *lworkl,
                       double *rwork, integer *info);

}

#endif // ARCOMP_H

#endif // PARPACKF_H
