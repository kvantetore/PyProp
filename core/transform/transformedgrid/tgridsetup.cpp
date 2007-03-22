#include <stdio.h>
#include <iostream>

#include "tools.h"
#include "../../utility/blitzblas.h"

#define FORTRAN_NAME(x) x

namespace TransformedGrid
{

extern "C"
{
   void FORTRAN_NAME(dgeevx)(char *BALANC,char *JOBVL,char *JOBVR,char *SENSE,int *N,
                 double *A,int *LDA,double *WR,double *WI,
                 double *VL,int *LDVL,double *VR,int *LDVR,int *ILO,int *IHI,
                 double *SCALE,double *ABNRM,
                 double *RCONDE,double *RCONDV,double *WORK,int *LWORK,
                 int *IWORK,int *INFO);
 
   void FORTRAN_NAME(dgesv)(int *LDA, int *NRHS, double *A, int *LDA2, int *IPIV,
                 double *B, int *LDB, int *INFO);
}


void setup (int N, const Parameter &param, Array<double, 1> &WR, Array<double,2> &ev, Array<double,2> &evInv)
{
   int Nm1 = N;
   int i;
   Array<double, 1> x;
   Array<double, 2> D;
   Array<double, 1> r;
   Array<double, 2> Dsec;
   Array<double, 1> XX;
   Array<double, 1> YY;
   Array<double, 2> A(N,N);
   Array<double, 2> B(N,N);
   Array<int, 1> IPIV(Nm1);
    char BALANC[1];
    char JOBVL[1];
    char JOBVR[1];
    char SENSE[1];
    int LDA;
    int LDVL;
    int LDVR;
    int NRHS;
    int LDB;
    int INFO;

	//resize output arrays
	WR.resize(N);
	ev.resize(N, N);
	evInv.resize(N, N);


// parameters for DGEEVX
    Array<double, 1>  WI(Nm1); // WR(Nm1),
           // The real and imaginary part of the eig.values
	Array<double, 2> VL(N, N);
    Array<double, 2> VR(Nm1,Nm1); //VR(Nm1,Nm1); 
           // The left and rigth eigenvectors
    int ILO, IHI;        // Info on the balanced output matrix
    Array<double, 1> SCALE(Nm1);     // Scaling factors applied for balancing
    double ABNRM;        // 1-Norm of the balanced matrix
    Array<double, 1> RCONDE(Nm1);  
           // the reciprocal cond. numb of the respective eig.val
    Array<double, 1> RCONDV(Nm1); 
           // the reciprocal cond. numb of the respective eig.vec
    int LWORK = (N+1)*(N+7); // Depending on SENSE            
    Array<double, 1> WORK(LWORK);
    Array<int, 1> IWORK(2*(N+1)-2);


// Compute the Chebyshev differensiation matrix and D*D

//   cheb(N, x, D);
   cheb(N, x, D);
   Dsec.resize(D.shape());
   MatrixMatrixMultiply(D, D, Dsec);


// Compute the 1. and 2. derivatives of the transformations

   XYmat(N, param, XX, YY, r);

   // Set up the full timepropagation matrix A
   // dy/dt = - i A y
   Range range(1, N); //Dsec and D have range 0, N+1. 
   					  //We don't want the edge points in A
   A = XX(tensor::i) * Dsec(range, range) + YY(tensor::i) * D(range, range);
   //Transpose A
   for (int i=0; i<A.extent(0); i++)
   {
		for (int j=0; j<i; j++)
		{
			double t = A(i,j);
			A(i,j) = A(j, i);
			A(j,i) = t;
		}
   }

// Add radialpart of non-time dependent potential here

// Compute eigen decomposition

   BALANC[0] ='B';
   JOBVL[0]  ='V';
   JOBVR[0]  ='V';
   SENSE[0]  ='B';
   LDA = Nm1;
   LDVL = Nm1;
   LDVR = Nm1;

   FORTRAN_NAME(dgeevx)(BALANC, JOBVL, JOBVR, SENSE, &Nm1,
            A.data(), &LDA, WR.data(), WI.data(),
            VL.data(), &LDVL, VR.data(), &LDVR, &ILO, &IHI,
            SCALE.data(), &ABNRM,
            RCONDE.data(), RCONDV.data(), WORK.data(), &LWORK,
            IWORK.data(), &INFO);

// Compute the inverse of the eigen vector matrix

    NRHS = Nm1;
	
    evInv = VR ;// VL;
    LDB = LDA;
    B = 0.0;
    for (i=0; i<Nm1; i++) B(i,i) = 1.0; 

    FORTRAN_NAME(dgesv)(&Nm1, &NRHS, evInv.data(), &LDA, IPIV.data(), B.data(), &LDB, &INFO);

	ev = VR(tensor::j, tensor::i);   //Transpose
	evInv = B(tensor::j, tensor::i); //Transpose
	
	//cout << "Eigenvectors (right): " << ev << endl;
	//cout << "Eigenvectors (inv): " << evInv << endl;
    //printf(" Done  inverse, INFO = %d \n", INFO);
} // done

};

