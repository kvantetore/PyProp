#ifndef SPHERICALHARMONICS_H 
#define SPHERICALHARMONICS_H

#include <blitz/array.h>

#include "common.h"

enum SphericalHarmonicsRank
{
	SphericalHarmonicsRankLm = 0,
	SphericalHarmonicsRankAngular = 1
};

enum AngularValuesRank
{
	AngularValuesRankIndex = 0,
	AngularValuesRankCoord = 1
};

enum AngularCoordinateIndex
{
	AngularCoordIndexTheta = 0,
	AngularCoordIndexPhi = 1
};

/** 
@author = Raymon Nepstad / Tore Birkeland
@remarks = Creates the associated legendre functions P_{l,m}(cos theta) for all m and l up to l = maxL. 
A two dimensional grid is returned, where 
- the first rank is the {l,m}-index: lm = (m + 2*l) + (l*l)
- the second rank is the angular index corresonding to the input array ThetaValues 
*/
blitz::Array<cplx, 2> CreateAssociatedLegendre(int maxL, blitz::Array<double, 1> &ThetaGrid);

/**
@author = Tore Birkeland
@remarks = Creates the spherical harmonics function Y_{l,m)(theta, phi) for all m and l up to l = maxL.
A two dimensional grid is where 
- the first rank is the {l,m}-index: lm = (m + 2*l) + (l*l)
- the second rank is the angular index corresonding to the input array Angular values, 
  Omega_j = {Theta_j, Phi_j}
  
Use the functions MapYlmIndex() and MapYlmReverse() to map to and from {l,m}-indices and
omega-indices. 
*/
blitz::Array<cplx, 2> CreateSphericalHarmonics(int maxL, blitz::Array<double, 2> &AngularGrid);

/**
@author = Tore Birkeland
@remarks = Creates an equi-spaced grid Omega_{j*maxL + i} = {Theta_i, Phi_j} 
for i, j: {0 -> maxL}
theta : {0 -> pi}
phi   : {0 -> 2 pi}
*/
blitz::Array<double, 2> CreateEquispacedAngularGrid(int maxL);

/* 
	Mapping to and from linear ylm-index and (l,m)-index 
*/
// Get linear ylm index from (l,m)
inline int MapYlmIndex(int l, int m)
{
	return (m + 2*l) + (l*l);
}

//Get the l-value from a linear ylm index
inline int GetYlmL(int index)
{
	return static_cast<int>(sqrt(index));
}

//Get (l,m) from linear ylm index
inline blitz::TinyVector<int,2> MapYlmReverse(int index)
{
	int l = GetYlmL(index);
	int m = index - 2 * l - l * l;
	return blitz::TinyVector<int, 2>(l, m);
}

//Get the m-value from a linear ylm index
inline int GetYlmM(int index)
{
	return MapYlmReverse(index)(1);
}

inline cplx AsImaginary(double value)
{
	return cplx(0.0, value);
}

//Make the index-functions available in blitz expressions
BZ_DECLARE_FUNCTION(GetYlmM);
BZ_DECLARE_FUNCTION_RET(AsImaginary, cplx);

#endif
