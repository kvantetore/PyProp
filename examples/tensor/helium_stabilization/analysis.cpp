#include <core/common.h>
#include <core/wavefunction.h>
#include <core/representation/coupledspherical/coupledrange.h>
#include <core/transform/spherical/shtools.h>

#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_coulomb.h>
#include <gsl/gsl_errno.h>


using namespace boost::python;
using namespace blitz;

typedef Array<cplx, 2> MatrixType;
typedef Array<cplx, 1> VectorType;

/*
 * Calculate the projection of the wavefunction on a product of single particle
 * radial wavefunctions for corresponding to a given angular momentum for 
 * each radial wavefunction. 
 *
 *	<r1(1) r2(2) | psi(1,2) >
 *
 * It takes in a list of angular indices, corresponding to which
 * angular momentum indices is corresponding to this l1, l2 pair.
 *
 * Remarks:
 * - No symmetrization is made on either psi or the radial functions
 * - Integration weights are assumed to have been applied to psi on beforehand
 *
 * A 3D array is returned
 * 		rank0: angular indices
 * 		rank1: radial function 1 indices
 * 		rank2: radial function 2 indices
 */
Array<cplx, 3> CalculateProjectionRadialProductStates(int l1, MatrixType V1, int l2, MatrixType V2, Array<cplx, 3> psiData, Array<int, 1> angularIndices)
{
	int count0 = angularIndices.extent(0);
	int count1 = V1.extent(1);
	int count2 = V2.extent(1);
	int rcount = V1.extent(0);

	blitz::Array<cplx, 3> proj(count0, count1, count2);
	blitz::Array<cplx, 2> tempProj(rcount, count2);

	for (int i0=0; i0<angularIndices.extent(0); i0++)
	{
		int angIdx = angularIndices(i0);
		MatrixType psiSlice = psiData(angIdx, Range::all(), Range::all());


		//project on V2 states
		tempProj = 0;
		for (int r1=0; r1<rcount; r1++)
		{
			for (int r2=0; r2<rcount; r2++)
			{
				cplx curPsi = psiSlice(r1, r2);
				for (int i2=0; i2<count2; i2++)
				{
					tempProj(r1, i2) += conj(V2(r2, i2)) * curPsi;
				}
			}
		}

		//project on V1 states
		for (int r1=0; r1<rcount; r1++)
		{
			for (int i1=0; i1<count1; i1++)
			{
				for (int i2=0; i2<count2; i2++)
				{
					proj(i0, i1, i2) += conj(V1(r1, i1)) * tempProj(r1, i2);
				}
			}
		}
	}

	return proj;
}


/*
 * Removes the projection of psi on the radial product states in V1, V2
 * The projection is calculated from psiData, and is removed from outputPsi
 *
 * integration weights should be multiplied to psiData, but not to outputPsi, 
 * see CalculateProjectionRadialProductStates above
 */
void RemoveProjectionRadialProductStates(int l1, MatrixType V1, int l2, MatrixType V2, Array<cplx, 3> psiData, Array<int, 1> angularIndices, Array<cplx, 3> outputPsi)
{
	int count0 = angularIndices.extent(0);
	int count1 = V1.extent(1);
	int count2 = V2.extent(1);
	int rcount = V2.extent(0);

	//Calculate projection
	Array<cplx, 3> proj = CalculateProjectionRadialProductStates(l1, V1, l2, V2, psiData, angularIndices);

	for (int i0=0; i0<count0; i0++)
	{
		int angIdx = angularIndices(i0);
		MatrixType psiSlice = outputPsi(angIdx, Range::all(), Range::all());

		for (int r1=0; r1<rcount; r1++)
		{
			for (int r2=0; r2<rcount; r2++)
			{
				for (int i1=0; i1<count1; i1++)
				{
					for (int i2=0; i2<count2; i2++)
					{
						psiSlice(r1, r2) -= proj(i0, i1, i2) * V1(r1, i1) * V2(r2, i2);
					}
				}
			}
		}
	}
}


/* 
 * Calculates the population (abssqr of projection) of the wavefunction
 * on a product of radial states
 *
 * see CalculateProjectionRadialProductStates.
 *
 * Returns a list of tuples (i1, i2, pop) giving the population
 * in state i1, i2 for the given (l1, l2) pair. summed over angularIndices
 */
list CalculatePopulationRadialProductStates(int l1, MatrixType V1, int l2, MatrixType V2, Array<cplx, 3> psiData, Array<int, 1> angularIndices)
{
	int count1 = V1.extent(1);
	int count2 = V2.extent(1);

	Array<cplx, 3> proj = CalculateProjectionRadialProductStates(l1, V1, l2, V2, psiData, angularIndices);
	proj *= conj(proj);
	
	list popList;
	for (int i1=0; i1<count1; i1++)
	{
		for (int i2=0; i2<count2; i2++)
		{
			double pop = 2 * real(sum(proj(Range::all(), i1, i2)));
			popList.append(make_tuple(i1, i2, pop));
		}
	}
	return popList;
}


// ------------------------------------------------------------------------------------------------------- 
//                                            Coulomb Wavefunction
// ------------------------------------------------------------------------------------------------------- 

/* 
 * Gets the Coulomb phase sigma_l = arg(gamma(l + 1 + i*eta))
 */
double GetCoulombPhase(int l, double eta)
{
	gsl_sf_result absval, argval;
	if (gsl_sf_lngamma_complex_e(1.0+l, eta, &absval, &argval) == GSL_ELOSS)
	{
		cout << "Overflow error in gsl_sf_lngamma_complex_e, l=" << l << ", eta=" << eta << endl;
	}
	return argval.val;
}

/*
 * Sets the radial Coulomb wave F_l(k*r, eta), with eta = Z/k into data for all radial
 * grid points specified by r
 */
void SetRadialCoulombWave(int Z, int l, double k, blitz::Array<double, 1> r, blitz::Array<double, 1> data)
{
	double eta = Z / k;

	for (int i=0; i<r.size(); i++)
	{	
		double x = k * r(i);
		gsl_sf_result F, Fp, G, Gp;
		double exp_F, exp_G;
		int error = gsl_sf_coulomb_wave_FG_e(eta, x, (double)l, 0., &F, &Fp, &G, &Gp, &exp_F, &exp_G);
		if (error == GSL_EOVRFLW)
		{
			cout << "WARNING: Overflow in SetCoulombWave(" << Z << ", " << l << ", " << k << ", r=" << r(i) << ");" << endl;
			cout << "         exp_F = " << exp_F << ", exp_G = " << exp_G << endl;
		}

		data(i) = F.val;
	}
}


// ------------------------------------------------------------------------------------------------------- 
//                                            Angular Distribution
// ------------------------------------------------------------------------------------------------------- 

/*
 * adds exp(-i pi / 2 (l1 + l2) + i (sigma_l1 + sigma_l2)) Y^{L,M}_{l1,l2}(theta1, theta2) radialProj(l1, l2, L, M, i1, i2) 
 * to the existing angularData(theta1, theta2, i1, i2)
 *
 * radialProj(l1, l2, L, M, E1, E2) should calculated from CalculateProjectionRadialProductStates
 *
 * angularData
 * 		rank0 - theta1
 * 		rank1 - theta2
 * 		rank2 - i1
 * 		rank3 - i2
 */
void AddDoubleAngularProjectionAvgPhi(Array<cplx, 4> angularData, Array<cplx, 2> sphericalHarmonics, 
		Array<cplx, 2> radialProj, int l1, int l2, int m, int M)
{
	//create a view to the current spherical harmonic and radialProjection arrays
	Array<cplx, 1> sph1 = sphericalHarmonics( MapLmIndex(l1, m), Range::all());
	Array<cplx, 1> sph2 = sphericalHarmonics( MapLmIndex(l2, M-m), Range::all());
	
	//Add data
	angularData +=    sph1(tensor::i) * sph2(tensor::j)
					* radialProj(tensor::k, tensor::l);
}


void AddSingleAngularProjectionAvgPhi(Array<cplx, 3> angularData, Array<cplx, 2> sphericalHarmonics,
		Array<cplx, 2> radialProj, int l2, int m)
{
	//create a view to the current spherical harmonic and radialProjection arrays
	Array<cplx, 1> sph2 = sphericalHarmonics( MapLmIndex(l2, m), Range::all());
	
	//Add data
	angularData +=    sph2(tensor::j)
					* radialProj(tensor::i, tensor::k);
}

// ------------------------------------------------------------------------------------------------------- 
//                                            Symmetrization
// ------------------------------------------------------------------------------------------------------- 

/*
 * Maps the given wavefunction to one where the particles are exchanged
 * psi(1,2) -> psi(2,1)
 */
Wavefunction<3>::Ptr GetWavefunctionParticleExchange(Wavefunction<3>::Ptr psi, list angularSymmetrizationPairs)
{
	typedef Array<cplx, 3> ArrayType;
	ArrayType data = psi->GetData();

	int countr = data.extent(1);
	typedef stl_input_iterator<tuple> Iterator;
	Iterator begin(angularSymmetrizationPairs);
	Iterator end;

	Wavefunction<3>::Ptr exchgPsi = psi->Copy();
	ArrayType exchgData = exchgPsi->GetData();

	for (Iterator i=begin; i!=end; i++)
	{
		int a1 = extract<int>((*i)[0]);
		int a2 = extract<int>((*i)[1]);

		for (int r1=0; r1<countr; r1++)
		{
			for (int r2=0; r2<countr; r2++)
			{
				exchgData(a1, r1, r2) = data(a2, r2, r1);
			}
		}
	}
	
	return exchgPsi;
}


