#include <core/wavefunction.h>
#include <core/potential/dynamicpotentialevaluator.h>

using namespace blitz;

template<int Rank>
class KineticEnergyPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double mass;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("mass", mass);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		return - 1. / (2. * mass);
	}

};


template<int Rank>
class DipoleLaserPotentialLength : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	int angularRank;
	int radialRank;
	double Charge;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("angular_rank", angularRank);
		config.Get("radial_rank", radialRank);
		config.Get("charge", Charge);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double theta = pos(angularRank);
		double r = pos(radialRank);
		return -Charge * r * cos(theta);
	}
};

template<int Rank>
class DipoleLaserPotentialVelocityRadialDerivative : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	int angularRank;
	int radialRank;
	double Charge;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("angular_rank", angularRank);
		config.Get("radial_rank", radialRank);
		config.Get("charge", Charge);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double theta = pos(angularRank);
		return Charge * cplx(0.,1.)*cos(theta);
	}
};

template<int Rank>
class DipoleLaserPotentialVelocityAngularDerivative : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	int angularRank;
	int radialRank;
	double Charge;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("angular_rank", angularRank);
		config.Get("radial_rank", radialRank);
		config.Get("charge", Charge);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r = pos(radialRank);
		//double theta = pos(angularRank);
		//return -cplx(0.,1.)*sin(theta) / r;
		return -Charge * cplx(0.,1.) / r;
	}
};

template<int Rank>
class DipoleLaserPotentialVelocity: public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	int angularRank;
	int radialRank;
	double Charge;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("angular_rank", angularRank);
		config.Get("radial_rank", radialRank);
		config.Get("charge", Charge);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r = pos(radialRank);
		double theta = pos(angularRank);
		return -Charge * cplx(0.,1.) * cos(theta) / r;
	}
};

template<int Rank>
class AngularKineticEnergyPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	int angularRank;
	int radialRank;
	double mass;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("angular_rank", angularRank);
		config.Get("radial_rank", radialRank);
		config.Get("mass", mass);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double l = pos(angularRank);
		double r = pos(radialRank);
		return l*(l+1.0) / (2 * mass * r*r);
	}
};


template<int Rank>
class SingleActiveElectronPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	int radialRank;

	//Potential parameters
	double Z;
	double a1;
	double a2;
	double a3;
	double a4;
	double a5;
	double a6;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("z", Z);
		config.Get("a1", a1);
		config.Get("a2", a2);
		config.Get("a3", a3);
		config.Get("a4", a4);
		config.Get("a5", a5);
		config.Get("a6", a6);

		config.Get("radial_rank", radialRank);
	}

	/*
	 * Called for every grid point at every time step. 
	 *
	 * Some general tips for max efficiency:
	 * - If possible, move static computations to ApplyConfigSection.
	 * - Minimize the number of branches ("if"-statements are bad)
	 * - Minimize the number of function calls (sin, cos, exp, are bad)
	 * - Long statements can confuse the compiler, consider making more 
	 *   simpler statements
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r = std::fabs(pos(radialRank));
		
		double V = -(Z + a1 * exp(-a2 * r) + a3 * r * exp(-a4 * r)
			+ a5 * exp(-a6 * r)) / r;

		return V;
	}
};


template<int Rank>
class ComplexAbsorbingPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	int radialRank;
	double scalingReal;
	double scalingImag;
	double factorReal;
	double factorImag;
	double absorberStart;
	double absorberLength;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("radial_rank", radialRank);
		config.Get("absorber_start", absorberStart);
		config.Get("absorber_length", absorberLength);
		config.Get("scaling_real", scalingReal);
		config.Get("scaling_imag", scalingImag);
		config.Get("factor_real", factorReal);
		config.Get("factor_imag", factorImag);
	}

	/*
	 * Called for every grid point 
	 */
	inline cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r = pos(radialRank);
		cplx V = 0;
		if (absorberStart < r && r <= absorberStart+absorberLength)
		{
			double curLength = (r - absorberStart) / absorberLength;
			double Vr = factorReal * std::pow(curLength, scalingReal);
			double Vi = factorImag * std::pow(curLength, scalingImag);
			V = cplx(Vr , Vi);
		}

		return V;
	}
};


#include <core/representation/combinedrepresentation.h>
#include <core/representation/reducedspherical/reducedsphericalharmonicrepresentation.h>
#include "laserhelper.h"

/* First part of the linearly polarized laser in the velocity gauge
 * expressed in spherical harmonics
 *
 * <Ylm | - \frac{1}{r} \sin \theta \partialdiff{}{\theta} 
 *	      - \frac{\cos \theta}{r} | Yl'm'>
 */
template<int Rank>
class CustomPotential_LaserVelocity1_ReducedSpherical
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

private:
	BasisPairList AngularBasisPairs;
	int AngularRank;
	int RadialRank;

public:
	CustomPotential_LaserVelocity1_ReducedSpherical() {}
	virtual ~CustomPotential_LaserVelocity1_ReducedSpherical() {}

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("radial_rank", RadialRank);
		config.Get("angular_rank", AngularRank);
	}

	virtual void SetBasisPairs(int rank, const BasisPairList &basisPairs)
	{
		if (rank != AngularRank)
		{
			throw std::runtime_error("Only angular rank supports basis pairs");
		}
		AngularBasisPairs.reference(basisPairs.copy());
	}

	BasisPairList GetBasisPairList(int rank)
	{
		if (rank == AngularRank)
			return AngularBasisPairs;
		else
			return BasisPairList();
	}

	virtual void UpdatePotentialData(typename blitz::Array<cplx, Rank> data, typename Wavefunction<Rank>::Ptr psi, cplx timeStep, double curTime)
	{
		using namespace ReducedSpherical;

		typedef CombinedRepresentation<Rank> CmbRepr;
		typename CmbRepr::Ptr repr = boost::static_pointer_cast< CmbRepr >(psi->GetRepresentation());
		ReducedSphericalHarmonicRepresentation::Ptr angRepr = boost::static_pointer_cast< ReducedSphericalHarmonicRepresentation >(repr->GetRepresentation(AngularRank));
	
		int rCount = data.extent(RadialRank);
		int angCount = data.extent(AngularRank);

		blitz::Array<double, 1> localr = psi->GetRepresentation()->GetLocalGrid(RadialRank);
		BasisPairList angBasisPairs = GetBasisPairList(AngularRank);

		if (psi->GetRepresentation()->GetDistributedModel()->IsDistributedRank(AngularRank)) throw std::runtime_error("Angular rank can not be distributed");
		if (data.extent(RadialRank) != rCount) throw std::runtime_error("Invalid r size");
		if (data.extent(AngularRank) != angBasisPairs.extent(0)) throw std::runtime_error("Invalid ang size");

		cplx IM(0,1.0);

		data = 0;
		blitz::TinyVector<int, Rank> index;
	
		for (int angIndex=0; angIndex<angCount; angIndex++)
		{
			index(AngularRank) = angIndex;

			int leftIndex = angBasisPairs(angIndex, 0);
			int rightIndex = angBasisPairs(angIndex, 1);
	
			//"Left" quantum numbers
			int l = leftIndex;
			int m = 0;
			
			//"Right" quantum numbers (Mp = M)
			int lp = rightIndex;
			int mp = 0;

			double C = LaserHelper::C(lp, mp) * LaserHelper::kronecker(l, lp+1);
			double D = LaserHelper::D(lp, mp) * LaserHelper::kronecker(l, lp-1);
			double E = LaserHelper::E(lp, mp) * LaserHelper::kronecker(l, lp+1);
			double F = LaserHelper::F(lp, mp) * LaserHelper::kronecker(l, lp-1);

			double coupling = -(C + D) - (E + F);

			for (int ri=0; ri<rCount; ri++)
			{
				index(RadialRank) = ri;
				double r = localr(ri);

				data(index) = - IM * coupling/r;
			}
		}
	}
};

/* Second part of the linearly polarized laser in the velocity gauge
 * expressed in spherical harmonics.
 *
 * Should be used with first order differentiation in r
 *
 * <Ylm | \frac{\partial}{\partial r} \cos \theta | Yl'm'>
 */
template<int Rank>
class CustomPotential_LaserVelocity2_ReducedSpherical
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

private:
	BasisPairList AngularBasisPairs;
	int AngularRank;
	int RadialRank;

public:
	CustomPotential_LaserVelocity2_ReducedSpherical() {}
	virtual ~CustomPotential_LaserVelocity2_ReducedSpherical() {}

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("radial_rank", RadialRank);
		config.Get("angular_rank", AngularRank);
	}

	virtual void SetBasisPairs(int rank, const BasisPairList &basisPairs)
	{
		if (rank != AngularRank)
		{
			throw std::runtime_error("Only angular rank supports basis pairs");
		}
		AngularBasisPairs.reference(basisPairs.copy());
	}

	BasisPairList GetBasisPairList(int rank)
	{
		if (rank == AngularRank)
			return AngularBasisPairs;
		else
			return BasisPairList();
	}

	virtual void UpdatePotentialData(typename blitz::Array<cplx, Rank> data, typename Wavefunction<Rank>::Ptr psi, cplx timeStep, double curTime)
	{
		using namespace ReducedSpherical;

		typedef CombinedRepresentation<Rank> CmbRepr;
		typename CmbRepr::Ptr repr = boost::static_pointer_cast< CmbRepr >(psi->GetRepresentation());
		ReducedSphericalHarmonicRepresentation::Ptr angRepr = boost::static_pointer_cast< ReducedSphericalHarmonicRepresentation >(repr->GetRepresentation(AngularRank));
	
		int rCount = data.extent(RadialRank);
		int angCount = data.extent(AngularRank);

		blitz::Array<double, 1> localr = psi->GetRepresentation()->GetLocalGrid(RadialRank);
		BasisPairList angBasisPairs = GetBasisPairList(AngularRank);

		if (psi->GetRepresentation()->GetDistributedModel()->IsDistributedRank(AngularRank)) throw std::runtime_error("Angular rank can not be distributed");
		if (data.extent(RadialRank) != rCount) throw std::runtime_error("Invalid r size");
		if (data.extent(AngularRank) != angBasisPairs.extent(0)) 
		{
			cout << "Angular Rank = " << AngularRank << ", " << data.extent(AngularRank) << " != " << angBasisPairs.extent(0) << endl;
			throw std::runtime_error("Invalid ang size");
		}

		data = 0;
		blitz::TinyVector<int, Rank> index;

		cplx IM(0,1.0);
	
		for (int angIndex=0; angIndex<angCount; angIndex++)
		{
			index(AngularRank) = angIndex;

			int leftIndex = angBasisPairs(angIndex, 0);
			int rightIndex = angBasisPairs(angIndex, 1);
	
			//"Left" quantum numbers
			int l = leftIndex;
			int m = 0;
			
			//"Right" quantum numbers (Mp = M)
			int lp = rightIndex;
			int mp = 0;

			double E = LaserHelper::E(lp, mp) * LaserHelper::kronecker(l, lp+1);
			double F = LaserHelper::F(lp, mp) * LaserHelper::kronecker(l, lp-1);

			double coupling = (E + F);

			for (int ri=0; ri<rCount; ri++)
			{
				index(RadialRank) = ri;
				double r = localr(ri);

				data(index) =  - IM * coupling;
			}
		}
	}
};


/*
 * Radial Mask potential, used to filter out the part of the wavefunction outside a 
 * box in the radial grid
 */
template<int Rank>
class RadialMaskPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double CutoffDistance;
	double MaskStart;
	double MaskEnd;
	int RadialRank;


	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("radial_rank", RadialRank);
		config.Get("mask_start", MaskStart);
		config.Get("mask_end", MaskEnd);
	}

	/*
	 * Called for every grid point.
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r = pos(RadialRank);

		double V1;
		if ( MaskStart <= r && r < MaskEnd )
		{
			V1 = 1;
		}
		else
		{
			V1 = 0;
		}

		return V1;
	}
};


#include <gsl/gsl_sf_coulomb.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>

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

double GetCoulombNormalization(double Z, int l, double k)
{
	double eta = Z / k;
	gsl_sf_result C;
	int error = gsl_sf_coulomb_CL_e(l, eta, &C);
	return C.val;
}

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


