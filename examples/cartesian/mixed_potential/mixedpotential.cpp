#include <core/wavefunction.h>
#include <core/representation/cartesianrepresentation.h>
#include <core/potential/dynamicpotentialevaluator.h>

/* 
 * Example potential for the MixedPotential momentum 
 * evaluator. One rank (specified in the config
 * parameter fourierrank) is transformed into fourier space
 * whereas all the other ranks are in grid space
 */


/*
 * Base for all potentials in this project. It contains the ApplyConfigSection method, 
 * which ensures that all potentials get the same parameters
 */
class PotentialBase
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double omega;
	double A0;
	double q1,q2,m1,m2;

	double static const c = 137.0;

	//Calculated pararmeters
	double Q;
	double qTilde;
	double qMerket;
	double my;
	double M;
	double E0;	

	virtual void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("omega", omega);
		config.Get("a0", A0);
		config.Get("q1", q1);
		config.Get("q2", q2);
		config.Get("m1", m1);
		config.Get("m2", m2);	

		Q = q1 + q2;
		M = m1 + m2;
		my = m1 * m2 / (m1 + m2);
		qTilde = my * (q1/m1 - q2/m2);
		qMerket = (q1*m2*m2 + q2*m1*m1) / (M*M);
		E0 = A0 * omega;
	}
};
double const PotentialBase::c;


/*
 * MixedPotentialX is the mixed potential where the x-rank (rank 0) is transformed
 * into fourier space 
 */
template<int Rank>
class MixedPotentialX : public PotentialBase
{
public:
	virtual void ApplyConfigSection(const ConfigSection &config)
	{
		PotentialBase::ApplyConfigSection(config);
	}

	double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		//Kinetic energy
		double p_x = pos(0);

		double kineticEnergy = p_x * p_x / (2.0 * M);
		
		return kineticEnergy; 
	}
};


/*
 * MixedPotentialY is the mixed potential where the y-rank (rank 1) is transformed
 * into fourier space, while the other dimensions are in grid space
 */
template<int Rank>
class MixedPotentialY : public PotentialBase
{
public:
	void ApplyConfigSection(const ConfigSection &config)
	{
		PotentialBase::ApplyConfigSection(config);
	}

	cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		//Kinetic energy
		double x = pos(0);
		double p_y = pos(1);

		double kineticEnergy = p_y * p_y / (2.0 * M);
		cplx fieldEnergy  = I * qTilde/my * A0 * sin(omega * CurTime) * p_y;
		cplx fieldEnergy2 = qMerket / (my * c) * I * p_y * E0 * cos(omega * CurTime) * x;
		
		return kineticEnergy + fieldEnergy + fieldEnergy2;
	}
};

/*
 * GridPotential is the potential where all ranks are in grid space
 */
template<int Rank>
class GridPotential : public PotentialBase
{
public:
	double SoftRadius;

	void ApplyConfigSection(const ConfigSection &config)
	{
		PotentialBase::ApplyConfigSection(config);
		config.Get("soft_radius", SoftRadius);
	}

	cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double x = pos(0);
		double y = pos(1);
		double r = sqrt(x*x + y*y + SoftRadius*SoftRadius);
		//r = max(r, minimumR);
		
		double coloumbEnergy = q1 * q2 / r ;
		
		double fieldStrength = A0 * E0 * sin(omega * CurTime) * cos(omega * CurTime) * x;
		double fieldEnergy = (1.0/c) * qTilde * (Q/M + qMerket/my) * fieldStrength;
		
		return coloumbEnergy + fieldEnergy;
	}
};


