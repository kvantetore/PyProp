#ifndef EXAMPLEPOTENTIALS_H
#define EXAMPLEPOTENTIALS_H

/* 
Example dynamic potential. for a Electric puslse convoluted by a cosine pulse
in the dipole approximation ( E = E0 * x *cos(conv(t)) * sin(w*t) ) , 
*/
template<int Rank>
class DipoleElectricPulse
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;
	
	//Potential parameters
	double FieldStrength;
	double Frequency;
	double ConvolutionLength;
	int PolarizationAxis;
	
	//Calculated parameters
	double ConvAngFreq;
	double AngularFrequency;

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("field_strength", FieldStrength);
		config.Get("frequency", Frequency);
		config.Get("convolution_time", ConvolutionLength);
		config.Get("polarization_axis", PolarizationAxis);
		
		ConvAngFreq = M_PI / ConvolutionLength;
		AngularFrequency = M_PI * Frequency;
	}

	inline cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double z = pos(PolarizationAxis);
		
		double convolution = cos(CurTime * ConvAngFreq);
		double field = sin(CurTime * M_PI * AngularFrequency);
		
		return FieldStrength * z * convolution * field;
	}
};

/*
 * Dynamic potential for a harmonic oscillator
 */
template<int Rank>
class HarmonicOscillatorPotential
{
public:
	double CurTime;
	cplx TimeStep;

	double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double pot = 0;
		for (int i=0;i<Rank;i++)
		{
			pot += sqr(pos(i));
		}
		return 0.5 * pot;
	}
};

/* 
Dynamic potential for evaluation of the kinetic energy potential for CartesianFFTEvaluator
*/
template<int Rank>
class CartesianKineticEnergyPotential
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;
	
	void ApplyConfigSection(const ConfigSection &config)
	{
	}

	inline cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &momentum)
	{
		double kineticPotential = 0.0;
		for (int i=0;i<Rank;i++)
		{
			kineticPotential += sqr(momentum(i));
		}
		return kineticPotential/2.0;
	}
};

#endif

