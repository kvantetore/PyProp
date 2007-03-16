#include <core/wavefunction.h>
#include <core/potential/sphericaldynamicpotentialevaluator.h>

template<int Rank>
class H2pPotential
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double Charge;
	double Softing;
	double cosOrient;
	double sinOrient;

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("charge", Charge);
		config.Get("softing", Softing);
		
		double orientation = 0;
		config.Get("nuclear_orientation", orientation);
		sinOrient = sin(orientation);
		cosOrient = cos(orientation);

		cout << "Orient " << cosOrient << " " << sinOrient << endl;
	}

	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		//Coordinates
		double R = fabs(pos(0));
		double r = fabs(pos(1));
		double theta = pos(2);
		double phi = pos(3);

		double r2 = sqr(r) + sqr(R) / 4 + sqr(Softing);
		double z = r * cos(theta);
		double x = r * sin(theta) * cos(phi);

		//cout << "R, r, theta, phi = "<< R << ", " << r << ", " << theta << ", " << phi << endl;
	
		double angDep = R * (x * sinOrient + z * cosOrient); 
		double V1 = 1 / sqrt(r2 + angDep); 
		double V2 = 1 / sqrt(r2 - angDep); 

		double V3 = 1.0 / sqrt(sqr(R) + sqr(Softing));

		return Charge * (V1 + V2) + V3;
	}
};




