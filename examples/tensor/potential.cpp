#include <core/wavefunction.h>
#include <core/potential/dynamicpotentialevaluator.h>

#include <core/transform/bspline/bspline.h>
#include <core/utility/blitztricks.h>
#include <core/utility/blitzblas.h>

#include <core/transform/reducedspherical/reducedsphericaltools.h>

using namespace blitz;

template<int Rank>
cplx VectorDotProduct(Wavefunction<Rank> &psi1, Wavefunction<Rank> &psi2)
{
	return VectorInnerProduct(psi1.GetData(), psi2.GetData());
}

template<class TBase, int Rank>
void RepresentPotentialInBasisBSpline( BSpline::BSpline::Ptr bsplineObject, Array<TBase, Rank> source, Array<TBase, Rank> dest, Array<int, 2> indexPair, std::string storageId, int rank, int differentiation )
{
	typedef Array<TBase, 3> Array3D;
	typedef Array<TBase, 1> Array1D;
	Array3D source3d = MapToRank3(source, rank, 1);
	Array3D dest3d = MapToRank3(dest, rank, 1);

	int preCount = source3d.extent(0);
	int postCount = source3d.extent(2);
	int pairCount = indexPair.extent(0);

	bool isHermitian = storageId == "Herm";
	double scaling = 1;

	for (int preIndex=0; preIndex<preCount; preIndex++)
	{
		for (int pairIndex=0; pairIndex<pairCount; pairIndex++)
		{
			int rowIndex = indexPair(pairIndex, 0);
			int colIndex = indexPair(pairIndex, 1);
			if (isHermitian && rowIndex == colIndex)
			{
				scaling = 0.5;
			}
			else
			{
				scaling = 1.0;
			}

			for (int postIndex=0; postIndex<postCount; postIndex++)
			{
				Array1D sourceSlice = source3d(preIndex, Range::all(), postIndex);
				dest3d(preIndex, pairIndex, postIndex) = scaling * bsplineObject->BSplineGlobalOverlapIntegral(sourceSlice, differentiation, rowIndex, colIndex);
			}
		}
	}
}

typedef ReducedSpherical::ReducedSphericalTools::Ptr ReducedSphericalToolsPtr;

template<class TBase, int Rank>
void RepresentPotentialInBasisReducedSphericalHarmonic( ReducedSphericalToolsPtr obj, Array<TBase, Rank> source, Array<TBase, Rank> dest, Array<int, 2> indexPair, std::string storageId, int rank, int differentiation )
{
	typedef Array<TBase, 3> Array3D;
	typedef Array<TBase, 1> Array1D;
	Array3D source3d = MapToRank3(source, rank, 1);
	Array3D dest3d = MapToRank3(dest, rank, 1);

	int preCount = source3d.extent(0);
	int thetaCount = source3d.extent(1);
	int postCount = source3d.extent(2);
	int pairCount = indexPair.extent(0);

	
	Array<double, 2> leftAssocLegendre = obj->GetAssociatedLegendrePolynomial();
	Array<double, 2> rightAssocLegendre;
	if (differentiation == 0)
	{
		rightAssocLegendre.reference(obj->GetAssociatedLegendrePolynomial());
	}
	else if (differentiation == 1)
	{
		rightAssocLegendre.reference(obj->GetAssociatedLegendrePolynomialDerivative());
	}
	else 
	{
		throw std::runtime_error("Reduced spherical harmonics does not support higher derivatives than 1");
	}
	Array<double, 1> weights = obj->GetWeights();
	
	dest3d = 0;
	bool isHermitian = storageId == "Herm";
	double scaling = 1;

	for (int preIndex=0; preIndex<preCount; preIndex++)
	{
		for (int pairIndex=0; pairIndex<pairCount; pairIndex++)
		{
			int rowIndex = indexPair(pairIndex, 0);
			int colIndex = indexPair(pairIndex, 1);
			if (isHermitian && rowIndex == colIndex)
			{
				scaling = 0.5;
			}
			else
			{
				scaling = 1.0;
			}

			for (int thetaIndex=0; thetaIndex<thetaCount; thetaIndex++)
			{
				double leftLegendre = leftAssocLegendre(thetaIndex, rowIndex);
				double rightLegendre = rightAssocLegendre(thetaIndex, colIndex);
				double weight = weights(thetaIndex);
		
				for (int postIndex=0; postIndex<postCount; postIndex++)
				{
					dest3d(preIndex, pairIndex, postIndex) += leftLegendre * source3d(preIndex, thetaIndex, postIndex) * rightLegendre * weight;
				}
			}
		}
	}
}


template<class TBase>
void MultiplyTensorPotential_1(blitz::Array<TBase, 1> potentialData, double scaling, blitz::Array<int, 2> pair0, blitz::Array<cplx, 1> source, blitz::Array<cplx, 1> dest, int algo)
{
	for (int i=0; i<potentialData.extent(0); i++)
	{
		int row = pair0(i, 0);
		int col = pair0(i, 1);
		dest(row) += potentialData(i) * scaling * source(col);
	}
}

template<class TBase>
void MultiplyTensorPotential_2(blitz::Array<TBase, 2> potentialData, double scaling, blitz::Array<int, 2> pair0, blitz::Array<int, 2> pair1, blitz::Array<cplx, 2> source, blitz::Array<cplx, 2> dest, int algo)
{
	if (algo == 0)
	{
		for (int i=0; i<potentialData.extent(0); i++)
		{
			int row0 = pair0(i, 0);
			int col0 = pair0(i, 1);
			for (int j=0; j<potentialData.extent(1); j++)
			{
				int row1 = pair1(j, 0);
				int col1 = pair1(j, 1);
		
				dest(row0, row1) += potentialData(i, j) * scaling * source(col0, col1);
			}
		}
	}

	if (algo == 1)
	{
		int k = 8;
		
		for (int i=0; i<potentialData.extent(0); i++)
		{
			int row0 = pair0(i, 0);
			int col0 = pair0(i, 1);
		
			int N1 = source.extent(1);

			cplx* __restrict__ potentialPtr = &potentialData(i, 0);
			cplx* __restrict__ destPtr = &dest(row0, 0);
			for (int row1=0; row1<N1; row1++)
			{
				int startCol1 = std::max(0, row1-k+1);
				int endCol1 = std::min(N1, row1+k);
				cplx* __restrict__ sourcePtr = &source(col0, startCol1);
				for (int col1=startCol1; col1<endCol1; col1++)
				{
					(*destPtr) += (*potentialPtr++) * scaling * (*sourcePtr++);
				}
				destPtr++;
			}
		}
	}

	if (algo == 2)
	{
		int k = 8;
		
		for (int i=0; i<potentialData.extent(0); i++)
		{
			int row0 = pair0(i, 0);
			int col0 = pair0(i, 1);
		
			int N1 = source.extent(1);

			cplx* __restrict__ potentialPtr = &potentialData(i, 0);
			cplx* __restrict__ destPtr = &dest(row0, 0);
			
			cplx* __restrict__ sourceStartPtr = &source(col0, 0);
			for (int row1=0; row1<k; row1++)
			{
				int startCol1 = 0;
				int endCol1 = row1+k;
				//cout << "Start -> Stop (0)= " << startCol1 << " -> " << endCol1 << endl;
				cplx* __restrict__ sourcePtr = sourceStartPtr;
				for (int col1=startCol1; col1<endCol1; col1++)
				{
					(*destPtr) += (*potentialPtr++) * scaling * (*sourcePtr++);
				}
				destPtr++;
			}
			for (int row1=k; row1<N1-k; row1++)
			{
				int startCol1 = row1-k+1;
				int endCol1 = row1+k;
				//cout << "Start -> Stop (1)= " << startCol1 << " -> " << endCol1 << endl;
				cplx* __restrict__ sourcePtr = sourceStartPtr;
				sourceStartPtr++;
				for (int col1=startCol1; col1<endCol1; col1++)
				{
					(*destPtr) += (*potentialPtr++) * scaling * (*sourcePtr++);
				}
				destPtr++;
			}
			for (int row1=N1-k; row1<N1; row1++)
			{
				int startCol1 = row1-k+1;
				int endCol1 = N1;
				//cout << "Start -> Stop (2)= " << startCol1 << " -> " << endCol1 << endl;
				cplx* __restrict__ sourcePtr = sourceStartPtr;
				sourceStartPtr++;
				for (int col1=startCol1; col1<endCol1; col1++)
				{
					(*destPtr) += (*potentialPtr++) * scaling * (*sourcePtr++);
				}
				destPtr++;
			}
			//exit(-1);
		}
	}

/*
	if (algo == 3)
	{
		int k = 8;
		
		for (int i=0; i<potentialData.extent(0); i++)
		{
			int row0 = pair0(i, 0);
			int col0 = pair0(i, 1);
		
			int N1 = source.extent(1);

			cplx* __restrict__ potentialPtr = &potentialData(i, 0);
			cplx* __restrict__ destPtr = &dest(row0, 0);
			for (int band1=-k+1; band1<k; band1++)
			{
				int start1 = 0;
				int end1 = N1;
				if (band1 < 0)
				{
					start1 = -band1;
				}
				else
				{
					end1 = N1 - band1;
				}

				int j = band1+k-1;
				for (int index1=start1; index1<end1; index1++)
				{
					int col1 = index1+band1+k;
					int row1 = index1;
					if (row1 != pair1(j, 0) || col1 != pair1(j, 1))
					{
						cout << "Error (" << row1 << ", " << col1 << ") != (" << pair1(j, 0) << ", " <<  pair1(j, 1) << ")" << endl; 
					}
					dest(row0, row1) += potentialData(i, j) * source(col0, col1);
					j += 2*k - 1;
				}
			}
			exit(-1);
		}
	}
*/
}



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
class DipoleLaserPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	int angularRank;
	int radialRank;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("angular_rank", angularRank);
		config.Get("radial_rank", radialRank);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double theta = pos(angularRank);
		double r = pos(radialRank);
		return r * cos(theta);
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

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("angular_rank", angularRank);
		config.Get("radial_rank", radialRank);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double theta = pos(angularRank);
		return cplx(0.,1.)*cos(theta);
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

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("angular_rank", angularRank);
		config.Get("radial_rank", radialRank);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r = pos(radialRank);
		double theta = pos(angularRank);
		//return -cplx(0.,1.)*sin(theta) / r;
		return -cplx(0.,1.) / r;
		//return cos(theta) / r;
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

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("angular_rank", angularRank);
		config.Get("radial_rank", radialRank);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r = pos(radialRank);
		double theta = pos(angularRank);
		return - cplx(0.,1.) * cos(theta) / r;
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
class RadialHarmonicPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	int angularRank;
	int radialRank;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("angular_rank", angularRank);
		config.Get("radial_rank", radialRank);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r = pos(radialRank);
		return 0.5*r*r;
	}
};

template<int Rank>
class RadialCoulombPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	int angularRank;
	int radialRank;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("angular_rank", angularRank);
		config.Get("radial_rank", radialRank);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r = pos(radialRank);
		return - 1.0 / r;
	}
};



