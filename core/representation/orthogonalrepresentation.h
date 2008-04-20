#ifndef ORTHOGONALREPRESENTATION_H
#define ORTHOGONALREPRESENTATION_H

#include <iostream>
#include "../common.h"
#include "representation.h"


/** A partial representation for orthogonal basises. It sets up the 
 * OverlapMatrix (which is diagonal) from the GetGlobalWeights()
 */
typedef Representation<1> Representation1D;

class OrthogonalRepresentation : public Representation1D
{
public:
	typedef shared_ptr<OrthogonalRepresentation> Ptr;

	//Constructors:
	OrthogonalRepresentation() {}
	virtual ~OrthogonalRepresentation() {}

	virtual std::complex<double> InnerProduct(const Wavefunction<1>& w1, const Wavefunction<1>& w2)
	{
		throw std::runtime_error("AngularRepresentation::InnerProduct is not implemented");
	}

	virtual OverlapMatrix::Ptr GetGlobalOverlapMatrix(int rank)
	{
		if (GlobalOverlapMatrix == 0)
		{
			OverlapMatrixEvaluatorDiagonal evaluator(this->GetGlobalWeights(rank));
			int size = this->GetFullShape()(0);
			GlobalOverlapMatrix = OverlapMatrix::Ptr(new OverlapMatrix(size, 0, evaluator));
		}
		return GlobalOverlapMatrix;
	}

private:
	OverlapMatrix::Ptr GlobalOverlapMatrix;

};

#endif



