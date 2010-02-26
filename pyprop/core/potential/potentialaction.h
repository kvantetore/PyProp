#ifndef POTENTIAL_ACTION_H
#define POTENTIAL_ACTION_H

// Action classes define the operations to be performed on the data

// Apply the Exp of the Potential directly to the data
template <int Rank>
class ApplyPotentialClass
{
public:
	typedef typename blitz::Array<cplx,Rank>::iterator iterator;
		
	void ApplyAction(iterator &it, const cplx &potentialValue, const cplx &timeStep)
	{
		*it *= exp( - I * potentialValue * timeStep );
	}
};

//Multiply the the Potential directly to the data
template <int Rank>
class MultiplyPotentialClass
{
public:
	typedef typename blitz::Array<cplx,Rank>::iterator iterator;
	iterator DestIterator;
		
	void ApplyAction(iterator &it, const cplx &potentialValue, const cplx &timeStep)
	{
		*DestIterator += (*it) * potentialValue;
		DestIterator++;
	}
};

// Update the Potential data with the potential to be applied
template <int Rank>
class UpdatePotentialClass
{
public:
	typedef typename blitz::Array<cplx,Rank>::iterator iterator;
	
	void ApplyAction(iterator &it, const cplx &potentialValue, const cplx &timeStep)
	{
		*it = exp( - I * potentialValue * timeStep );
	}
};

// Get the potential values in the data
template <int Rank>
class GetPotentialClass
{
public:
	typedef typename blitz::Array<cplx,Rank>::iterator iterator;
	
	void ApplyAction(iterator &it, const cplx &potentialValue, const cplx &timeStep)
	{
		*it = potentialValue;
	}

};



#endif

