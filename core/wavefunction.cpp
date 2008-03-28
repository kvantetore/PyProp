#include "wavefunction.h"
#include "representation/representation.h"

template<int Rank>
void Wavefunction<Rank>::AllocateData()
{
	blitz::TinyVector<int, Rank> shape = Repr->GetInitialShape();

	int name = AllocateData(shape);
	SetActiveBuffer(name);
}

/* Returns the number of bytes currently allocated by this wavefunction */
template<int Rank>
size_t Wavefunction<Rank>::GetMemoryFootprint() const
{
	size_t size = 0;
	for (size_t i = 0; i < WavefunctionData.size(); i++)
	{
		size += WavefunctionData[i]->size();
	}
	return size * sizeof(cplx);
}

/* Allocates a new wavefunction buffer of size specifed by shape,
 * and returns the name of the data buffer
 */
template<int Rank>
int Wavefunction<Rank>::AllocateData(blitz::TinyVector<int, Rank> shape)
{
	unsigned long byteCount = (unsigned long)blitz::product(shape) * sizeof(cplx);
	if (byteCount/(1024*1024) > 2048)
	{
		cout << "Trying to allocate array larger than 2GB ("
		     << byteCount/(1024*1024*1024) << "GB). "
		     << "This is most likely an error, aborting." << endl;
		cout << "Shape = " << ToString(shape) << endl;
		throw std::runtime_error("Wavefunction too large");
	}

	std::cout 
		<< "Creating wavefunctions of shape " << shape
		<< " (~ " << byteCount / (1024*1024) << "MB)"
		<< std::endl;
	
	DataArrayPtr data = DataArrayPtr( new DataArray(shape) );
	int newName = WavefunctionData.size();
	WavefunctionData.push_back(data);
	return newName;
}

/* Frees a previously allocated data buffer specified by it name
 */
template<int Rank>
void Wavefunction<Rank>::FreeData(int bufferName) 
{
	(*WavefunctionData[bufferName]).resize(0);
}

/* Returns the name of the active data buffer */
template<int Rank>
int Wavefunction<Rank>::GetActiveBufferName() const
{
	return ActiveBufferName;
}

/* Set the currently active data buffer and returns the name of the previous
 * active databuffer
 */
template<int Rank>
int Wavefunction<Rank>::SetActiveBuffer(int bufferName)
{
	int oldActiveBufferName = GetActiveBufferName();
	Data.reference(*WavefunctionData[bufferName]);
	ActiveBufferName = bufferName;
	return oldActiveBufferName;
}

/* Get a reference to a buffer */
template<int Rank>
typename Wavefunction<Rank>::DataArray& Wavefunction<Rank>::GetData(int bufferName)
{
	return *WavefunctionData[bufferName];
}


template<int Rank>
const typename Wavefunction<Rank>::DataArray& Wavefunction<Rank>::GetData(int bufferName) const
{
	return *WavefunctionData[bufferName];
}


template<int Rank>
void Wavefunction<Rank>::SetData(Wavefunction<Rank>::DataArray &newData)
{
	WavefunctionData[GetActiveBufferName()]->reference(newData);
	Data.reference(newData);
}


template<int Rank>
Wavefunction<Rank>* 
Wavefunction<Rank>::Copy() const
{
	/* Set up representations and stuff */
	Wavefunction<Rank>* newPsi = new Wavefunction();
	newPsi->SetRepresentation(this->Repr->Copy());
	
	/* Allocate data */
	int bufferName = newPsi->AllocateData(Data.shape());
	newPsi->SetActiveBuffer(bufferName);

	/* Copy data */
	newPsi->Data = this->Data;

	return newPsi;
}

template<int Rank>
Wavefunction<Rank>* 
Wavefunction<Rank>::CopyDeep() const
{
	/* Set up representations and stuff */
	Wavefunction<Rank>* newPsi = new Wavefunction();
	newPsi->SetRepresentation(this->Repr->Copy());
	
	/* Allocate data */
	for (size_t i = 0; i < this->WavefunctionData.size(); i++)
	{
		//Allocate data buffer in new wavefunction
		DataArray oldData ( GetData(i) );
		int bufferName = newPsi->AllocateData(oldData.shape());
		if (bufferName != (int)i)
		{
			throw std::runtime_error("What! something is wrong in Wavefunction::CopyDeep()");
		}

		//Copy data buffer to new wavefunction
		DataArray newData ( newPsi->GetData(bufferName) );
		newData = oldData;
	}

	//Set active buffer on the new wavefunction
	newPsi->SetActiveBuffer(this->GetActiveBufferName());

	return newPsi;
}

template<int Rank>
double Wavefunction<Rank>::GetNorm() const
{
	double localNorm = GetLocalNorm();
	return sqrt(GetRepresentation()->GetDistributedModel()->GetGlobalSum(localNorm));
}

template<int Rank>
double Wavefunction<Rank>::GetLocalNorm() const
{
	return GetRepresentation()->InnerProduct(*this, *this).real();
}

template<int Rank>
cplx Wavefunction<Rank>::InnerProduct(const Wavefunction<Rank> &psi) const
{
	cplx localInnerProd = LocalInnerProduct(psi);
	return GetRepresentation()->GetDistributedModel()->GetGlobalSum(localInnerProd);
}

template<int Rank>
cplx Wavefunction<Rank>::LocalInnerProduct(const Wavefunction<Rank> &psi) const
{
	return GetRepresentation()->InnerProduct(psi, *this);
}



template class Wavefunction<1>;
template class Wavefunction<2>;
template class Wavefunction<3>;
template class Wavefunction<4>;


