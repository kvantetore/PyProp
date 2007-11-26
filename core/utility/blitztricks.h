#ifndef BLITZTRICKS_H
#define BLITZTRICKS_H

/* 
 * Maps a blitz array down to a new array of rank 2. The
 * two resulting arrays will have the same size, and share
 * the same buffer, but different shapes.
 * The parameter firstRankCount specifies the number of ranks in
 * data which will be mapped to the first rank in the output array
 * 
 * i.e. given the array
 * Array<T, 4> A(4, 4, 4, 4)
 * 
 * MapToRank2(A, 1) will return a 4 x 64 array
 * MapToRank2(A, 2) will return a 16 x 16 array
 * MapToRank2(A, 3) will return a 64 x 4 array
 *
 * WARNING: Under the current implementation, the lifetime of the returned array
 * is determined by the lifetime of the input array. In effect: be sure that the 
 * returned array is only used when data is stil allocated
 */
template<int Rank, class T>
blitz::Array<T, 2> MapToRank2(blitz::Array<T, Rank> &data, int firstRankCount)
{
	using namespace blitz;

    //Check that the array is in C-order
	for (int i=0; i<Rank; i++)
	{
		if (data.ordering(i) != (Rank - 1 - i))
		{
			throw std::runtime_error("data is not in C-style ordering");
		}
	}
	
	//the size of the first rank is first firstRankCount rows combined
	int firstRankSize = 1;
	for (int i=0; i<firstRankCount; i++)
	{
		firstRankSize *= data.extent(i);
	}
	//The size of the second rank is the rest of the array
	int secondRankSize = data.size() / firstRankSize;

	//Set up shape and stride vectors
	TinyVector<int, 2> shape(firstRankSize, secondRankSize);
	TinyVector<int, 2> stride(secondRankSize, 1);

	//Create an array
	//TODO: fix this so it doesn't reference the underlying vector, but rather
	//the reference counted memory pool. That way we will have deallocation 
	//as expected
	Array<cplx, 2> ret(data.data(), shape, stride, neverDeleteData);

	return ret;
}



template<class T, int Rank>
blitz::Array<T, 3> MapToRank3(blitz::Array<T, Rank> &array, int firstRankCount, int secondRankCount)
{
	using namespace blitz;

    //Check that the array is in C-order
	for (int i=0; i<Rank; i++)
	{
		if (array.ordering(i) != (Rank - 1 - i))
		{
			throw std::runtime_error("data is not in C-style ordering");
		}
	}
	
	int firstRankSize = 1;
	for (int i=0; i<firstRankCount;i++)
	{
		firstRankSize *= array.extent(i);
	}

	int secondRankSize = 1;
	for (int i=0; i<secondRankCount; i++)
	{
		int curRank = i + firstRankCount;
		secondRankSize *= array.extent(curRank);
	}

	int thirdRankSize = array.size() / (secondRankSize * firstRankSize);
	
	//Set up shape and stride vectors
	TinyVector<int, 3> shape(firstRankSize, secondRankSize, thirdRankSize);
	TinyVector<int, 3> stride(secondRankSize*thirdRankSize, thirdRankSize, 1);

	//Create an array
	//TODO: fix this so it doesn't reference the underlying vector, but rather
	//the reference counted memory pool. That way we will have deallocation 
	//as expected
	//The current implementation requires that the lifetime of the input array
	//exceeds the lifetime of the output array
	Array<cplx, 3> ret(array.data(), shape, stride, neverDeleteData);

	return ret;
}


template <class T, int Rank>
blitz::Array<T, 1> MapToRank1(blitz::Array<T, Rank> array)
{
	using namespace blitz;

	TinyVector<int, 1> shape( array.size() );
	TinyVector<int, 1> stride( 1 );

	Array<T, 1> ret(array.data(), shape, stride, neverDeleteData);

	return ret;
}

template <class T, int Rank>
void OuterProduct( const blitz::TinyVector< blitz::Array<T, 1>, Rank> &in, blitz::Array<T, Rank> &out)
{
	cout << "Error: OuterProduct is not implemented for rank = " << Rank << endl;
	throw std::runtime_error("Outer Product not implemented for this rank");
}


/*
 * OuterProduct is the N-dimensional extension of the numpy method outer.
 * It takes N 1-dimensional arrays in(0)...in(n) and forms the outer product
 *
 * out(x0, ..., xn) = in(0)(x0) * ... * in(n)(xn)
 *
 * the shape of in and out should be such that
 * extent(i) = in(i).size()
 */

template <class T>
void OuterProduct( const blitz::TinyVector< blitz::Array<T, 1>, 1> &in, blitz::Array<T, 1> &out)
{
	out = in(0)(blitz::tensor::i);
}

template <class T>
void OuterProduct( const blitz::TinyVector< blitz::Array<T, 1>, 2> &in, blitz::Array<T, 2> &out)
{
	out = in(0)(blitz::tensor::i) 
		* in(1)(blitz::tensor::j);
}

template <class T>
void OuterProduct( const blitz::TinyVector< blitz::Array<T, 1>, 3> &in, blitz::Array<T, 3> &out)
{
	out = in(0)(blitz::tensor::i) 
		* in(1)(blitz::tensor::j)
		* in(2)(blitz::tensor::k);
}

template <class T>
void OuterProduct( const blitz::TinyVector< blitz::Array<T, 1>, 4> &in, blitz::Array<T, 4> &out)
{
	out = in(0)(blitz::tensor::i) 
		* in(1)(blitz::tensor::j)
		* in(2)(blitz::tensor::k)
		* in(3)(blitz::tensor::k);
}

#endif

