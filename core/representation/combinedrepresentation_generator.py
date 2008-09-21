"""

Generates C++ code to calculate the inner product efficiently for different ranks.
my friend temlate programming failed me in this case, and i had to use python to generate 
c++ code

The inner product of two r-dimensional vectors |psi1> and |psi2>, with overlap matrix S_i for
each rank i is in matrix notation

innerProduct = <psi1| S_0 (X) S_1 (X) ... (X) S_r |psi2>

where S_0 (X) S_1 is the tensor product between S_0 and S_1

the overlap matrices are in general assumed to be hermitian and band limited with bandwidth k_i. 
For orthogonal basises (currently all other basises than BSplines) k_i = 1, and the weights
can be placed on the diagonal. It is therefore extemely important to do this case efficiently. 

Using summation notation, we can write the inner product as

innerProduct = sum_{i0, ...,ir} sum_{band0, ..., band1} 
					psi1(i0, ..., ir) * psi2(i0+band0, ..., ir+band_r) 
					* S_0(i0, band0) * ... * S_r(ir, bandr)

Its a sum over all indices of the wavefunction and for each index it is a sum over all
bands. This means it will be a 2*r nested loops. These can be organized in a number of ways. 
Currently we will only consider two cases
1) The loop over all indices is the outer loop, while the loop over all bands is the inner loop
2) The loop over all bands is the outer loop, while the loop over indices is the inner loop
In effect this is to switch the ordering of the two summations above.

1) will be fastest when there are a lot of bands, and the vectors are large, as we will reuse the 
vectors in cache more easily than with 2). However, when there are few bands, the overhead connected
with the innermost loops will slow down the computation.

We will use extensive pointer-foo instead of indexing blitz arrays, as it enables us to exploit the 
forward-only way we are accessing the data

"""

def PrettyPrint(str):
	outStr = ""
	lines = [line.strip() for line in str.split("\n")]
	
	indent = 0
	for line in lines:
		if line == "}":
			indent -= 1

		outStr += "\t" * indent	+ line + "\n"

		if line == "{":
			indent += 1

	return outStr


#---------------------------------------------------------------------------------------------------------
#    Algorithm 1
#---------------------------------------------------------------------------------------------------------

def GenerateAlgorithm1(rank):
	str = """
		// Generated function. Do Not Touch!
		template<> cplx CombinedRepresentation<%(rank)i>::InnerProductImpl_Algo1(DataArray d1, DataArray d2)
		{
			cplx innerProduct = 0;
	""" % { "rank": rank }

	str += RecurseAlgorithm1OuterLoop(rank, 0)

	str += """
			return innerProduct;
		}
	"""
	return str


def RecurseAlgorithm1OuterLoop(rank, curRank):
	indexStr  = ", ".join(["i%i" % (i) for i in range(rank)])
	indexStartBandStr  = ", ".join(["i%i + startBand%i" % (i, i) for i in range(rank)])

	str = ""

	#if this is the outermost loop
	if curRank == 0:
		for i in range(rank):
			str += """
		
				blitz::Array<double, 2> overlapFull%(rank)i = this->GetGlobalOverlapMatrix(%(rank)i)->GetOverlapFullCol();
				int bw%(rank)i = (overlapFull%(rank)i.extent(1) + 1) / 2;
			""" % { "rank" : i }


	str += """
	for (int i%(rank)i=0; i%(rank)i<d1.extent(%(rank)i); i%(rank)i++)
	{
	""" % { "rank": curRank }

	#If this is the innermost loop
	if curRank == rank-1:
		str += """
			cplx curValue = conj(d1(%(index)s));
		""" % { "index": indexStr }
		for i in range(rank):
			str += """
				int startBand%(rank)i = std::max(-bw%(rank)i+1, -i%(rank)i);
				int endBand%(rank)i = std::min(bw%(rank)i, d1.extent(%(rank)i)-i%(rank)i);
				double* __restrict__ overlapStartPtr%(rank)i = & overlapFull%(rank)i(i%(rank)i, startBand%(rank)i+bw%(rank)i-1);

			""" % { "rank": i }

		str += """
			cplx* __restrict__ d2Ptr0 = & d2(%(indexStartBand)s);
		""" % { "indexStartBand": indexStartBandStr }

		str += RecurseAlgorithm1InnerLoop(rank, 0)
	else:
		str += RecurseAlgorithm1OuterLoop(rank, curRank+1)

	str += """
	}
	"""

	return str


def RecurseAlgorithm1InnerLoop(rank, curRank):
	str = """
	double* __restrict__ overlapPtr%(rank)i = overlapStartPtr%(rank)i;
	for (int band%(rank)i=startBand%(rank)i; band%(rank)i<endBand%(rank)i; band%(rank)i++)
	{
	""" % { "rank": curRank }

	if curRank == 0:
		str += """
			double overlap%(rank)i = (*overlapPtr0++);
		""" % { "rank": curRank }
	else:
		str += """
			double overlap%(rank)i = overlap%(prevRank)i * (*overlapPtr%(rank)i++);
		""" % { "rank": curRank, "prevRank": curRank-1 }

	if curRank != rank-1:
		str += """
			cplx* __restrict__ d2Ptr%(nextRank)i = d2Ptr%(rank)i;
			d2Ptr%(rank)i += d2.stride(%(rank)i);
		""" % { "rank": curRank, "nextRank": curRank+1 }

		str += RecurseAlgorithm1InnerLoop(rank, curRank+1)

	else:
		str += """
			innerProduct += curValue * (*d2Ptr%(rank)i++) * overlap%(rank)i;
		""" % { "rank": curRank }

	str += """
	}
	"""

	return str


	
#---------------------------------------------------------------------------------------------------------
#    Algorithm 2
#---------------------------------------------------------------------------------------------------------

def GenerateAlgorithm2(rank):
	str = """
		// Generated function. Do Not Touch!
		template<> cplx CombinedRepresentation<%(rank)i>::InnerProductImpl_Algo2(DataArray d1, DataArray d2)
		{
			cplx innerProduct = 0;
	""" % { "rank": rank }

	str += RecurseAlgorithm2OuterLoop(rank, 0)

	str += """
			return innerProduct;
		}
	"""
	return str


def RecurseAlgorithm2OuterLoop(rank, curRank):
	index1  = ", ".join(["start%i" % (i) for i in range(rank)])
	index2  = ", ".join(["start%i + band%i" % (i, i) for i in range(rank)])

	str = ""

	#if this is the outermost loop
	if curRank == 0:
		for i in range(rank):
			str += """
				blitz::Array<double, 2> overlapFull%(rank)i = this->GetGlobalOverlapMatrix(%(rank)i)->GetOverlapFullRow();
				int bw%(rank)i = (overlapFull%(rank)i.extent(0) + 1) / 2;
			""" % { "rank" : i }


	str += """
	for (int band%(rank)i=-bw%(rank)i+1; band%(rank)i<bw%(rank)i; band%(rank)i++)
	{
		int start%(rank)i=0;
		int end%(rank)i=d1.extent(%(rank)i);
		if (band%(rank)i < 0)
		{
			start%(rank)i = -band%(rank)i;
		}
		else
		{
			end%(rank)i -= band%(rank)i;
		}
		double* __restrict__ overlapStartPtr%(rank)i = & overlapFull%(rank)i(band%(rank)i+bw%(rank)i-1, start%(rank)i);
	""" % { "rank": curRank }

	#If this is the innermost loop
	if curRank == rank-1:
		str += """
			cplx* __restrict__ d1Ptr0 = & d1(%(index1)s);
			cplx* __restrict__ d2Ptr0 = & d2(%(index2)s);
		""" % { "index1": index1, "index2": index2 }

		str += RecurseAlgorithm2InnerLoop(rank, 0)
	else:
		str += RecurseAlgorithm2OuterLoop(rank, curRank+1)

	str += """
	}
	"""

	return str


def RecurseAlgorithm2InnerLoop(rank, curRank):
	str = """
	double* __restrict__ overlapPtr%(rank)i = overlapStartPtr%(rank)i;
	for (int i%(rank)i=start%(rank)i; i%(rank)i<end%(rank)i; i%(rank)i++)
	{
	""" % { "rank": curRank }

	if curRank == 0:
		str += """
			double overlap%(rank)i = (*overlapPtr0++);
		""" % { "rank": curRank }
	else:
		str += """
			double overlap%(rank)i = overlap%(prevRank)i * (*overlapPtr%(rank)i++);
		""" % { "rank": curRank, "prevRank": curRank-1 }

	if curRank != rank-1:
		str += """
			cplx* __restrict__ d1Ptr%(nextRank)i = d1Ptr%(rank)i;
			d1Ptr%(rank)i += d1.stride(%(rank)i);

			cplx* __restrict__ d2Ptr%(nextRank)i = d2Ptr%(rank)i;
			d2Ptr%(rank)i += d2.stride(%(rank)i);
		""" % { "rank": curRank, "nextRank": curRank+1 }

		str += RecurseAlgorithm2InnerLoop(rank, curRank+1)

	else:
		str += """
			innerProduct += conj(*d1Ptr%(rank)i++) * (*d2Ptr%(rank)i++) * overlap%(rank)i;
		""" % { "rank": curRank }

	str += """
	}
	"""

	return str


#---------------------------------------------------------------------------------------------------------
#    Algorithm 
#---------------------------------------------------------------------------------------------------------

def GenerateAlgorithm3(rank):
	def GenerateRankOneBody(storageId):
		methodBody = ""
		for i in range(rank):
			methodBody += """
				if (rank == %(curRank)s)
				{
					TensorPotential::TensorPotentialMultiply_%(storageIdString)s_Wrapper(potential, scaling, source, dest);
				}
			""" % { \
				"rank": rank, \
				"curRank": i, \
				"storageIdString": "_".join(["Ident"]*i + [storageId] + ["Ident"]*(rank-i-1)), \
			}
		return methodBody

	str = ""
	for storageId in ["Band", "BandNH", "Dense"]:
		str += """
		void TensorPotentialMultiply_Rank1_%(storageId)s(int rank, blitz::Array<cplx, %(rank)i> potential, double scaling, blitz::Array<cplx, %(rank)i> &source, blitz::Array<cplx, %(rank)i> &dest)
		{
			%(body)s
		}
		""" % { \
			"rank": rank, \
			"storageId": storageId, \
			"body": GenerateRankOneBody(storageId), \
		}

	#Rank-1 wrapper for Distributed tensors. This by itself because it requires two extra parameters, globalSize and bands
	str += """
	void TensorPotentialMultiply_Rank1_Distr(int rank, blitz::Array<cplx, %(rank)i> potential, double scaling, blitz::Array<cplx, %(rank)i> &source, blitz::Array<cplx, %(rank)i> &dest, int globalSize, int bands)
	{
	""" % { "rank": rank }
	for i in range(rank):
		str += """
			if (rank == %(curRank)s)
			{
				TensorPotential::TensorPotentialMultiply_%(storageIdString)s_Wrapper(potential, scaling, source, dest, globalSize, bands);
			}
		""" % { \
			"rank": rank, \
			"curRank": i, \
			"storageIdString": "_".join(["Ident"]*i + ["Distr"] + ["Ident"]*(rank-i-1)), \
		}
	str += "}\n"


	return str
		

for i in range(1,4+1):
	print '#include "combinedrepresentation.h"'
	print '#include "../tensorpotential/tensorpotentialmultiply_wrapper.h"'
	print PrettyPrint(GenerateAlgorithm1(i))
	print PrettyPrint(GenerateAlgorithm2(i))
	print PrettyPrint(GenerateAlgorithm3(i))

