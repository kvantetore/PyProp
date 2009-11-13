#------------------------------------------------------------------------------------
#                       BSpline
#------------------------------------------------------------------------------------


class GeometryInfoBSplineBanded(GeometryInfoBase):
	"""
	Geometry information for BSpline geometries
	using a very straight forward index-pair map.

	This is only provided to compare with the GeometryInfoBSplineBandedBlas
	which uses BLAS to obtain superior efficiency.

	"""
	def __init__(self, bsplineObject):
		#Set member variables 
		self.BSplineObject = bsplineObject	

	def UseGridRepresentation(self):
		return True
	
	def GetGlobalBasisPairCount(self):
		N = self.BSplineObject.NumberOfBSplines
		k = self.BSplineObject.MaxSplineOrder
		return 2 * N * k - (k) * (k-1) - N

		
	def GetBasisPairs(self):
		count = self.GetGlobalBasisPairCount()

		N = self.BSplineObject.NumberOfBSplines
		k = self.BSplineObject.MaxSplineOrder
		pairs = zeros((self.GetGlobalBasisPairCount(), 2), dtype=int32)
		pairs[:,0] = 0
		pairs[:,1] = N-1

		index = 0
		for i in xrange(N):
			for j in xrange(max(0, i-k+1), min(i+k, N)):
				pairs[index, 0] = i
				pairs[index, 1] = j
				index+=1

		if index != pairs.shape[0]:
			raise Execption()

		return pairs
	
	def GetStorageId(self):
		return "Simp"

	def GetMultiplyArguments(self, psi):
		if not hasattr(self, "BasisPairs"):
			self.BasisPairs = self.GetBasisPairs()
		return [self.BasisPairs]

class GeometryInfoBSplineBandedBlas(GeometryInfoBase):
	"""
	Geometry information for BSpline geometries, the potential
	is stored in the BLAS hermitian banded format such that 
	multiplication can be done with blas calls
	"""
	def __init__(self, bsplineObject):
		#Set member variables 
		self.BSplineObject = bsplineObject	

	def UseGridRepresentation(self):
		return True
	
	def GetGlobalBasisPairCount(self):
		N = self.BSplineObject.NumberOfBSplines
		k = self.BSplineObject.MaxSplineOrder
		return N * k 

		
	def GetBasisPairs(self):
		count = self.GetGlobalBasisPairCount()

		N = self.BSplineObject.NumberOfBSplines
		k = self.BSplineObject.MaxSplineOrder

		pairs = zeros((self.GetGlobalBasisPairCount(), 2), dtype=int32)
		pairs[:,0] = 0
		pairs[:,1] = N-1

		index = 0
		for i in xrange(N):
			for j in xrange(k):
				#for j in xrange(i, min(i+k, N)):
				if j+i < N:
					pairs[index, 0] = i
					pairs[index, 1] = i+j
				else:
					pairs[index, 0] = 0
					pairs[index, 1] = N-1
				index+=1

		if index != pairs.shape[0]:
			raise Exeption()

		return pairs
	
	def GetStorageId(self):
		return "Band"

	def GetMultiplyArguments(self, psi):
		return []


class GeometryInfoBsplineBandedDistributed(GeometryInfoDistributedBase):
	"""
	Geometry for a banded B-spline matrix with support for parallel
	matvec through Epetra.
	"""
	def __init__(self, repr, useGrid):
		GeometryInfoDistributedBase.__init__(self)
		#Set member variables 
		self.UseGrid = useGrid
		self.BaseRank = repr.GetBaseRank()
		self.Representation = repr
		self.BSplineObject = repr.GetBSplineObject()	

	def UseGridRepresentation(self):
		return self.UseGrid

	def GetGlobalBasisPairCount(self):
		N = self.BSplineObject.NumberOfBSplines
		k = self.BSplineObject.MaxSplineOrder
		return 2 * N * k - (k) * (k-1) - N

	def SetupBasisPairs(self):
		distr = self.Representation.GetDistributedModel()
		rank = self.Representation.GetBaseRank()
		localRange = distr.GetLocalIndexRange(self.GetGlobalBasisPairCount(), rank)

		count = self.GetGlobalBasisPairCount()

		N = self.BSplineObject.NumberOfBSplines
		k = self.BSplineObject.MaxSplineOrder
		pairs = zeros((self.GetGlobalBasisPairCount(), 2), dtype=int32)
		pairs[:,0] = 0
		pairs[:,1] = N-1

		index = 0
		for i in xrange(N):
			for j in xrange(max(0, i-k+1), min(i+k, N)):
				pairs[index, 0] = i
				pairs[index, 1] = j
				index+=1

		if index != pairs.shape[0]:
			raise Execption()

		#store pairs
		self.GlobalIndexPairs = pairs
		self.LocalBasisPairIndices = r_[localRange]
		self.LocalIndexPairs = pairs[localRange]

		return pairs

	def GetStorageId(self):
		return "Distr"

	def GetMultiplyArguments(self, psi):
		raise Exception("GetMultiplyArguments not implemented, use Epetra for matvec!")

	def SetupTempArrays(self, psi):
		raise Exception("SetupTempArrays not implemented, use Epetra for matvec!")


class BasisfunctionBSpline(BasisfunctionBase):
	"""
	Basisfunction class for BSplines
	see BasisfunctionBase for details
	"""

	def ApplyConfigSection(self, configSection):
		self.SetupBasis( InitBSpline(configSection) )

	def SetupBasis(self, repr):
		self.Representation = repr
		#self.BSplineObject = bsplineObject
		self.BSplineObject = repr.GetBSplineObject()

	def GetGridRepresentation(self):
		repr = core.BSplineGridRepresentation()
		repr.SetupRepresentation(self.BSplineObject)
		return repr

	def GetBasisRepresentation(self):
		repr = core.BSplineRepresentation()
		repr.SetupRepresentation(self.BSplineObject)
		return repr

	def GetGeometryInfo(self, geometryName):
		geom = geometryName.lower().strip()

		BasisSize = self.BSplineObject.NumberOfBSplines
		BandCount = self.BSplineObject.MaxSplineOrder-1

		if geom == "identity":
			PrintOut( "WARNING: Identity geometry might not do what you expect for BSplines" )
			return GeometryInfoCommonIdentity(True)
		elif geom == "banded-nonhermitian":
			return GeometryInfoCommonBandedNonHermitian(BasisSize, BandCount, True)
		elif geom == "banded-packed":
			return GeometryInfoBSplineBanded(self.BSplineObject)
		elif geom == "banded-hermitian":
			return GeometryInfoBSplineBandedBlas(self.BSplineObject)
		elif geom == "dense":
			return GeometryInfoCommonDense(BasisSize, True)
		elif geom == "hermitian":
			return GeometryInfoCommonDense(BasisSize, True)
		elif geom == "banded-bspline-distributed":
			return GeometryInfoBsplineBandedDistributed(self.Representation, True)
		else:
			raise UnsupportedGeometryException("Geometry '%s' not supported by BasisfunctionBSpline" % geometryName)

	def RepresentPotentialInBasis(self, source, dest, rank, geometryInfo, differentiation):
		pairs = geometryInfo.GetGlobalBasisPairs()
		storageId = geometryInfo.GetStorageId()

		if (storageId == "Band" or storageId == "Herm") and differentiation % 1 == 1:
			raise Exception("Cannot use hermitian storage and first order differentiation matrix (antihermitian)")

		if storageId == "Ident":
			sourceSlice = [slice(0, None, None)]*len(source.shape)
			sourceSlice[rank] = slice(0, 1, None)
			dest[:] = source[sourceSlice]
		else:
			core.RepresentPotentialInBasisBSpline(self.BSplineObject, source, dest, pairs, storageId, rank, differentiation)



