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
	
	def GetBasisPairCount(self):
		N = self.BSplineObject.NumberOfBSplines
		k = self.BSplineObject.MaxSplineOrder
		return 2 * N * k - (k) * (k-1) - N

		
	def GetBasisPairs(self):
		count = self.GetBasisPairCount()

		N = self.BSplineObject.NumberOfBSplines
		k = self.BSplineObject.MaxSplineOrder
		pairs = zeros((self.GetBasisPairCount(), 2), dtype=int32)
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
	
	def GetBasisPairCount(self):
		N = self.BSplineObject.NumberOfBSplines
		k = self.BSplineObject.MaxSplineOrder
		return N * k 

		
	def GetBasisPairs(self):
		count = self.GetBasisPairCount()

		N = self.BSplineObject.NumberOfBSplines
		k = self.BSplineObject.MaxSplineOrder

		pairs = zeros((self.GetBasisPairCount(), 2), dtype=int32)
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


class BasisfunctionBSpline(BasisfunctionBase):
	"""
	Basisfunction class for BSplines
	see BasisfunctionBase for details
	"""

	def ApplyConfigSection(self, configSection):
		self.SetupBasis( InitBSpline(configSection) )

	def SetupBasis(self, bsplineObject):
		self.BSplineObject = bsplineObject

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
		else:
			raise UnsupportedGeometryException("Geometry '%s' not supported by BasisfunctionBSpline" % geometryName)

	def RepresentPotentialInBasis(self, source, dest, rank, geometryInfo, differentiation):
		pairs = geometryInfo.GetBasisPairs()
		storageId = geometryInfo.GetStorageId()

		if (storageId == "Band" or storageId == "Herm") and differentiation % 1 == 1:
			raise Exception("Cannot use hermitian storage and first order differentiation matrix (antihermitian)")

		if storageId == "Ident":
			sourceSlice = [slice(0, None, None)]*len(source.shape)
			sourceSlice[rank] = slice(0, 1, None)
			dest[:] = source[sourceSlice]
		else:
			core.RepresentPotentialInBasisBSpline(self.BSplineObject, source, dest, pairs, storageId, rank, differentiation)



