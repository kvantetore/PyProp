from numpy import int32
from numpy import array
from numpy import asarray
from numpy import zeros
from numpy import ones
import re

execfile("SimpleDistributed.py")

def PrintOut(args):
	if pyprop.ProcId == 0:
		print args

#------------------------------------------------------------------------------------
#                       Interfaces
#------------------------------------------------------------------------------------

class UnsupportedGeometryException(Exception):
	def __init__(self, string):
		Exception.__init__(self, string)

class GeometryInfoBase(object):
	"""
	base for all Geometry information objects
	a GeometryInfo object gives information about how the geometry behaves
	"""

	def UseGridRepresentation(self):
		"""
		Whether the grid representation or the basis representation should
		be used to evaluate the potential
		"""
		raise NotImplementedException()

	def GetBasisPairCount(self):
		"""
		Returns the number of stored matrix elements
		"""
		raise NotImplementedException()

	def GetBasisPairs(self):
		"""
		Returns the row-col index pairs for the matrix elements as a 
		two dimensional array. Each row is an index pair, with 
		the first column being the row index, and the second being
		the col-index
		"""
		raise NotImplementedException()

	def GetStorageId(self):
		"""
		Returns the string defining how the index pairs are ordered
		"""
		raise NotImplementedException()	

	def GetMultiplyArguments(self):
		"""
		Returns the arguments needed for TensorMultiply
		"""
		raise NotImplementedException()

	def HasParallelMultiply(self):
		"""
		Whether this storage supports TensorMultiply when the rank 
		is distributed
		"""
		return False


class BasisfunctionBase(object):
	"""
	Abstract class which bases all basis function classes.
	Its interface provides 
	"""

	def ApplyConfigSection(self, config):
		raise NotImplementedException()

	def GetGridRepresentation(self):
		"""
		Creates a one dimensional representation for this rank
		in the grid representation
		"""
		raise NotImplementedException()

	def GetBasisRepresentation(self):
		"""
		Creates a one dimensional representation for this rank
		in the basis representation
		"""
		raise NotImplementedException()

	def GetGeometryInfo(self, geometryName):
		raise NotImplementedException()

	def RepresentPotentialInBasis(self, source, dest, rank, geometryInfo, differentiation):
		raise NotImplementedException


#------------------------------------------------------------------------------------
#                       Common geometries
#------------------------------------------------------------------------------------

class GeometryInfoCommonDense(GeometryInfoBase):
	"""
	General geometry information for dense matrices
	"""
	def __init__(self, rankCount, useGrid):
		#Set member variables 
		self.RankCount = rankCount
		self.UseGrid = useGrid

	def UseGridRepresentation(self):
		return self.UseGrid
	
	def GetBasisPairCount(self):
		return self.RankCount**2 
		
	def GetBasisPairs(self):
		pairs = zeros((self.GetBasisPairCount(), 2), dtype=int32)
		index = 0
		for i in xrange(self.RankCount):
			for j in xrange(self.RankCount):
				pairs[index, 0] = i
				pairs[index, 1] = j
				index+=1

		return pairs

	def GetStorageId(self):
		return "Simp"

	def GetMultiplyArguments(self):
		if not hasattr(self, "BasisPairs"):
			self.BasisPairs = self.GetBasisPairs()
		return [self.BasisPairs]

class GeometryInfoCommonDenseHermitian(GeometryInfoBase):
	"""
	General geometry information for dense matrices
	"""
	def __init__(self, rankCount, useGrid):
		#Set member variables 
		self.RankCount = rankCount	
		self.UseGrid = useGrid

	def UseGridRepresentation(self):
		return self.UseGrid
	
	def GetBasisPairCount(self):
		N = self.RankCount
		return N * (N + 1) / 2
		
	def GetBasisPairs(self):
		pairs = zeros((self.GetBasisPairCount(), 2), dtype=int32)
		index = 0
		for i in xrange(self.RankCount):
			for j in xrange(i, self.RankCount):
				pairs[index, 0] = i
				pairs[index, 1] = j
				index+=1

		return pairs

	def GetStorageId(self):
		return "Herm"

	def GetMultiplyArguments(self):
		if not hasattr(self, "BasisPairs"):
			self.BasisPairs = self.GetBasisPairs()
		return [self.BasisPairs]

class GeometryInfoCommonDiagonal(GeometryInfoBase):
	"""
	General geometry information for diagonal matrices.

	The matrix is diagonal in the sense that all index pairs
	except where row==col, gives a zero contribution. In this case,
	only the row==col elements are stored

	Matrices are usually diagonal when the basis functions are the
	eigenvectors of the potential. They will then also set UseGridRepresentation
	to false. This is the case for the angular momentum
	which is diagonal in the spherical harmonic representation

	V(l, r) =  l*(l+1) / (2 * m * r*r)

	"""
	def __init__(self, rankCount, useGrid):
		#Set member variables 
		self.RankCount = rankCount	
		self.UseGrid = useGrid

	def UseGridRepresentation(self):
		return self.UseGrid
	
	def GetBasisPairCount(self):
		return self.RankCount 
		
	def GetBasisPairs(self):
		pairs = zeros((self.GetBasisPairCount(), 2), dtype=int32)
		index = 0
		for i in xrange(self.RankCount):
			pairs[index, 0] = i
			pairs[index, 1] = i
			index+=1

		return pairs

	def GetStorageId(self):
		return "Diag"

	def GetMultiplyArguments(self):
		return []

	def HasParallelMultiply(self):
		return True


class GeometryInfoCommonBandedDistributed(GeometryInfoBase):
	"""
	Geometry for a banded matrix with support for parallel
	TensorMultiply
	"""
	def __init__(self, rankCount, bandCount, useGrid):
		#Set member variables 
		self.RankCount = rankCount	
		self.BandCount = bandCount
		self.UseGrid = useGrid

	def UseGridRepresentation(self):
		return self.UseGrid
	
	def GetBasisPairCount(self):
		return self.RankCount * (2 * self.BandCount + 1)
		
	def GetBasisPairs(self):
		pairs = zeros((self.GetBasisPairCount(), 2), dtype=int32)
		index = 0
		for i in xrange(self.RankCount):
			for j in xrange(i-self.BandCount, i+self.BandCount+1):
				if 0 <= j < self.RankCount:
					packedRow = j
					packedCol = self.BandCount - j + i

					index = packedRow * (2 * self.BandCount + 1) + packedCol
					pairs[index, 0] = i
					pairs[index, 1] = j

		return pairs

	def GetStorageId(self):
		return "Distr"

	def GetMultiplyArguments(self):
		return [self.RankCount, self.BandCount]

	def HasParallelMultiply(self):
		return True

class GeometryInfoCommonIdentity(GeometryInfoBase):
	"""
	General geometry information for Identity matrices. 

	By Identity in this setting, we mean that it is independent
	of the variable in this rank, i.e it will be Identity in 
	the y rank in the following example

	V(x,y) = f(x)

	For Identity to work properly, the basis should be orthogonal. 
	In this case, the above function will give 
				{  f(x) if i'==i
	V(x,i'i) =  {
				{  0 if i'!=i
	if the basis is not orthogonal, we will get the overlap matrix
	S_{i',i} if i!=i, which will ruin this scheme

	To make less special cases, this class will return 1 
	index pair, (0, 0)
	"""
	def __init__(self, useGrid):
		#Set member variables 
		self.UseGrid = useGrid

	def UseGridRepresentation(self):
		return self.UseGrid
	
	def GetBasisPairCount(self):
		return 1 
		
	def GetBasisPairs(self):
		pairs = zeros((1, 2), dtype=int32)
		pairs[0, 0] = 0
		pairs[0, 1] = 0

		return pairs

	def GetStorageId(self):
		return "Ident"

	def GetMultiplyArguments(self):
		return []


class GeometryInfoCommonBandedNonHermitian(GeometryInfoBase):
	"""
	Geometry information for general banded matrices, the potential
	is stored in the BLAS general banded format such that 
	multiplication can be done with blas calls
	"""
	def __init__(self, rankCount, bandCount, useGrid):
		"""
		rankCount is the number of rows in the matrix
		bandCount is the number of super diagonals.

		The matrix is assumed to have equally many super 
		and sub diagonals, making the total number of diagonals
		2 * bandCount + 1
		"""
		#Set member variables 
		self.RankCount = rankCount
		self.BandCount = bandCount
		self.UseGrid = useGrid

	def UseGridRepresentation(self):
		return self.UseGrid
	
	def GetBasisPairCount(self):
		return self.RankCount * (self.BandCount * 2 + 1) 

	def GetBasisPairs(self):
		count = self.GetBasisPairCount()

		N = self.RankCount
		k = self.BandCount

		pairs = zeros((self.GetBasisPairCount(), 2), dtype=int32)
		pairs[:,0] = 0
		pairs[:,1] = N-1

		index = 0
		for i in xrange(N):
			for j in xrange(i-k, i+k+1):
				if 0 <= j < N:
					blasJ = j
					blasI = k + i - j

					index = blasJ * (2*k+1) + blasI
					pairs[index, 0] = i
					pairs[index, 1] = j

		return pairs
	
	def GetStorageId(self):
		return "BandNH"

	def GetMultiplyArguments(self):
		return []



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

	def GetMultiplyArguments(self):
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

	def GetMultiplyArguments(self):
		return []


from pyprop import InitBSpline

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
		repr = pyprop.core.BSplineGridRepresentation()
		repr.SetupRepresentation(self.BSplineObject)
		return repr

	def GetBasisRepresentation(self):
		repr = pyprop.core.BSplineRepresentation()
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
		elif geom == "banded-old":
			return GeometryInfoBSplineBanded(self.BSplineObject)
		elif geom == "banded":
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
			RepresentPotentialInBasisBSpline(self.BSplineObject, source, dest, pairs, storageId, rank, differentiation)



#------------------------------------------------------------------------------------
#                       Reduced Spherical Harmonic
#------------------------------------------------------------------------------------

class GeometryInfoReducedSphHarmSelectionRule(GeometryInfoBase):
	"""
	Geometry information for Reduced spherical harmonic geometries
	with delta l = +- 1 symetry, that is for potentials on the form
	P(x, theta) = f(x) * cos(theta)
	"""
	def __init__(self, sphericalHarmonicObject, useGrid):
		#Set member variables 
		self.SphericalHarmonicObject = sphericalHarmonicObject	
		self.UseGrid = useGrid

	def UseGridRepresentation(self):
		return self.UseGrid
	
	def GetBasisPairCount(self):
		return self.SphericalHarmonicObject.GetLMax()*2 
		
	def GetBasisPairs(self):
		count = self.SphericalHarmonicObject.GetLMax()+1
		
		pairs = zeros((self.GetBasisPairCount(), 2), dtype=int32)
		index = 0
		for i in xrange(count):
			if i > 0:
				pairs[index, 0] = i
				pairs[index, 1] = i-1
				index+=1
			if i < count-1:
				pairs[index, 0] = i
				pairs[index, 1] = i+1
				index+=1

		return pairs

	def GetStorageId(self):
		return "Simp"

	def GetMultiplyArguments(self):
		if not hasattr(self, "BasisPairs"):
			self.BasisPairs = self.GetBasisPairs()
		return [self.BasisPairs]



class GeometryInfoReducedSphHarmSelectionRuleHermitian(GeometryInfoBase):
	"""
	Geometry information for Reduced spherical harmonic geometries
	with delta l = +- 1 symetry, that is for potentials on the form
	P(x, theta) = f(x) * cos(theta)

	This is the same as SphHarmSelectionRule, only in addition
	exploiting the hermiticity of the potential.
	"""
	def __init__(self, sphericalHarmonicObject, hermitian):
		#Set member variables 
		self.SphericalHarmonicObject = sphericalHarmonicObject	

	def UseGridRepresentation(self):
		return True
	
	def GetBasisPairCount(self):
		return self.SphericalHarmonicObject.GetLMax()
		
	def GetBasisPairs(self):
		basisSize = self.SphericalHarmonicObject.GetLMax()+1
		
		pairs = zeros((self.GetBasisPairCount(), 2), dtype=int32)
		index = 0
		for i in xrange(basisSize-1):
			pairs[index, 0] = i
			pairs[index, 1] = i+1
			index+=1

		return pairs

	def GetStorageId(self):
		return "Herm"

	def GetMultiplyArguments(self):
		if not hasattr(self, "BasisPairs"):
			self.BasisPairs = self.GetBasisPairs()
		return [self.BasisPairs]



class BasisfunctionReducedSphericalHarmonic(BasisfunctionBase):
	"""
	Basisfunction class for BSplines
	see BasisfunctionBase for details
	"""

	def ApplyConfigSection(self, configSection):
		self.SetupBasis(configSection.lmax)

	def SetupBasis(self, lmax):
		self.LMax = lmax
		self.SphericalHarmonicObject = pyprop.core.ReducedSphericalTools()
		self.SphericalHarmonicObject.Initialize(self.LMax)

	def GetGridRepresentation(self):
		repr = pyprop.core.ThetaRepresentation()
		repr.SetupRepresentation(self.LMax)
		return repr

	def GetBasisRepresentation(self):
		repr = pyprop.core.ReducedSphericalHarmonicRepresentation()
		repr.SetupRepresentation(self.LMax)
		return repr

	def GetGeometryInfo(self, geometryName):
		geom = geometryName.lower().strip()

		#Use grid unless the geometry is prefixed with "custom-"
		useGrid = True
		customStr = "custom-"
		if geom.startswith(customStr):
			useGrid = False
			geom = geom[len(customStr):]

		#Select geometry info
		if geom == "identity":
			return GeometryInfoCommonIdentity(useGrid)
		elif geom == "diagonal":
			return GeometryInfoCommonDiagonal(self.LMax+1, False) #diagonal potentials are always in basis repr
		elif geom == "dense":
			return GeometryInfoCommonDense(self.LMax+1, useGrid)
		elif geom == "dipoleselectionrule":
			return GeometryInfoReducedSphHarmSelectionRule(self.SphericalHarmonicObject, useGrid)
		elif geom == "bandeddistributed":
			return GeometryInfoCommonBandedDistributed(self.LMax+1, 1, useGrid)
		else:
			raise UnsupportedGeometryException("Geometry '%s' not supported by BasisfunctionReducedSpherical" % geometryName)

	def RepresentPotentialInBasis(self, source, dest, rank, geometryInfo, differentiation):
		pairs = geometryInfo.GetBasisPairs()
		storageId = geometryInfo.GetStorageId()
		if (storageId == "Band" or storageId == "Herm") and differentiation % 1 == 1:
			raise Exception("Cannot use hermitian storage and first order differentiation matrix (antihermitian)")
		RepresentPotentialInBasisReducedSphericalHarmonic(self.SphericalHarmonicObject, source, dest, pairs, storageId, rank, differentiation)


#------------------------------------------------------------------------------------
#                       Coupled Spherical Harmonic
#------------------------------------------------------------------------------------

class GeometryInfoCoupledSphericalHarmonic(GeometryInfoBase):
	"""
	Geometry information for coupled spherical harmonic geometries
	The selection CoupledSphericalHarmonicRepresentation and a 
	CoupledSphericalSelectionRule is specified.

	"""
	def __init__(self, representation, selectionRule):
		self.SelectionRule = selectionRule
		self.Representation = representation
		self.IndexPairs = None

	def UseGridRepresentation(self):
		return False
	
	def GetBasisPairCount(self):
		return self.GetBasisPairs().shape[0] 
		
	def GetBasisPairs(self):
		if self.IndexPairs == None:
			self.IndexPairs = self.SelectionRule.GetBasisPairs(self.Representation)
		return self.IndexPairs

	def GetStorageId(self):
		return "Simp"

	def GetMultiplyArguments(self):
		return [self.GetBasisPairs()]


class GeometryInfoCoupledSphericalHarmonicDistributed(GeometryInfoBase):
	"""
	Geometry information for coupled spherical harmonic geometries using distributed matvec.
	The CoupledSphericalHarmonicRepresentation and a 
	CoupledSphericalSelectionRule is specified.

	"""
	def __init__(self, representation, selectionRule):
		self.SelectionRule = selectionRule
		self.Representation = representation
		self.LocalIndexPairs = None
		self.MultiplyArguments = None

	def UseGridRepresentation(self):
		return False
	
	def GetBasisPairCount(self):
		return self.GetBasisPairs().shape[0] 
		
	def GetBasisPairs(self):
		if self.LocalIndexPairs == None:
			self.SetupLocalBasisPairs()
		return self.LocalIndexPairs

	def GetStorageId(self):
		return "SimpD"

	def GetMultiplyArguments(self):
		if self.MultiplyArguments == None:
			self.SetupLocalBasisPairs()
		return self.MultiplyArguments
	
	def HasParallelMultiply(self):
		return True

	def SetupLocalBasisPairs(self):
		indexPairs = self.SelectionRule.GetBasisPairs(self.Representation)
		distrib = self.Representation.GetDistributedModel()
		rank = self.Representation.GetBaseRank()
		globalSize = self.Representation.GetFullShape()[0]
	
		distribIndexList = SetupDistributedIndexList(globalSize, indexPairs, distrib, rank)
		stepList = SetupStepList(globalSize, indexPairs, distribIndexList, distrib, rank)

		self.MultiplyArguments = [int(globalSize)] + list(StepListToArray(stepList))
		self.LocalIndexPairs = array([[step.GlobalRow, step.GlobalCol] for step in stepList if step.LocalMatrixIndex!=-1], dtype=int32)




class BasisfunctionCoupledSphericalHarmonic(BasisfunctionBase):
	"""
	Basisfunction class for coupled spherical harmonics Y{L,M,l1,l1} 

	The CoupledSphericalHarmonic representation is never represented 
	on the grid, but rather the matrix elements should be calculated
	analytically, and passed to the TensorGenerator using a custom
	potential evaluator.

	A custom potential evaluator is a class implementing the following 
	methods: 
		UpdatePotential(potentialData, psi, dt, t)
		SetBasisPairs(rank, basisPairs)

	"""

	def ApplyConfigSection(self, configSection):
		self.SetupBasis(configSection)
		
	def SetupBasis(self, basisRepresentation):
		self.BasisRepresentation = basisRepresentation
		self.BasisSize = self.BasisRepresentation.Range.Count()

	def GetGridRepresentation(self):
		raise Exception("Coupled Spherical Harmonics does not support grid representations")

	def GetBasisRepresentation(self):
		return self.BasisRepresentation
		
	def GetGeometryInfo(self, geometryName):
		geom = geometryName.lower().strip()
		if geom == "identity":
			return GeometryInfoCommonIdentity(False)
		elif geom == "diagonal":
			#return GeometryInfoCommonDiagonal(self.BasisSize, False)
			selectionRule = CoupledSphericalSelectionRuleDiagonal()
			return GeometryInfoCoupledSphericalHarmonicDistributed(self.BasisRepresentation, selectionRule)
		elif geom == "dense":
			return GeometryInfoCommonDense(self.BasisSize, False)
		elif geom.startswith("selectionrule_"):
			selectionRuleName = geom[len("selectionrule_"):]
			if selectionRuleName == "r12":
				selectionRule = CoupledSphericalSelectionRuleR12()
			elif selectionRuleName == "linearpolarizedfield":
				selectionRule = CoupledSphericalSelectionRuleLinearPolarizedField()
			else:	
				raise Exception("Unkonwn selection rule %s" % selectionRuleName)
				
			return GeometryInfoCoupledSphericalHarmonicDistributed(self.BasisRepresentation, selectionRule)
		else:
			raise UnsupportedGeometryException("Geometry '%s' not supported by BasisfunctionReducedSpherical" % geometryName)

	def RepresentPotentialInBasis(self, source, dest, rank, geometryInfo, differentiation):
		raise Exception("Coupled Spherical Harmonics are already in a basis and should not be integrated")



#------------------------------------------------------------------------------------
#                       Finite Difference Basis Function
#------------------------------------------------------------------------------------


class BasisfunctionFiniteDifference(BasisfunctionBase):
	"""
	Basisfunction class for finite differences. 

	Finite differences is not really a basis expansion, as it is
	always represented on a grid, but it's practical to have it 
	in this framework, as it makes it easier to switch between 
	FD and B-splines
	"""

	def ApplyConfigSection(self, configSection):
		self.DifferenceOrder = configSection.difference_order
		self.BandCount = (self.DifferenceOrder - 1) / 2

	def SetupBasis(self, repr):
		self.Representation = repr
		baseRank = repr.GetBaseRank()
		self.GridSize = len(self.Representation.GetGlobalGrid(baseRank))
		self.DifferenceOrder = 1
		self.BandCount = (self.DifferenceOrder - 1) / 2
		
	def GetGridRepresentation(self):
		return self.Representation

	def GetBasisRepresentation(self):
		return self.Representation
		
	def GetGeometryInfo(self, geometryName):
		geom = geometryName.lower().strip()
		if geom == "identity":
			return GeometryInfoCommonIdentity(True)
		elif geom == "diagonal":
			return GeometryInfoCommonDiagonal(self.GridSize, True)
		else:
			diffOrderSearch =  re.search("-\d+", geom)
			if not diffOrderSearch:
				raise  UnsupportedGeometryException("BasisfunctionFiniteDifference requires specification of difference order, you said %s" % geometryName)

			self.DifferenceOrder =  eval(diffOrderSearch.group()[1:])
			self.BandCount = (self.DifferenceOrder - 1) / 2
			if geom.startswith == "dense":
				return GeometryInfoCommonDense(self.GridSize, True)
			elif re.search("banded-nonhermitian", geom) or re.search("banded", geom):
				return GeometryInfoCommonBandedNonHermitian(self.GridSize, self.BandCount, True)
			elif re.search("bandeddistributed", geom):
				return GeometryInfoCommonBandedDistributed(self.GridSize, self.BandCount, True)
			else:
				raise UnsupportedGeometryException("Geometry '%s' not supported by BasisfunctionFiniteDifference" % geometryName)

	def RepresentPotentialInBasis(self, source, dest, rank, geometryInfo, differentiation):
		if differentiation == 0:
			diffMatrix = ones((self.GridSize, 1), dtype=complex)
		elif differentiation == 2:
			fd = pyprop.core.FiniteDifferenceHelper()
			grid = self.Representation.GetGlobalGrid(self.Representation.GetBaseRank())
			fd.Setup(grid, self.DifferenceOrder)
			diffMatrix = fd.SetupLaplacianBlasBanded().copy()
		else:
			raise Exception("Finite Difference currently only supports diff of order 2")

		indexPairs = geometryInfo.GetBasisPairs()
		RepresentPotentialInBasisFiniteDifference(diffMatrix, source, dest, indexPairs, rank)



#------------------------------------------------------------------------------------
#                       TensorPotentialGenerator main
#------------------------------------------------------------------------------------

from pyprop import CreateInstanceRank
from pyprop import CreateDistribution
from pyprop import CreateWavefunctionInstance
from pyprop import CreatePotentialInstance

def CreateBasisFromRepresentation(representation):
	if representation.__class__ == pyprop.core.BSplineRepresentation:
		basis = BasisfunctionBSpline()
		basis.SetupBasis(representation.GetBSplineObject())
	
	elif representation.__class__ == pyprop.core.ReducedSphericalHarmonicRepresentation:
		basis = BasisfunctionReducedSphericalHarmonic()
		basis.SetupBasis(representation.Range.MaxL)
	
	elif representation.__class__ == pyprop.core.CoupledSphericalHarmonicRepresentation:
		basis = BasisfunctionCoupledSphericalHarmonic()
		basis.SetupBasis(representation)

	elif representation.__class__ == pyprop.core.CustomGridRepresentation:
		basis = BasisfunctionFiniteDifference()
		basis.SetupBasis(representation)

	else:
		raise NotImplementedException("Unknown representation %s" % representation)

	return basis


class TensorPotentialGenerator(object):

	def __init__(self, **args):
		self.BasisList = []

		if "config" in args:
			config = args["config"]
			self.Rank = config.Representation.rank
			self.SetupBasisFromConfig(config)

		elif "representation" in args:
			representation = args["representation"]	
			self.Rank = len(representation.GetFullShape())
			self.SetupBasisFromRepresentation(representation)

		else:
			raise Exception("Specify either 'config' or 'representation'")

	def SetupBasisFromConfig(self, config):
		self.BasisList = []
		for i in range(self.Rank):
			#Get the config section for this rank
			rankSection = config.GetSection("Rank%i" % i)
			#Create basisfunction object
			basisObject = rankSection.basis()
			rankSection.Apply(basisObject)
			self.BasisList.append(basisObject)

	def SetupBasisFromRepresentation(self, representation):
		self.BasisList = []
		for i in range(self.Rank):
			subRepr = representation.GetRepresentation(i)
			basisObject = CreateBasisFromRepresentation(subRepr)
			self.BasisList.append(basisObject)

	def GeneratePotential(self, configSection):
		"""
		0) Check the validity of the geometries of the potential
		1) Setup CombinedRepresentation
		2) Create wavefunction
		3) Create potential evaluator
		4) Evaluate potential on grid
		5) Represent potential in basis
		"""

		#Check geometry for each basis
		geometryList = self.GetGeometryList(configSection)

		#1) Create Representation:
		repr = self.SetupRepresentation(geometryList)

		#2) Create a wavefunction in the potential representation. 
		#This will be used to hold the and representation in order 
		#to use the existing potential evaluator framework in pyprop
		psi = CreateWavefunctionInstance(repr)

		#3) Create potential evaluator
		potentialEvaluator = self.SetupPotentialEvaluator(configSection, psi)
		
		#4) Initialize potential with basis pairs for custom potentials
		isCustomPotential = hasattr(potentialEvaluator, "SetBasisPairs")
		potentialShape = list(psi.GetRepresentation().GetInitialShape())
		if isCustomPotential:
			for rank, geometryInfo in enumerate(geometryList):
				if not geometryInfo.UseGridRepresentation():
					#We currently can not be distributed in a basis rank
					if repr.GetDistributedModel().IsDistributedRank(rank) and not geometryInfo.HasParallelMultiply():
						raise Exception("Representation %i can not be distributed during potential evaluation")

					#Get basis pairs from geometry info, and pass it to the potential
					potentialShape[rank] = geometryInfo.GetBasisPairCount()
					potentialEvaluator.SetBasisPairs(rank, geometryInfo.GetBasisPairs())

		#4) Evaluate the potential on the grid
		potentialData = zeros(potentialShape, dtype=complex)
		potentialEvaluator.UpdatePotentialData(potentialData, psi, 0, 0)
		
		# Save a copy of the potential if we're debugging
		debugPotential = getattr(configSection, "debug_potential", False)
		if debugPotential:	
			self.OriginalPotential = potentialData.copy()
			
		#for parallelization
		#We only support the first rank initially distributed
		assert(repr.GetDistributedModel().GetDistribution()[0] == 0)
		assert(len(repr.GetDistributedModel().GetDistribution()) == 1)
		fullShape = repr.GetFullShape()
		if repr.GetDistributedModel().IsSingleProc():
			distribution = []
			origDistribution = []
		else:
			origDistribution = list(repr.GetDistributedModel().GetDistribution())
			distribution = list(repr.GetDistributedModel().GetDistribution())
		transpose = repr.GetDistributedModel().GetTranspose() 

		for distribRank in distribution:
			geomInfo = geometryList[distribRank]
			if not geometryList[distribRank].HasParallelMultiply():
				raise Exception("Distributed Rank %i (%s), has no support for parallel multiplication" % (distribRank, geomInfo.GetStorageId()))

		#5) Represent the potential in the bases
		source = potentialData
		for rank in reversed(range(self.Rank)):
			basis = self.BasisList[rank]
			geometryInfo = geometryList[rank]

			#If we should not use the grid representation for this rank, it means that the potential
			#is already expressed in the basis representation, and thus we dont need to do anything
			#for this rank
			if geometryInfo.UseGridRepresentation():
				differentiation = 0
				if hasattr(configSection, "differentiation%i" % rank):
					differentiation = configSection.Get("differentiation%i" % rank)

				#Check if this rank is distributed. If it is, we must redistribute
				if rank in distribution:
					PrintOut( "Rank %i is distributed" % rank )
					assert(rank != self.Rank-1)
					assert(not rank+1 in distribution)

					#Create new distribution
					distribIndex = distribution.index(rank)
					newDistribution = list(distribution)
					newDistribution[distribIndex] = rank+1

					#Create shape of transposed function and allocate dest buffer
					transposedShape = transpose.CreateDistributedShape(fullShape, array(newDistribution, dtype=int32))
					dest = zeros(transposedShape, dtype=complex)

					#Transpose
					transpose.Transpose(fullShape, source, array(distribution, int), dest, array(newDistribution, int))

					#Use the new buffer
					source = dest
					distribution = newDistribution

				#Calculate dest shape
				destShape = array(source.shape)
				destShape[rank] = geometryInfo.GetBasisPairCount()

				#Update full shape
				fullShape[rank] = geometryInfo.GetBasisPairCount()
			
				#Allocate the destination array
				dest = zeros(destShape, dtype=complex)

				#Represent this rank in the basis
				basis.RepresentPotentialInBasis(source, dest, rank, geometryInfo, differentiation) 

				#Use the destination from this rank as the source to the next
				source = dest

		#done!

		#Transpose back to the original distribution
		if distribution != origDistribution:
			#Create shape of transposed function and allocate dest buffer
			transposedShape = transpose.CreateDistributedShape(fullShape, array(origDistribution, dtype=int32))
			dest = zeros(transposedShape, dtype=complex)

			#Transpose
			transpose.Transpose(fullShape, source, array(distribution, int), dest, array(origDistribution, int))

			#Use the new buffer
			source = dest
			distribution = origDistribution


		return source


		
	def GetGeometryList(self, configSection):
		"""
		Creates a list of geometry infos for each basos on the
		the given potential (configSection)
		"""
		geometryList = []
		for i, basis in enumerate(self.BasisList):
			geometryName = configSection.Get("geometry%i" % i)
			geometryInfo = basis.GetGeometryInfo(geometryName)
			geometryList.append(geometryInfo)

		return geometryList


	def SetupRepresentation(self, geometryList):
		"""
		Creates a CombinedRepresentation with subrepresentations to match the
		current potential
		"""

		#Create Combined Representation
		repr = CreateInstanceRank("core.CombinedRepresentation", self.Rank)
	
		#Set distributed model
		distrib = CreateDistribution(None, rank=self.Rank)
		repr.SetDistributedModel(distrib)
	
		#Setup sub representations
		for i, basis in enumerate(self.BasisList):
			geometryInfo = geometryList[i]
	
			#Create representation
			if geometryInfo.UseGridRepresentation():
				subRepr = basis.GetGridRepresentation()
			else:
				subRepr = basis.GetBasisRepresentation()
		
			#Set distributed model
			subDistrib = distrib.CreateSubDistributedModel()
			subRepr.SetDistributedModel(subDistrib)
			subRepr.SetBaseRank(i)
			repr.SetRepresentation(i, subRepr)


		return repr

	def SetupPotentialEvaluator(self, configSection, psi):
		evaluatorPrefix = "core.DynamicPotentialEvaluator"
		classname = configSection.classname
		potentialEvaluator = CreatePotentialInstance(classname, self.Rank, evaluatorPrefix)
		configSection.Apply(potentialEvaluator)
		if hasattr(potentialEvaluator, "Setup"):
			potentialEvaluator.Setup(psi)
			
		return potentialEvaluator
		

#------------------------------------------------------------------------------------
#                       Propagation classes
#------------------------------------------------------------------------------------

from pyprop import PotentialWrapper

class TensorPotential(PotentialWrapper):
	"""
	Potential wrapper for TensorPotentials. See PotentialWrapper for more information on the
	PotentialWrapper interface

	A TensorPotential is potential expressed in the basis functions of the wavefunction.
	This makes it possible to apply the potential without transforming the wavefunction to the
	grid. The downside is that potentials are usually diagonal in the grid basis, and not in the
	basisfunction basis. This means that we must do a matrix vector product to apply the potential.

	in 1D
	out(i) = sum_j V(i, j) psi(j) 
	in 2D
	out(i', j') = sum_{i,j} V(i', i, j', j) psi(i', j')

	in the tensor potential V has the same rank as the wavefunction, where i',i and j',j are compressed
	into one rank each. This is done in order to easier be able to exploit symmetries in the potentials and
	bandedness in the basises (such as the dipole selection rule V = 0 for l' != l +/- 1).

	
	"""

	def __init__(self, psi):
		self.GeometryList = None
		self.PotentialData = None
		self.Name = None
		self.psi = psi
		self.MultiplyAlgorithm = 4

	def ApplyConfigSection(self, configSection):
		#Check wheter this is a time dependent potential
		self.IsTimeDependent = False
		if hasattr(configSection, "time_function"):
			self.IsTimeDependent = True
			self.TimeFunction = lambda t: configSection.time_function(configSection, t)
			self.OriginalTimeFunction = configSection.time_function

	def SetupStep(self, timestep):
		self.BasisPairs = [geom.GetBasisPairs() for geom in self.GeometryList]

		multiplyFuncName = "pyprop.core.TensorPotentialMultiply_" + "_".join([geom.GetStorageId() for geom in self.GeometryList])
		self.MultiplyFunction = eval(multiplyFuncName)
		
	def AdvanceStep(self, t, timestep):
		raise NotImplementedException("TensorPotentials can not be exponentiated directly")

	def MultiplyPotential(self, srcPsi, destPsi, t, timestep):
		rank = srcPsi.GetRank()
		
		source = srcPsi.GetData()
		dest = destPsi.GetData()

		timeScaling = 1.0
		if self.IsTimeDependent:
			timeScaling = self.TimeFunction(t)
		
		#TODO: Implement support for parallelization. 

		#Construct argument list
		#Default parameters
		argList = [self.PotentialData, timeScaling, source, dest]
		#Parameters for each storage
		for i, geom in enumerate(self.GeometryList):
			argList += geom.GetMultiplyArguments()

		#Perform multiplication
		self.MultiplyFunction(*argList)

	#def GetExpectationValue(self, tmpPsi, t, timeStep):
	#	self.GetExpectationValue(self.psi, tmpPsi, t, timeStep)
	
	def GetExpectationValue(self, tmpPsi, t, timeStep):
		tmpPsi.Clear()
		self.MultiplyPotential(self.psi, tmpPsi, t, timeStep)

		#Solve for all overlap matrices
		repr = self.psi.GetRepresentation()
		repr.SolveOverlap(tmpPsi)
		
		return self.psi.InnerProduct(tmpPsi)


from pyprop import PropagatorBase

class BasisPropagator(PropagatorBase):
	"""
	Propagator that does not transform the wavefunction to the grid basis, but rather applies the potentials
	(TensorPotentials) directly in the basis.

	BasisPropagator assumes that the wavefunction has a CombinedRepresentation.

	AdvanceStep is not supported by BasisPropagator. Rather, one should use a propagator utilizing MultplyHamiltonian
	such as the ODE and Krylov propagator.
	"""

	__Base = PropagatorBase
	
	def __init__(self, psi):
		self.__Base.__init__(self, psi)
		self.Rank = psi.GetRank()

		#Create a tensor potential generator 		
		#We will use this to construct tensor potentials from potentials specified on the grid
		self.TensorPotentialGenerator = TensorPotentialGenerator(representation = psi.GetRepresentation())

		#We need the temp array for solving overlap matrix eqns
		self.TempPsi2 = self.psi.Copy()

	def ApplyConfig(self, config):
		#Create any precalculated potentials 
		self.__Base.ApplyConfig(self, config)

		self.GeneratePotentials(config)
		self.ConsolidatePotentials()


	def GeneratePotentials(self, config):
		"""
		Genereate TensorPotentials from potentials specified on the grid in the
		configuration file
		"""

		#Potentials we should create on the fly
		if hasattr(config.Propagation, "grid_potential_list"):
			potentials = config.Propagation.grid_potential_list
			generator = self.TensorPotentialGenerator

			for potentialName in potentials:
				#Find the corresponding config section
				configSection = config.GetSection(potentialName)

				#Use TensorPotentialGenerator to construct potential in basis
				geometryList = generator.GetGeometryList(configSection)
				potentialData = generator.GeneratePotential(configSection)
				originalPotential = getattr(generator, "OriginalPotential", None)

				#Create PotentialWrapper for TensorPotential
				potential = TensorPotential(self.psi)
				configSection.Apply(potential)
				potential.GeometryList = geometryList
				potential.PotentialData = potentialData
				potential.OriginalPotential = originalPotential
				potential.Name = potentialName

				#Add potential to potential list
				self.PotentialList.append(potential)


	def ConsolidatePotentials(self):
		"""
		Try to consolidate potentials having the same geometry into
		one potential to save evaulations. Potentials containing a time
		dependent part needs to be treated separately, and thus is not considered
		"""
		
		PrintOut( "Consolidating similar potentials: " )
		PrintOut( "Starting with potentials:" )
		for pot in self.PotentialList:
			PrintOut( "    %s" % pot.Name )

		#only non timedependent potentials are considered
		potentials = list(self.PotentialList) #[pot for pot in self.PotentialList if not pot.IsTimeDependent]
		removePotentials = []

		#Consolidate til we're at the last potential
		i = 0
		while(i < len(potentials)-1):
			curPot = potentials[i]

			#Loop over all potentials after curPot
			for otherPot in list(potentials[i+1:]):
				#We can consolidate curPot and otherPot if all the index pairs are the same,
				#And curPot and otherPot has the same time dependency
				canConsolidate = True

				#If one of the potentials is time dependent the other must also be
				if canConsolidate and curPot.IsTimeDependent != otherPot.IsTimeDependent:
					canConsolidate = False

				#If they are both time dependent, they must have the same time function
				if canConsolidate and curPot.IsTimeDependent and otherPot.IsTimeDependent:
					if otherPot.OriginalTimeFunction != curPot.OriginalTimeFunction:
						canConsolidate = False
			
				#Both potentials must have the same storage 
				if canConsolidate:
					for rank in range(self.Rank):
						if curPot.GeometryList[rank].GetStorageId() != otherPot.GeometryList[rank].GetStorageId():
							canConsolidate = False
							break
						else:
							#Potentials must have the same basis pairs
							curBasisPairs = asarray(curPot.GeometryList[rank].GetBasisPairs())
							otherBasisPairs = asarray(otherPot.GeometryList[rank].GetBasisPairs())
							
							if len(curBasisPairs) != len(otherBasisPairs):
								canConsolidate = False
								break
							elif not (curBasisPairs == otherBasisPairs).all():
								canConsolidate = False
								break

				#Add otherPot to curPot
				if canConsolidate:
					curPot.PotentialData[:] += otherPot.PotentialData
					curPot.Name += "+" + otherPot.Name
					potentials.remove(otherPot)
					removePotentials.append(otherPot)

			#Next potential
			i += 1

		#Remove consolidated potentials from the main potential list
		for pot in removePotentials:
			self.PotentialList.remove(pot)

		PrintOut( "Ended up with potentials:" )
		for pot in self.PotentialList:
			PrintOut( "    %s" % pot.Name )


	def ApplyConfigSection(self, configSection): 
		self.__Base.ApplyConfigSection(self, configSection)

	def SetupStep(self, dt):
		self.__Base.SetupStep(self, dt)

	def MultiplyHamiltonian(self, srcPsi, destPsi, t, dt):
		#Multiply potentials
		destPsi.GetData()[:] = 0
		self.MultiplyPotential(srcPsi, destPsi, t, dt)

		#Solve for all overlap matrices
		repr = srcPsi.GetRepresentation()
		if not all([repr.IsOrthogonalBasis(i) for i in range(self.Rank)]):
			repr.SolveOverlap(destPsi)

	def MultiplyHamiltonianBalancedOverlap(self, srcPsi, destPsi, t, dt):
		#Store input psi
		def StorePsi():
			self.TempPsi2.GetData()[:] = srcPsi.GetData()

		#Solve for all overlap matrices
		def SolveOverlap1():
			repr = srcPsi.GetRepresentation()
			repr.SolveSqrtOverlap(False, srcPsi)
		
		#Multiply potentials
		def MultiplyPotential():
			destPsi.GetData()[:] = 0
			self.MultiplyPotential(srcPsi, destPsi, t, dt)

		#Solve for all overlap matrices
		def SolveOverlap2():
			repr = destPsi.GetRepresentation()
			repr.SolveSqrtOverlap(True, destPsi)

		#Restore input psi back to its original state
		def RestorePsi():
			srcPsi.GetData()[:] = self.TempPsi2.GetData()

		StorePsi()
		SolveOverlap1()
		MultiplyPotential()
		SolveOverlap2()
		RestorePsi()

	def AdvanceStep(self, t, dt):
		raise NotImplementedException("BasisPropagator does not support AdvanceStep. Use it as a base for explicit propagators")

	def GetBasisFunction(self, rank, basisIndex):
		raise NotImplementedException("Implement GetBasisFunction...")



