import re
from numpy import ones
from pyprop import GetClassLogger, core


class GeometryInfoFiniteDifferenceBandedDistributed(GeometryInfoDistributedBase):
	"""
	Geometry for a banded finite difference matrix with support for parallel
	matvec through Epetra.
	"""
	def __init__(self, repr, bandCount, useGrid):
		GeometryInfoDistributedBase.__init__(self)
		#Set member variables 
		self.UseGrid = useGrid
		self.BaseRank = repr.GetBaseRank()
		self.Representation = repr
		self.BandCount = bandCount

	def UseGridRepresentation(self):
		return self.UseGrid

	def GetGlobalBasisPairCount(self):
		N = int(self.Representation.GetFullShape()[0])
		k = self.BandCount
		return (2 * k + 1) * N - (k*(k+1))

	def SetupBasisPairs(self):
		distr = self.Representation.GetDistributedModel()
		rank = self.Representation.GetBaseRank()
		localRange = distr.GetLocalIndexRange(self.GetGlobalBasisPairCount(), rank)

		count = self.GetGlobalBasisPairCount()

		N = int(self.Representation.GetFullShape()[0])
		k = self.BandCount
		pairs = zeros((count, 2), dtype=int32)
		pairs[:,0] = 0
		pairs[:,1] = N-1

		index = 0
		for i in xrange(N):
			for j in xrange(max(0, i-k), min(i+k+1, N)):
				pairs[index, 0] = i
				pairs[index, 1] = j
				index+=1

		if index != pairs.shape[0]:
			print "Fullshape = %s" % N
			print pairs
			raise Exception("All basis pairs where not assigned! index \
					!= pairs.shape[0] (%s, %s)" % (index, pairs.shape[0]))

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
		self.BoundaryScaling = 1.0
		self.Logger = GetClassLogger(self)

	def SetupBasis(self, repr):
		self.Representation = repr
		baseRank = repr.GetBaseRank()
		self.GridSize = len(self.Representation.GetGlobalGrid(baseRank))
		self.DifferenceOrder = 1
		self.BandCount = (self.DifferenceOrder - 1) / 2
		self.BoundaryScaling = 1.0
		self.Logger = GetClassLogger(self)
		
	def GetGridRepresentation(self):
		return self.Representation

	def GetBasisRepresentation(self):
		return self.Representation
		
	def GetGeometryInfo(self, geometryName):
		geom = geometryName.lower().strip()
		if geom == "identity":
			return GeometryInfoCommonIdentity(True)
		elif geom == "diagonal":
			return GeometryInfoCommonDiagonal(self.Representation, True)
		else:
			diffOrderSearch =  re.search("-\d+", geom)
			if not diffOrderSearch:
				raise  UnsupportedGeometryException("BasisfunctionFiniteDifference requires specification of difference order, you said %s" % geometryName)

			self.DifferenceOrder =  eval(diffOrderSearch.group()[1:])
			self.BandCount = (self.DifferenceOrder - 1) / 2
			if geom.startswith == "dense":
				return GeometryInfoCommonDense(self.GridSize, True)
			elif re.search("bandeddistributed", geom):
				#return GeometryInfoCommonBandedDistributed(self.Representation, self.BandCount, True)
				return GeometryInfoFiniteDifferenceBandedDistributed(self.Representation, self.BandCount, True)
			elif re.search("banded-nonhermitian", geom) or re.search("banded", geom):
				return GeometryInfoCommonBandedNonHermitian(self.GridSize, self.BandCount, True)
			else:
				raise UnsupportedGeometryException("Geometry '%s' not supported by BasisfunctionFiniteDifference" % geometryName)

	def RepresentPotentialInBasis(self, source, dest, rank, geometryInfo, \
			differentiation, configSection):
		if differentiation == 0:
			diffMatrix = ones((self.GridSize, 1), dtype=complex)
		elif differentiation == 2:
			self.BoundaryScaling = getattr(configSection, "boundary_scaling%i" % \
				rank, self.BoundaryScaling)
			self.Logger.debug("Using boundary condition scaling for rank %i: %s" %\
				(rank, self.BoundaryScaling))
			fd = core.FiniteDifferenceHelperCustomBoundary()
			grid = self.Representation.GetGlobalGrid(self.Representation.GetBaseRank())
			fd.Setup(grid, self.DifferenceOrder, self.BoundaryScaling)
			diffMatrix = fd.SetupLaplacianBlasBanded().copy()
		else:
			raise Exception("Finite Difference currently only supports diff of order 2")

		indexPairs = geometryInfo.GetGlobalBasisPairs()
		core.RepresentPotentialInBasisFiniteDifference(diffMatrix, source, dest, indexPairs, rank)

