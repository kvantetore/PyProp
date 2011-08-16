from numpy import ones, array, zeros, int32, r_
import re

import pyprop.core as core
import pyprop.tensorpotential.geometryinfo as geometryinfo
from pyprop.tensorpotential.tensorpotentialgenerator import RegisterBasisfunction
from pyprop.pyproplogging import GetClassLogger

import libfinitedifference

class GeometryInfoFiniteDifferenceBandedDistributed(geometryinfo.GeometryInfoDistributedBase):
	"""
	Geometry for a banded finite difference matrix with support for parallel
	matvec through Epetra.
	"""
	def __init__(self, repr, bandCount, useGrid):
		geometryinfo.GeometryInfoDistributedBase.__init__(self)
		#Set member variables
		self.UseGrid = useGrid
		self.BaseRank = repr.GetBaseRank()
		self.Representation = repr
		self.BandCount = bandCount

	def UseGridRepresentation(self):
		return self.UseGrid

	def GetGlobalBasisPairCount(self):
		"""Return number of global basis pairs

		Padding is added at start and beginning of matrix to simpliy
		parallelization, that is

		  N_global = numberOfRows * rowSize

		"""
		numberOfRows = int(self.Representation.GetFullShape()[0])
		k = self.BandCount
		rowSize = 2*k + 1
		return numberOfRows * rowSize

	def SetupBasisPairs(self):
		"""Setup global and local basis pairs

		"""

		distr = self.Representation.GetDistributedModel()
		rank = self.Representation.GetBaseRank()

		#Get rowsize and colsize
		numberOfRows = int(self.Representation.GetFullShape()[0])
		k = self.BandCount
		rowSize = 2*k + 1

		#Get local rows from wavefunction distribution
		psiRange = distr.GetLocalIndexRange(numberOfRows, rank)

		#Check that proc number for this rank divides number of rows
		if distr.IsDistributedRank(rank):
			procRankIdx = distr.GetDistribution()[rank]
			numProcRank = distr.GetTranspose().GetProcGridShape()[procRankIdx]
			infoStr = """
				Number of procs %i for rank %i does not divide number of
				matrix rows (psi shape) %s
				""".replace("\n", "").replace("\t", "") \
				% (numProcRank, rank, numberOfRows)
			assert (numberOfRows % numProcRank == 0), infoStr

		index = 0
		count = self.GetGlobalBasisPairCount()
		pairs = zeros((count, 2), dtype=int32)
		localIdx = []
		for i in xrange(numberOfRows):
			#compute start and end column indices
			colStart = max(0, i-k)
			colEnd = min(i+k+1, numberOfRows)

			#iterate one full row
			for j in xrange(colStart, colStart+rowSize):

				#Check for padded element, if so, set row/col to -1
				if j >= colEnd:
					pairs[index, 0] = -1
					pairs[index, 1] = -1
				else:
					pairs[index, 0] = i
					pairs[index, 1] = j

				#get index if current row is local
				if i in r_[psiRange]:
					localIdx += [index]
				index+=1

		#check that all pairs where assigned
		if index != pairs.shape[0]:
			print "Fullshape = %s" % numberOfRows
			print pairs
			raise Exception("All basis pairs where not assigned! index \
					!= pairs.shape[0] (%s, %s)" % (index, pairs.shape[0]))

		#store global and local index pairs
		self.GlobalIndexPairs = pairs
		self.LocalBasisPairIndices = array(localIdx, dtype=int32)
		self.LocalIndexPairs = pairs[localIdx]


	def GetStorageId(self):
		return "BandDistr"

	def GetMultiplyArguments(self, psi):
		raise Exception("GetMultiplyArguments not implemented, use Epetra for \
			matvec!")

	def SetupTempArrays(self, psi):
		raise Exception("SetupTempArrays not implemented, use Epetra for \
			matvec!")


#-------------------------------------------------------------------------------
#                       Finite Difference Basis Function
#-------------------------------------------------------------------------------

@RegisterBasisfunction(core.CustomGridRepresentation)
class BasisfunctionFiniteDifference(geometryinfo.BasisfunctionBase):
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
		self.BoundaryScaling = \
				array([.0] * (self.BandCount * (self.BandCount+1) / 2))
		self.Offset = 0
		self.Logger = GetClassLogger(self)

	def SetupBasis(self, repr):
		self.Representation = repr
		baseRank = repr.GetBaseRank()
		self.GridSize = len(self.Representation.GetGlobalGrid(baseRank))
		self.DifferenceOrder = 1
		self.BandCount = (self.DifferenceOrder - 1) / 2
		self.BoundaryScaling = \
				array([.0] * (self.BandCount * (self.BandCount+1) / 2))
		self.Offset = 0
		self.Logger = GetClassLogger(self)

	def GetGridRepresentation(self):
		return self.Representation

	def GetBasisRepresentation(self):
		return self.Representation

	def GetGeometryInfo(self, geometryName):
		geom = geometryName.lower().strip()
		if geom == "identity":
			return geometryinfo.GeometryInfoCommonIdentity(True)
		elif geom == "diagonal":
			return geometryinfo.GeometryInfoCommonDiagonal(self.Representation, True)
		else:
			diffOrderSearch =  re.search("-\d+", geom)
			if not diffOrderSearch:
				raise  geometryinfo.UnsupportedGeometryException("BasisfunctionFiniteDifference requires specification of difference order, you said %s" % geometryName)

			self.DifferenceOrder =  eval(diffOrderSearch.group()[1:])
			self.BandCount = (self.DifferenceOrder - 1) / 2
			if geom.startswith == "dense":
				return geometryinfo.GeometryInfoCommonDense(self.GridSize, True)
			elif re.search("banded-nonhermitian", geom) or re.search("banded", geom):
				return geometryinfo.GeometryInfoCommonBandedNonHermitian(self.GridSize, self.BandCount, True)
			elif re.search("bandeddistributed", geom):
				return geometryinfo.GeometryInfoCommonBandedDistributed(self.Representation, self.BandCount, True)
			else:
				raise geometryinfo.UnsupportedGeometryException("Geometry '%s' not supported by BasisfunctionFiniteDifference" % geometryName)

	def RepresentPotentialInBasis(self, source, dest, rank,
	                              geometryInfo, differentiation, configSection):
		if differentiation == 0:
			diffMatrix = ones((self.GridSize, 1), dtype=complex)
		elif differentiation == 2:
			customBoundary = getattr(configSection, "custom_boundary", False)
			if customBoundary:
				self.BoundaryScaling = getattr(configSection, "boundary_scaling%i"\
					%  rank, self.BoundaryScaling)
				self.Offset = getattr(configSection, "offset%i"\
					%  rank, self.Offset)
				self.Logger.debug("Using boundary condition scaling for rank %i: \
						%s" % (rank, self.BoundaryScaling))
				fd = libfinitedifference.FiniteDifferenceHelperCustomBoundary()
				baseRank = self.Representation.GetBaseRank()
				grid = self.Representation.GetGlobalGrid(baseRank)
				fd.Setup(grid, self.DifferenceOrder, self.BoundaryScaling, \
						self.Offset)
			else:
				fd = libfinitedifference.FiniteDifferenceHelper()
				baseRank = self.Representation.GetBaseRank()
				grid = self.Representation.GetGlobalGrid(baseRank)
				fd.Setup(grid, self.DifferenceOrder)
			diffMatrix = fd.SetupLaplacianBlasBanded().copy()

		else:
			raise Exception("Finite Difference currently only supports diff \
				of order 2")

		indexPairs = geometryInfo.GetGlobalBasisPairs()
		core.RepresentPotentialInBasisFiniteDifference(diffMatrix, source, dest, indexPairs, rank)
