from numpy import ones
import re

import pyprop.core as core
import pyprop.tensorpotential.GeometryInfo as geometryinfo


#------------------------------------------------------------------------------------
#                       Finite Difference Basis Function
#------------------------------------------------------------------------------------


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

	def RepresentPotentialInBasis(self, source, dest, rank, geometryInfo, differentiation):
		if differentiation == 0:
			diffMatrix = ones((self.GridSize, 1), dtype=complex)
		elif differentiation == 2:
			fd = core.FiniteDifferenceHelper()
			grid = self.Representation.GetGlobalGrid(self.Representation.GetBaseRank())
			fd.Setup(grid, self.DifferenceOrder)
			diffMatrix = fd.SetupLaplacianBlasBanded().copy()
		else:
			raise Exception("Finite Difference currently only supports diff of order 2")

		indexPairs = geometryInfo.GetGlobalBasisPairs()
		core.RepresentPotentialInBasisFiniteDifference(diffMatrix, source, dest, indexPairs, rank)




