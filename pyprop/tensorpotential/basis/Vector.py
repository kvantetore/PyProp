import pyprop.tensorpotential.GeometryInfo as geometryinfo

#------------------------------------------------------------------------------------
#                       Vector
#------------------------------------------------------------------------------------

class BasisfunctionVector(geometryinfo.BasisfunctionBase):
	"""
	Basisfunction class for Vector
	see BasisfunctionBase for details
	"""

	def ApplyConfigSection(self, configSection):
		self.SetupBasis(configSection)
		
	def SetupBasis(self, basisRepresentation):
		self.BasisRepresentation = basisRepresentation
		self.BasisSize = self.BasisRepresentation.VectorSize

	def GetGridRepresentation(self):
		raise Exception("Vector basis does not support grid representations")

	def GetBasisRepresentation(self):
		return self.BasisRepresentation
		
	def GetGeometryInfo(self, geometryName):
		geom = geometryName.lower().strip()

		if geom == "identity":
			return geometryinfo.GeometryInfoCommonIdentity(False)
		elif geom == "diagonal":
			return geometryinfo.GeometryInfoCommonDiagonal(self.BasisRepresentation, False)
		elif geom.startswith("banded-nonhermitian"):
			BandCount = int(geom.split("-")[2])
			return geometryinfo.GeometryInfoCommonBandedNonHermitian(self.BasisSize, BandCount, False)
		elif geom == "dense":
			return geometryinfo.GeometryInfoCommonDense(self.BasisSize, False)
		elif geom == "hermitian":
			return geometryinfo.GeometryInfoCommonDense(self.BasisSize, False)
		else:
			raise geometryinfo.UnsupportedGeometryException("Geometry '%s' not supported by BasisfunctionVector" % geometryName)

	def RepresentPotentialInBasis(self, source, dest, rank, geometryInfo, differentiation):
		raise Exception("Vector basis already in a basis and should not be integrated")

