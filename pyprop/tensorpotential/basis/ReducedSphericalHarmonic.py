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
		self.SphericalHarmonicObject = core.ReducedSphericalTools()
		self.SphericalHarmonicObject.Initialize(self.LMax)

	def GetGridRepresentation(self):
		repr = core.ThetaRepresentation()
		repr.SetupRepresentation(self.LMax)
		return repr

	def GetBasisRepresentation(self):
		repr = core.ReducedSphericalHarmonicRepresentation()
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
		core.RepresentPotentialInBasisReducedSphericalHarmonic(self.SphericalHarmonicObject, source, dest, pairs, storageId, rank, differentiation)


