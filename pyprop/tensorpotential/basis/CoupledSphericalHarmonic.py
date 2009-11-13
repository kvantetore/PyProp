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
	
	def GetGlobalBasisPairCount(self):
		return self.GetBasisPairs().shape[0]
		
	def GetBasisPairs(self):
		if self.IndexPairs == None:
			self.IndexPairs = self.SelectionRule.GetBasisPairs(self.Representation)
		return self.IndexPairs

	def GetStorageId(self):
		return "Simp"

	def GetMultiplyArguments(self, psi):
		return [self.GetBasisPairs()]


class GeometryInfoCoupledSphericalHarmonicDistributed(GeometryInfoDistributedBase):
	"""
	Geometry information for coupled spherical harmonic geometries using distributed matvec.
	The CoupledSphericalHarmonicRepresentation and a 
	CoupledSphericalSelectionRule is specified.

	"""
	def __init__(self, representation, selectionRule):
		GeometryInfoDistributedBase.__init__(self)

		self.SelectionRule = selectionRule
		self.Representation = representation
		self.StepArguments = None
		self.TempArrays = None
		self.MultiplyArguments = None

	def UseGridRepresentation(self):
		return False
	
	def GetStorageId(self):
		return "SimpD"

	def GetMultiplyArguments(self, psi):
		if self.MultiplyArguments == None:
			if self.StepArguments == None:
				self.SetupLocalBasisPairs()
			if self.TempArrays == None:
				self.SetupTempArrays(psi)
			self.MultiplyArguments = self.StepArguments + self.TempArrays

		return self.StepArguments + self.TempArrays
	
	def SetupBasisPairs(self):
		indexPairs = self.SelectionRule.GetBasisPairs(self.Representation)
		self.GlobalIndexPairs = array(indexPairs)
		distrib = self.Representation.GetDistributedModel()
		rank = self.Representation.GetBaseRank()
		globalSize = self.Representation.GetFullShape()[0]
	
		distribIndexList = SetupDistributedIndexList(globalSize, indexPairs, distrib, rank)
		self.LocalBasisPairIndices = array(distribIndexList[distrib.ProcId])
		stepList = SetupStepList(globalSize, indexPairs, distribIndexList, distrib, rank)
		self.MaxRecvCount = max([len(step.RecvProcList) for step in stepList])

		self.StepArguments = [int(globalSize)] + list(StepListToArray(stepList))
		self.LocalIndexPairs = array([[step.GlobalRow, step.GlobalCol] for step in stepList if step.LocalMatrixIndex!=-1], dtype=int32)

	def SetupTempArrays(self, psi):
		dataShape = psi.GetData().shape
		rank = self.Representation.GetBaseRank()

		recvTempShape = list(dataShape)
		recvTempShape[rank] = self.MaxRecvCount
		recvTemp = zeros(recvTempShape, dtype=complex)

		sendTempShape = list(dataShape)
		sendTempShape[rank] = 2
		sendTemp = zeros(sendTempShape, dtype=complex)

		self.TempArrays = [recvTemp, sendTemp]


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

		#Check if rank is distributed; if so, set appropriate Geometry function
		baseRank = self.BasisRepresentation.GetBaseRank()
		distr = self.BasisRepresentation.GetDistributedModel()
		self.GeometryFunction = GeometryInfoCoupledSphericalHarmonic
		if distr.IsDistributedRank(baseRank):
			self.GeometryFunction = GeometryInfoCoupledSphericalHarmonicDistributed

	def GetGridRepresentation(self):
		raise Exception("Coupled Spherical Harmonics does not support grid representations")

	def GetBasisRepresentation(self):
		return self.BasisRepresentation
		
	def GetGeometryInfo(self, geometryName):
		geom = geometryName.lower().strip()
		if geom == "identity":
			return GeometryInfoCommonIdentity(False)
		elif geom == "diagonal":
			selectionRule = core.CoupledSphericalSelectionRuleDiagonal()
			return self.GeometryFunction(self.BasisRepresentation, selectionRule)
		elif geom == "dense":
			return GeometryInfoCommonDense(self.BasisSize, False)
		elif geom.startswith("selectionrule_"):
			selectionRuleName = geom[len("selectionrule_"):]
			if selectionRuleName == "r12":
				selectionRule = core.CoupledSphericalSelectionRuleR12Old()
			elif selectionRuleName[:3] == "r12":
				multipoleCutoff = int(selectionRuleName[4:])
				selectionRule = core.CoupledSphericalSelectionRuleR12(multipoleCutoff)
			elif selectionRuleName == "linearpolarizedfield":
				selectionRule = core.CoupledSphericalSelectionRuleLinearPolarizedField()
			else:	
				raise Exception("Unkonwn selection rule %s" % selectionRuleName)

			return self.GeometryFunction(self.BasisRepresentation, selectionRule)
		else:
			raise UnsupportedGeometryException("Geometry '%s' not supported by BasisfunctionReducedSpherical" % geometryName)

	def RepresentPotentialInBasis(self, source, dest, rank, geometryInfo, differentiation):
		raise Exception("Coupled Spherical Harmonics are already in a basis and should not be integrated")

