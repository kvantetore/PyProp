

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
		self.Rank = psi.GetRank()
		self.OriginalPotential = None

	def ApplyConfigSection(self, configSection):
		self.DebugPotential = False
		if hasattr(configSection, "debug_potential"):
			self.DebugPotential = configSection.debug_potential

		#Check wheter this is a time dependent potential
		self.IsTimeDependent = False
		if hasattr(configSection, "time_function"):
			self.IsTimeDependent = True
			self.TimeFunction = lambda t: configSection.time_function(configSection, t)
			self.OriginalTimeFunction = configSection.time_function


		#Use TensorPotentialGenerator to generate potential data
		generator = TensorPotentialGenerator(representation=self.psi.GetRepresentation())
		geometryList = generator.GetGeometryList(configSection)
		self.GeometryList = geometryList
		self.Name = configSection.name

		#If filename is specified, load it from disk, otherwise, generate it
		if hasattr(configSection, "filename"):
			filename = configSection.filename
			datasetPath = configSection.dataset
			distr = self.psi.GetRepresentation().GetDistributedModel()

			#setup data buffer
			potShape = [geomInfo.GetBasisPairCount() for geomInfo in geometryList]
			dataBuffer = CreateInstanceRank("core.DataBuffer", self.Rank)
			dataBuffer.ResizeArray(array(potShape))
			self.PotentialDataBuffer = dataBuffer
			self.PotentialData = dataBuffer.GetArray()

			#Load from disk
			serialization.LoadTensorPotential(filename, datasetPath, self, distr)

		else:
			potentialDataBuffer = generator.GeneratePotential(configSection)
			self.PotentialDataBuffer = potentialDataBuffer
			self.PotentialData = potentialDataBuffer.GetArray()
			self.OriginalPotential = getattr(generator, "OriginalPotential", None)
		

	def SetupStep(self, timestep):
		self.BasisPairs = [geom.GetBasisPairs() for geom in self.GeometryList]

		multiplyFuncName = "core.TensorPotentialMultiply_" + "_".join([geom.GetStorageId() for geom in self.GeometryList])
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
			argList += geom.GetMultiplyArguments(srcPsi)

		#Perform multiplication
		self.MultiplyFunction(*argList)

	#def GetExpectationValue(self, tmpPsi, t, timeStep):
	#	self.GetExpectationValue(self.psi, tmpPsi, t, timeStep)
	
	def GetExpectationValue(self, psi, tmpPsi, t, timeStep):
		tmpPsi.Clear()
		self.MultiplyPotential(psi, tmpPsi, t, timeStep)

		#Solve for all overlap matrices
		repr = self.psi.GetRepresentation()
		repr.SolveOverlap(tmpPsi)
		
		return self.psi.InnerProduct(tmpPsi)


	def CanConsolidate(self, otherPot):
		"""
		Checks wheter this potential can be consolidated with otherPot
		"""
		#We can consolidate self and otherPot if all the index pairs are the same,
		#And self and otherPot has the same time dependency
		canConsolidate = True

		#If one of the potentials is time dependent the other must also be
		if canConsolidate and self.IsTimeDependent != otherPot.IsTimeDependent:
			canConsolidate = False

		#If they are both time dependent, they must have the same time function
		if canConsolidate and self.IsTimeDependent and otherPot.IsTimeDependent:
			if otherPot.OriginalTimeFunction != self.OriginalTimeFunction:
				canConsolidate = False

		#don't consolidate debug potentials
		if canConsolidate and (self.DebugPotential or otherPot.DebugPotential):
			canConsolidate = False
			
		#Both potentials must have the same storage 
		if canConsolidate:
			for rank in range(self.Rank):
				if self.GeometryList[rank].GetStorageId() != otherPot.GeometryList[rank].GetStorageId():
					canConsolidate = False
					break

		if canConsolidate:
			for rank in range(self.Rank):
				if any(self.GeometryList[rank].GetBasisPairs() != otherPot.GeometryList[rank].GetBasisPairs()):
					canConsolidate = False
					break


		return canConsolidate		
