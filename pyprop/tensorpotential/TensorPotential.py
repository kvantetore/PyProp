

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



