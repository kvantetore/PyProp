
class EpetraPotential(PotentialWrapper):
	"""
	Potential wrapper for EpetraPotential. See PotentialWrapper for more information on the
	PotentialWrapper interface

	"""

	def __init__(self, psi):
		self.Name = None
		self.Psi = psi
		self.Rank = psi.GetRank()
		self.Potential = CreateInstanceRank("core.EpetraPotential", self.Rank)

	def ApplyConfigSection(self, configSection):
		#Setup Epetra potential
		self.Potential.Setup(self.Psi)

		self.Name = configSection.name

		#Check wheter this is a time dependent potential
		self.IsTimeDependent = False
		if hasattr(configSection, "time_function"):
			self.IsTimeDependent = True
			self.TimeFunction = lambda t: configSection.time_function(configSection, t)
			self.OriginalTimeFunction = configSection.time_function

	
	def MultiplyPotential(self, srcPsi, destPsi, t, timestep):
		self.Potential.Multiply(srcPsi, destPsi)

		if self.IsTimeDependent:
			timeDepFactor = self.TimeFunction(t)
			destPsi.GetData()[:] *= timeDepFactor

	
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
		#We can consolidate self and otherPot if they the same time dependency,
		#and if self is not Filled.
		canConsolidate = True

		#If FECrsMatrix is filled, we cannot consolidate
		if canConsolidate and self.Potential.Filled():
			canConsolidate = False

		#If one of the potentials is time dependent the other must also be
		if canConsolidate and self.IsTimeDependent != otherPot.IsTimeDependent:
			canConsolidate = False

		#If they are both time dependent, they must have the same time function
		if canConsolidate and self.IsTimeDependent and otherPot.IsTimeDependent:
			if otherPot.OriginalTimeFunction != self.OriginalTimeFunction:
				canConsolidate = False

		return canConsolidate


	def AddTensorPotentialData(self, data, basisPairs, cutOff):
		"""
		Add data in array 'data' to own Epetra matrix
		"""
		self.Potential.AddTensorPotentialData(data, basisPairs, cutOff)


	def GlobalAssemble(self):
		"""
		Assemble matrix elements across procs, and signal FillComplete
		"""
		self.Potential.GlobalAssemble()
