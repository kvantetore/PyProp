class TwoElectronPreconditioner:
	"""
	Preconditioner for iterative solvers of linear equations,

	(1)     (a H + b S) x = y  ->  x

	where H is the Hamiltonian, S is the overlap matrix
	and a and b are complex scaling factors (scalingH and scalingS).

	The preconditioner is given a set of tensor potentials which
	should approximate the Hamiltonian as well as 
	possible, as the system (1) is solved exactly for the given
	potentials
	"""

	def __init__(self, psi):
		self.Rank = psi.GetRank()
		self.psi = psi

	def ApplyConfigSection(self, conf):
		self.OverlapSection = conf.Config.GetSection("OverlapMatrixPotential")
		self.PotentialSections = [conf.Config.GetSection(s) for s in conf.potential_evaluation]

	def SetHamiltonianScaling(self, scalingH):
		self.HamiltonianScaling = scalingH

	def SetOverlapScaling(self, scalingS):
		self.OverlapScaling = scalingS

	def GetHamiltonianScaling(self):
		return self.HamiltonianScaling

	def GetOverlapScaling(self):
		return self.OverlapScaling

	def Setup(self, prop):
		"""
		Set up a tensor potential for overlap potential and all other potentials
		and add them together, assuming they have the same layout
		ending up with a potential containing S + scalingH * (P1 + P2 + ...)
		"""

		#Setup overlap potential
		pyprop.PrintMemoryUsage("Before Preconditioner Generate Potential (Overlap)")
		tensorPotential = prop.BasePropagator.GeneratePotential(self.OverlapSection)
		tensorPotential.PotentialData *= self.GetOverlapScaling()

		#Add all potentials to solver
		scalingH = self.GetHamiltonianScaling()
		for conf in self.PotentialSections:
			pyprop.PrintMemoryUsage("Before Preconditioner Generate Potential (%s)" % conf)
			#Setup potential in basis
			potential = prop.BasePropagator.GeneratePotential(conf)
			if not tensorPotential.CanConsolidate(potential):
				raise Exception("Cannot consolidate potential %s with overlap-potential" % (potential.Name))
		
			#Add potential
			potential.PotentialData *= scalingH
			tensorPotential.PotentialData += potential.PotentialData
			del potential
		pyprop.PrintMemoryUsage("After Preconditioner Generate Potentials")
	
		#Setup solvers
		tensorPotential.SetupStep(0.0)
		self.SetupSolvers(tensorPotential)

	def SetupSolvers(self, tensorPotential):
		raise NotImplementedException("Please Override")
		
	def Solve(self, psi):
		raise NotImplementedException("Please Override")



class TwoElectronPreconditionerIfpack(TwoElectronPreconditioner):
	"""
	Preconditioner using Ifpack (ILU) to 
	approximately factorize the blocks
	"""

	def __init__(self, psi):
		TwoElectronPreconditioner.__init__(self, psi)

	def ApplyConfigSection(self, conf):
		TwoElectronPreconditioner.ApplyConfigSection(self, conf)
		self.Cutoff = conf.cutoff

	def SetupSolvers(self, tensorPotential):
		pyprop.PrintMemoryUsage("Before Ifpack Setup")

		#Setup the ILU preconditioner
		vector = self.psi.GetData()[:]
		matrix = tensorPotential.PotentialData[:, :]
		solver = IfpackRadialPreconditioner_2()
		basisPairs = tensorPotential.BasisPairs
		solver.Setup(vector, matrix, basisPairs, self.Cutoff)
		self.Solver = solver
		pyprop.PrintMemoryUsage("After Ifpack Setup")

	def Solve(self, psi):
		data = psi.GetData()
		self.Solver.Solve(data[:])
