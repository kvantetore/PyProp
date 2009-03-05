class BasisPropagatorInteractionPicture(BasisPropagator):
	"""
	Specialization of BasisPropagator for interaction picture calculations.

	For an hamiltonian which can be written as (Schrodinger picture),

	    H = H0 + H'(t),

	the interaction picture is useful. Applying the unitary transformation,

	    Psi_I = exp(-1j * H0 * t) * Psi_S

	the interaction picture operator becomes

	    F_I = exp(-1j * H0 * t) * H_S * exp(1j * H0 * t)
	"""

	__Base = BasisPropagator
	
	def __init__(self, psi):
		self.__Base.__init__(self, psi)
		self.Rank = psi.GetRank()

		#We need the temp array for solving overlap matrix eqns
		self.TempPsi2 = self.psi.Copy()

	def ApplyConfig(self, config):
		#Create any precalculated potentials 
		self.__Base.ApplyConfig(self, config)
		self.ApplyConfigSection(config.Propagation)

		self.__Base.GeneratePotentials(self, config)
		self.__Base.ConsolidatePotentials(self)
		self.LoadEnergies()
	
	def LoadEnergies(self):
		#Load energies
		h5file = tables.openFile(self.FilenameEnergies)
		self.Energies = array(h5file.getNode(self.DatasetEnergies))
		h5file.close()

	def ApplyConfigSection(self, configSection): 
		self.__Base.ApplyConfigSection(self, configSection)
		self.FilenameEnergies = configSection.filename
		self.DatasetEnergies = configSection.dataset

	def MultiplyHamiltonian(self, srcPsi, destPsi, t, dt):
		#First multiply by exp of left eigenvalues and time
		self.TempPsi2.GetData()[:] = exp(-1j * self.Energies * t) * srcPsi.GetData()

		#Multiply potentials
		destPsi.GetData()[:] = 0
		self.MultiplyPotential(self.TempPsi2, destPsi, t, dt)

		#Then right eigenvalues
		destPsi.GetData()[:] *= exp(1j * self.Energies * t)

		#Solve for all overlap matrices
		repr = srcPsi.GetRepresentation()
		if not all([repr.IsOrthogonalBasis(i) for i in range(self.Rank)]):
			repr.SolveOverlap(destPsi)

	def MultiplyHamiltonianNoOverlap(self, srcPsi, destPsi, t, dt):
		self.MultiplyHamiltonian(srcPsi, destPsi, t, dt)

	def MultiplyHamiltonianBalancedOverlap(self, srcPsi, destPsi, t, dt):
		self.MultiplyHamiltonian(srcPsi, destPsi, t, dt)



