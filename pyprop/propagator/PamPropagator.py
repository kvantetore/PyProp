
class PamPropagator(PropagatorBase):
	__Base = PropagatorBase
	
	def __init__(self, psi):
		self.__Base.__init__(self, psi)
		self.Rank = psi.GetRank()
		
		self.PampWrapper = CreateInstanceRank("core.krylov_PampWrapper", self.Rank)

	def ApplyConfig(self, config):
		self.__Base.ApplyConfig(self, config)

		#Set up the base propagator (which will perform matrix vector multiplications)
		propagatorType = config.Propagation.base_propagator
		self.BasePropagator = propagatorType(self.psi)
		config.Apply(self.BasePropagator)
		#Set up the expokit propagator
		config.Apply(self.PampWrapper)

	def ApplyConfigSection(self, configSection): 
		self.__Base.ApplyConfigSection(self, configSection)

		#Set up base propagator
		configSection.Apply(self.BasePropagator)
		#Set up the expokit propagator
		configSection.Apply(self.PampWrapper)
	
	def SetupStep(self, dt):
		self.BasePropagator.SetupStep(dt)
		self.PampWrapper.Setup(self.psi);
		#We need an additional wavefunction to perform matrix-vector muls
		self.TempPsi = self.psi.CopyDeep()

		self.PsiDistrib = self.psi.GetRepresentation().GetDistributedModel()
		self.TempDistrib = self.TempPsi.GetRepresentation().GetDistributedModel()

		#Use balanced overlap if BasePropagator has the method MultiplyHamiltonianBalancedOverlap, 
		#or all representations are orthogonal
		self.UseBalancedOverlap = hasattr(self.BasePropagator, "MultiplyHamiltonianBalancedOverlap")
		repr = self.psi.GetRepresentation()
		if all([repr.IsOrthogonalBasis(i) for i in range(self.Rank)]):
			#self.IsOrthogonalBasis = True
			self.IsOrthogonalBasis = False
		else:
			self.IsOrthogonalBasis = False

	def MultiplyHamiltonian(self, destPsi, t, dt):
		self.BasePropagator.MultiplyHamiltonian(destPsi, t, dt)

	def AdvanceStep(self, t, dt):
		if self.UseBalancedOverlap:
			repr = self.psi.GetRepresentation()
			repr.MultiplySqrtOverlap(False, self.psi)
			self.PampWrapper.AdvanceStep(self.MatVecCallback, self.psi, self.TempPsi, dt, t, False)
			repr.SolveSqrtOverlap(False, self.psi)

		elif self.IsOrthogonalBasis:
			repr = self.psi.GetRepresentation()
			repr.MultiplySqrtOverlap(False, self.psi)
			self.PampWrapper.AdvanceStep(self.MatVecCallback, self.psi, self.TempPsi, dt, t, False)
			repr.SolveSqrtOverlap(False, self.psi)
		
		else:
			#this is slow
			self.PampWrapper.AdvanceStep(self.MatVecCallback, self.psi, self.TempPsi, dt, t, True)

	def MatVecCallback(self, psi, tempPsi, dt, t):
		self.psi.GetRepresentation().SetDistributedModel(self.PsiDistrib)
		tempPsi.GetRepresentation().SetDistributedModel(self.TempDistrib)
	
		if self.UseBalancedOverlap:
			#This should be implemented on propagators but its not.
			#TODO: Do not solve for overap matrices on orthogonal ranks that is wrong!
			self.BasePropagator.MultiplyHamiltonianBalancedOverlap(tempPsi, t, dt)

		elif self.IsOrthogonalBasis:
			#This only works on orthogonal representations (Not BSplines)
			repr = psi.GetRepresentation()
			repr.SolveSqrtOverlap(False, psi)
			self.BasePropagator.MultiplyHamiltonian(tempPsi, t, dt)
			repr.MultiplySqrtOverlap(False, psi)
			repr.MultiplySqrtOverlap(False, tempPsi)

		else:
			self.BasePropagator.MultiplyHamiltonian(tempPsi, t, dt)
		

