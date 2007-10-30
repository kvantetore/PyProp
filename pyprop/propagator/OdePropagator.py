
class OdePropagator(PropagatorBase):
	__Base = PropagatorBase
	
	def __init__(self, psi):
		self.__Base.__init__(self, psi)
		self.Rank = psi.GetRank()
		
		self.OdeWrapper = CreateInstanceRank("core.ODE_OdeWrapper", self.Rank)

	def ApplyConfig(self, config):
		self.__Base.ApplyConfig(self, config)

		#Set up the base propagator (which will perform matrix vector multiplications)
		propagatorType = config.Propagation.base_propagator
		self.BasePropagator = propagatorType(self.psi)
		config.Apply(self.BasePropagator)

		#Set up the ODE propagator
		config.Apply(self.OdeWrapper)

	def ApplyConfigSection(self, configSection): 
		self.__Base.ApplyConfigSection(self, configSection)

		#Set up base propagator
		configSection.Apply(self.BasePropagator)
		#Set up the ODE propagator
		configSection.Apply(self.OdeWrapper)
	
	def SetupStep(self, dt):
		self.BasePropagator.SetupStep(dt)
		self.OdeWrapper.Setup(self.psi);
		#We need an additional wavefunction to perform matrix-vector muls
		self.TempPsi = self.psi.CopyDeep()

		self.PsiDistrib = self.psi.GetRepresentation().GetDistributedModel()
		self.TempDistrib = self.TempPsi.GetRepresentation().GetDistributedModel()

	def MultiplyHamiltonian(self, destPsi, t, dt):
		self.BasePropagator.MultiplyHamiltonian(destPsi, t, dt)

	def AdvanceStep(self, t, dt):
		self.OdeWrapper.AdvanceStep(self.MatVecCallback, self.psi, self.TempPsi, dt, t)

	def MatVecCallback(self, psi, tempPsi, t):
		self.psi.GetRepresentation().SetDistributedModel(self.PsiDistrib)
		tempPsi.GetRepresentation().SetDistributedModel(self.TempDistrib)

		self.BasePropagator.MultiplyHamiltonian(tempPsi, t, 0)

