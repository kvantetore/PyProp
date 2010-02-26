
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

	def RestartPropagation(self, timestep, startTime, propagationTime):
		"""
		Set start time on OdeWrapper.
		"""
		self.OdeWrapper.SetStartTime(startTime)

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
		
		self.MatVecCount = 0

	def MultiplyHamiltonian(self, srcPsi, destPsi, t, dt):
		self.BasePropagator.MultiplyHamiltonian(srcPsi, destPsi, t, dt)

	def AdvanceStep(self, t, dt):
		self.OdeWrapper.AdvanceStep(self.MatVecCallback, self.psi, self.TempPsi, dt, t)

	def MatVecCallback(self, psi, tempPsi, t):
#		self.psi.GetRepresentation().SetDistributedModel(self.PsiDistrib)
#		tempPsi.GetRepresentation().SetDistributedModel(self.TempDistrib)
	#	if self.MatVecCount == -1:
	#		print "In Psi = ", psi.GetData()
	#		print "inpsi[0] = %.17g " % abs(psi.GetData()[0])**2
	#		print "inpsi[1] = %.17g " % abs(psi.GetData()[1])**2
	#		print "inpsi[1] = %.17g " % abs(psi.GetData()[2])**2

		self.BasePropagator.MultiplyHamiltonian(psi, tempPsi, t, 0)

	#	if self.MatVecCount == -1:
	#		print "t = %.17g, %.17g" % (t, self.OdeWrapper.GetPropagatedTime())
	#		print "Out Psi = ", tempPsi.GetData()
	#		print "psi[0] = %.17g " % tempPsi.GetData()[0]
	#		print "psi[1] = %.17g " % tempPsi.GetData()[1]
	#		print "psi[1] = %.17g " % tempPsi.GetData()[2]
	#		raise Exception("STOP!")
		self.MatVecCount += 1

	def GetBasisFunction(self, rank, basisIndex):
		return self.BasePropagator.GetBasisFunction(rank, basisIndex)

	def GetPotentialList(self):
		return self.BasePropagator.GetPotentialList()

