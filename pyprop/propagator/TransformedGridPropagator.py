#---------------------------------------------------------------------
# Very simple propagator for 1D transformed grid problems
#---------------------------------------------------------------------

class TransformedGridPropagator(PropagatorBase):
	__Base = PropagatorBase
	
	def __init__(self, psi):
		self.__Base.__init__(self, psi)
	
		rank = psi.GetRank()
		if rank != 1:
			raise "Only rank==1 is supported by TransformedGridPropagator"
	
		#instantiate transforms		
		self.Propagator = CreateInstanceRank("core.TransformedGridPropagator", 1)
		self.SplittingOrder = 1

	def ApplyConfig(self, config):
		self.__Base.ApplyConfig(self, config)
		
	def ApplyConfigSection(self, configSection): 
		self.__Base.ApplyConfigSection(self, configSection)
		
	def SetupStep(self, dt):
		print "        Setup Potential"
		if self.SplittingOrder == 1:
			self.SetupPotential(dt)
			
		elif self.SplittingOrder == 2:
			self.SetupPotential(dt/2.)

		else:
			raise "invalid splitting order"
	
		print "        Setup Propagator"
		gridParam = self.psi.GetRepresentation().Range.Param
		self.Propagator.Setup(gridParam, dt, self.psi, 0);
	
	def AdvanceStep(self, t, dt):
		if self.SplittingOrder == 1:
			#apply grid potential
			self.ApplyPotential(t, dt)

			#apply propagator
			self.Propagator.AdvanceStep(self.psi)
	
		elif self.SplittingOrder == 2:
			#apply grid potential
			self.ApplyPotential(t, dt/2.)

			#apply propagator
			self.Propagator.AdvanceStep(self.psi)
	
			#apply grid potential
			self.ApplyPotential(t, dt/2.)

		else:
			raise "invalid splitting order"

