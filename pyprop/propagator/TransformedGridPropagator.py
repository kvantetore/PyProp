#---------------------------------------------------------------------
# Very simple propagator for 1D transformed grid problems
#---------------------------------------------------------------------
import numpy.linalg

class TransformedGridPropagator(PropagatorBase):
	__Base = PropagatorBase
	
	def __init__(self, psi):
		self.__Base.__init__(self, psi)
	
		self.IncludePotential = False

		rank = psi.GetRank()
		if rank != 1:
			raise "Only rank==1 is supported by TransformedGridPropagator"
	
		self.SplittingOrder = 2

	def ApplyConfig(self, config):
		self.__Base.ApplyConfig(self, config)
		self.__config = config
		
	def ApplyConfigSection(self, configSection): 
		self.__Base.ApplyConfigSection(self, configSection)

		if hasattr(configSection, "include_potential"):
			if configSection.include_potential:
				print "Warning: Including potential into TransformedGridPropagator."
				print "This is only implented as a test, and should not be used for production runs. "
				print "Interface may change"
				self.IncludePotential = True

		#instantiate transforms		
		propagatorClass = "core.TransformedGridPropagator"
		if hasattr(configSection, "propagator_class"):
			propagatorClass = configSection.propagator_class
		self.Propagator = CreateInstanceRank(propagatorClass, 1)
		configSection.Apply(self.Propagator)
		
	def SetupStep(self, dt):
		print "        Setup Potential"
		if self.SplittingOrder == 1:
			effectiveDt = dt	
			
		elif self.SplittingOrder == 2:
			effectiveDt = dt/2.

		else:
			raise "invalid splitting order"
		self.SetupPotential(effectiveDt)
	
		print "        Setup Propagator"
		gridParam = self.psi.GetRepresentation().Range.Param
		self.Propagator.Setup(gridParam, dt, self.psi, 0);

		if self.IncludePotential:
			#Include the first potential in the propagator
			pot = self.PotentialList[0].GetPotential(effectiveDt)
			#Make sure we don't apply it later
			del self.PotentialList[0]

			#Set up propagation matrix
			diffMatrix = self.Propagator.GetDifferentiationMatrix()
			diffMatrix += diag(pot)
			L, X = numpy.linalg.eig(diffMatrix)
			L = exp(- 1.0j * dt * L)
			L = numpy.diag(L)
			Xinv = numpy.linalg.inv(X)

			dot = numpy.dot
			propMatrix = dot(dot(X, L), Xinv)
			self.Propagator.GetPropagationMatrix()[:] = propMatrix

	
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

	def MultiplyHamiltonian(self, dstPsi, t, dt):
		if self.SplittingOrder == 2:
			dt /= 2.0
		self.Propagator.ApplyDifferentiationMatrix(self.psi, dstPsi)
		self.MultiplyPotential(dstPsi, t, dt)
