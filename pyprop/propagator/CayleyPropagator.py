
class CayleyPropagator(PropagatorBase):
	__Base = PropagatorBase
	
	def __init__(self, psi):
		self.__Base.__init__(self, psi)
		self.Rank = psi.GetRank()
		
		self.Solver = CreateInstanceRank("core.krylov_GmresWrapper", self.Rank)

	def ApplyConfig(self, config):
		self.__Base.ApplyConfig(self, config)

		#Set up the base propagator (which will perform matrix vector multiplications)
		propagatorType = config.Propagation.base_propagator
		self.BasePropagator = propagatorType(self.psi)
		config.Apply(self.BasePropagator)

		#Let solver have config
		config.Apply(self.Solver)

		#Setup preconditioner
		preconditionerName = config.Propagation.preconditioner
		if preconditionerName:
			preconditionerSection = config.GetSection(preconditionerName)
			preconditioner = preconditionerSection.type(self.psi)
			preconditionerSection.Apply(preconditioner)
			self.Preconditioner = preconditioner
		else:
			self.Preconditioner = None


	def ApplyConfigSection(self, configSection): 
		self.__Base.ApplyConfigSection(self, configSection)

		#Set up base propagator
		configSection.Apply(self.BasePropagator)
		#Set up the solver 
		configSection.Apply(self.Solver)
	
	def SetupStep(self, dt):
		self.BasePropagator.SetupStep(dt)
		self.Solver.Setup(self.psi);

		#We need an additional wavefunction to perform matrix-vector muls
		self.TempPsi = self.psi.Copy()

		#Check if the representations is orthogonal
		repr = self.psi.GetRepresentation()
		if all([repr.IsOrthogonalBasis(i) for i in range(self.Rank)]):
			self.IsOrthogonalBasis = True
		else:
			self.IsOrthogonalBasis = False

		#Setup preconditioner
		if self.Preconditioner:
			self.Preconditioner.Setup(self, dt)

	def MultiplyHamiltonian(self, srcPsi, destPsi, t, dt):
		self.BasePropagator.MultiplyHamiltonian(srcPsi, destPsi, t, dt)

	def AdvanceStep(self, t, dt):
		"""
		Advance the solution one timestep by approximating the solution
		by the cayley propagator
		(S + i dt H(t+dt)/2)^-1 psi(t + dt) = (S - i dt H(t)/2) psi(t)
		"""

		#plot(abs(self.psi.GetData()), "r")

		#construct tempPsi = (S + i dt H) psi(t)
		self.MultiplySpH(self.psi, self.TempPsi, -1.0j*dt/2, t, dt)
		
		#plot(abs(self.TempPsi.GetData()), "g--")

		#solve (S - i dt H) psi(t+dt)
		if self.Preconditioner:
			self.Preconditioner.Solve(self.TempPsi)
		callback = lambda srcPsi, dstPsi: self.SolverCallback(srcPsi, dstPsi, t+abs(dt), dt)
		#self.psi.GetData()[:] = 0
		self.Solver.Solve(callback, self.TempPsi, self.psi, False)

		#plot(abs(self.psi.GetData()), "b--")

		#check convergence
		err = self.Solver.GetErrorEstimate()

		#psi2 = self.psi.Copy()
		#self.MultiplySpH(self.psi, psi2, 1.0*dt/2, t+dt, dt)
		#psi2.GetData()[:] -= self.TempPsi.GetData()
		#err2 = psi2.GetNorm()
		#
		if err > 1e-4 :
			print "Error = %s" % (err)
		print self.Solver.GetErrorEstimateList()

	def MultiplySpH(self, sourcePsi, destPsi, hFactor, t, dt):
		"""
		Multiply S + H

		Apply one side of the Cayley Form
		destPsi = (S + hFactor H) sourcePsi

		S is the overlap matrix for non orthogonal basises,
		and the identity operator for orthogonal basises
		"""
		repr = destPsi.GetRepresentation()

		#Multiply
		destPsi.GetData()[:] = sourcePsi.GetData()
		repr.MultiplyOverlap(destPsi)

		destPsi.GetData()[:] *= 1.0/hFactor
		self.BasePropagator.MultiplyHamiltonian(sourcePsi, destPsi, t, dt)
		destPsi.GetData()[:] *= hFactor
		

	def SolverCallback(self, sourcePsi, destPsi, t, dt):
		self.MultiplySpH(sourcePsi, destPsi, 1.0j*dt/2.0, t, dt)

		#Solve left preconditioner in-place
		if self.Preconditioner:
			self.Preconditioner.Solve(destPsi)

		tempPsi = sourcePsi.Copy()
		tempPsi.GetData()[:] -= destPsi.GetData()
		print tempPsi.GetNorm()


		

