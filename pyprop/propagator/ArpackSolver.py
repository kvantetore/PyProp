
class ArpackSolver:

	def __init__(self, prop):
		self.BaseProblem = prop
		self.Rank = prop.psi.GetRank()

		#Create a copy of the wavefunction to calculate H|psi>
		self.TempPsi = prop.psi.CopyDeep()

		self.Solver = CreateInstanceRank("core.krylov_ArpackPropagator", self.Rank)
		self.Solver.Setup(prop.psi)

	def Solve(self):
		psi = self.BaseProblem.psi;
		tempPsi = self.TempPsi

		self.Solver.Solve(self.__MatVecCallback, psi, tempPsi)

	def __MatVecCallback(self, psi, tempPsi):
		self.BaseProblem.Propagator.MultiplyHamiltonian(tempPsi, 0, 0)

		#outN = self.TempPsi.GetNorm()
		#if outN < 0.001:
		#	pylab.plot(abs(psi.GetData())**2)
		#	pylab.plot(abs(tempPsi.GetData())**2)
		#	print "inN = ", inN
		#	print "outN = ", outN
		#	assert False
	


