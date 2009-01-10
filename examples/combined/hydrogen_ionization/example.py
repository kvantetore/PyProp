#Import system modules
import sys
import os
from pylab import *
from numpy import *

#Load pyprop
sys.path.insert(1, os.path.abspath("./pyprop"))
import pyprop
pyprop = reload(pyprop)
pyprop.ProjectNamespace = globals()

#Load the project module
from libpotential import *

def SetupConfig(**args):
	#Decide which config file to use
	configFile = "config.ini"
	if "config" in args:
		configFile = args["config"]

	#Load the config file
	conf = pyprop.Load(configFile)

	#Modify the config
	if "imtime" in args:
		imtime = args["imtime"]
		propSection = conf.Propagation
		dt = abs(propSection.timestep)
		renormalize = False
		if imtime:
			dt = -1.0j * dt
			renormalize = True

		propSection.timestep = dt
		propSection.renormalization = renormalize

	if "amplitude" in args:
		amplitude = args["amplitude"]
		conf.DynamicPotential.amplitude = amplitude

	if "duration" in args:
		duration = args["duration"]
		conf.Propagation.duration = duration
		
	return conf

def SetupProblem(**args):
	conf = SetupConfig(**args)
	prop = pyprop.Problem(conf)
	prop.SetupStep()
	return prop

def FindGroundstate(**args):
	args['imtime'] = True
	prop = SetupProblem(**args)
	
	for t in prop.Advance(10):
		E = prop.GetEnergy()
		print "t = %f, E = %f" % (t, E)

	E = prop.GetEnergy()
	print "Ground State Energy = %f" % E

	SaveWavefunction("groundstate.h5", "/wavefunction", prop.psi)

	return prop

def SaveWavefunction(filename, dataset, psi):
	if pyprop.ProcId == 0:
		if os.path.exists(filename):
			os.unlink(filename)
	pyprop.serialization.SaveWavefunctionHDF(filename, dataset, psi)

def FindIonizationProbability(**args):
	args['imtime'] = False
	args['config'] = "propagation.ini"
	prop = SetupProblem(**args)

	initialPsi = prop.psi.Copy()

	timeList = []
	corrList = []

	for t in prop.Advance(50):
		norm = prop.psi.GetNorm()
		corr = abs(prop.psi.InnerProduct(initialPsi))**2
		timeList.append(t)
		corrList.append(corr)
		print "t = %f, Norm = %f, Corr = %f" % (t, norm, corr)

	norm = prop.psi.GetNorm()
	corr = abs(prop.psi.InnerProduct(initialPsi))**2
	print "Ionization Probability = %f" % norm
	print "Initial state correlation = %f" % corr

	prop.TimeList = array(timeList)
	prop.CorrelationList = array(corrList)

	return prop

def FindIonizationProbabilityAmplitude():
	amplitudeList = r_[0:2:0.2]
	ionizationList = zeros(len(amplitudeList), dtype=double)

	for i in range(len(amplitudeList)):
		prop = FindIonizationProbability(amplitude=amplitudeList[i])
		ionizationList[i] = 1.0 - prop.psi.GetNorm()

	plot(amplitudeList, ionizationList, label="Ionization Probability")
	xlabel("Electric Field Strength")
	ylabel("Ionization Probability")
	legend()

	return amplitudeList, ionizationList


def GetHamiltonMatrixSubspace(prop, l):
	size = prop.psi.GetData().shape[0]
	matrix = zeros((size, size), dtype=complex)
	tempPsi = prop.GetTempPsi()

	for i in range(size):
		prop.psi.GetData()[:] = 0
		prop.psi.GetData()[i,l] = 1

		tempPsi.GetData()[:] = 0
		prop.MultiplyHamiltonian(prop.psi, tempPsi)
		
		matrix[:, i] = tempPsi.GetData()[:,l]
		
	return matrix

def GetEigenstates(l):
	prop = SetupProblem(silent=True)
	M = GetHamiltonMatrixSubspace(prop, l)
	E, V = eig(M)

	#Filter out the eigenstates which are not 0 in origin
	if isinstance(prop.Propagator, pyprop.CombinedPropagator):
		originIndex = prop.Propagator.SubPropagators[0].OriginIndex
		goodIdx = where(abs(V[256, :]) < 1e-10)[0]
		E = real(E[goodIdx])
		V = V[:, goodIdx]

	idx = argsort(E)
	E = E[idx]
	V = V[:,idx]

	return E, V

def SetEigenstate(prop, V, l, idx):
	prop.psi.Clear()
	prop.psi.GetData()[:,l] = V[:,idx]

def PlotFourierRadial(prop):
	prop.Propagator.SubPropagators[0].TransformForward(prop.psi)

	k = prop.psi.GetRepresentation().GetGlobalGrid(0).copy()
	data = sum(abs(prop.psi.GetData())**2, axis=1)

	prop.Propagator.SubPropagators[0].TransformInverse(prop.psi)

	#plot(fftshift(k), fftshift(data))
	plot(k, data)
	


def GetGridLinear(conf, xmax=None, xmin=None):
	if xmin == None: xmin = conf.xmin
	if xmax == None: xmax = conf.xmax
	count = conf.count

	start = xmin
	end = xmax

	if not conf.include_left_boundary:
		count += 1
	if not conf.include_right_boundary:
		count += 1

	dx = (xmax - xmin) / float(count-1)
	if not conf.include_left_boundary:
		start += dx
	if conf.include_right_boundary:
		end += dx
	grid = r_[start:end:dx]

	return array(grid, dtype=double)

class FiniteDifferencePreconditioner:
	def __init__(self, psi):
		self.Rank = psi.GetRank()
		#self.Solver = pyprop.CreateInstanceRank("FiniteDifferenceSolver", self.Rank, globals=globals())
		self.Solver = FiniteDifferenceSolver_2()
		self.psi = psi

	def ApplyConfigSection(self, conf):
		self.PreconditionRank = conf.rank
		potentialNames = conf.potential_evaluation
		self.Potentials = [pyprop.CreatePotential(conf.Config, potName, self.psi) for potName in potentialNames]

	def Setup(self, prop, dt):
		#Add all potentials to solver
		for pot in self.Potentials:
			pot.SetupStep(0)
			if not isinstance(pot, pyprop.StaticPotentialWrapper):
				raise "Only static potentials supported"
			if pot.Storage != pyprop.StaticStorageModel.StorageValue:
				raise "Only StorageValue is supported"

			potentialData = asarray(pot.GetPotential(0), dtype=complex)
			self.Solver.AddPotential(potentialData)
		del self.Potentials

		#Setup solver
		mass = 1.0
		scalingI = 1.0
		scalingH = (1.0j*dt/2.)
		self.Solver.Setup(prop.psi, self.PreconditionRank, 3, 1.0, scalingH)
		print "Setting up!"

	def Solve(self, psi):
		self.Solver.Solve(psi)

