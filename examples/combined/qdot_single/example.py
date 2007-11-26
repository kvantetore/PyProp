#improt system modules
import sys
import os

#Make sure we use the correct pyprop library
sys.path.insert(1, os.path.abspath("../../.."))

#Load and reload pyprop in order to get recent changes
import pyprop
pyprop = reload(pyprop)
from libqdot import *

#numpy an pylab for good measure
import pylab
import numpy
from pylab import *
from numpy import *

def FindGroundstate(**args):
	prop = SetupProblem(**args)
	for t in prop.Advance(10):
		print "t = ", t, ", E = ", prop.GetEnergy()
	hold(False)
	pyprop.Plot1D(prop)
	return prop

def FindEigenvalues(**args):
	prop = SetupProblem(**args)

	solver = pyprop.ArpackSolver(prop)
	solver.Solve()

	print sort(solver.GetEigenvalues().real)

	evCount = len(solver.GetEigenvalues())
	eigenvectors = solver.GetEigenvectors()

	sortIndex = argsort(solver.GetEigenvalues())

	m = 0
	for i in sortIndex:
		shape = prop.psi.GetData().shape
		prop.psi.GetData()[:] = reshape(eigenvectors[i,:], shape)
		prop.psi.Normalize()
		print "E%i = %s" % (m, prop.GetEnergyExpectationValue())
		m += 1

	return solver;

def GetSubEnergy(solver, i):
	E0 = GetSubPropagatorEnergy(solver, 0, i)
	E1 = GetSubPropagatorEnergy(solver, 1, i)
	E2 = GetSubPropagatorEnergy(solver, -1, i)
	
	Esum = E0 + E1 + E2
	Etotal = solver.BaseProblem.GetEnergyExpectationValue()
	print "Radial Energy = %4.4f, \nAngular Energy = %4.4f, \nPotential Energy = %4.4f, \nSum = %4.4f" % (E0, E1, E2, Esum)
	print "Total = %4.4f" % Etotal

def GetSubPropagatorEnergy(solver, subpropagator, i):
	prop = solver.BaseProblem
	t = prop.PropagatedTime
	dt = prop.TimeStep

	#Get temp wavefunction
	tempPsi = prop.GetTempPsi()
	tempPsi.GetData()[:] = 0

	#Activate eigenstate
	GetEigenstate(solver, i)

	#Multiply operator with wavefunction
	if subpropagator >= 0:
		prop.Propagator.SubPropagators[subpropagator].MultiplyHamiltonian(tempPsi, t, dt)
	else:
		prop.Propagator.MultiplyPotential(tempPsi, t, dt)

	#Calculate expectation value by projecting on original state
	E = tempPsi.InnerProduct(prop.psi)
	if abs(E.imag) > 10e-10:
		raise "Complex energy!"
	return E.real

def GetEigenstate(solver, i):
	psi = solver.BaseProblem.psi
	eigenvectors = solver.GetEigenvectors()
	sortIndex = argsort(solver.GetEigenvalues())

	solver.SetEigenvector(psi, sortIndex[i])
	return psi
	
def PlotEigenstate(solver, i):
	psi = GetEigenstate(solver, i)
	pcolor(abs(psi.GetData())**2, shading='flat')

def InnerProduct2(psi1, psi2):
	w = [psi1.GetRepresentation().GetLocalWeights(i) for i in xrange(psi1.GetRank())]
	x = [psi1.GetRepresentation().GetLocalGrid(i) for i in xrange(psi1.GetRank())]

	shape = array([len(i) for i in w])
	weight = zeros(shape, dtype=double)
	pyprop.OuterProduct(w, weight)

	result = conj(psi1.GetData()) * psi2.GetData() *  weight 
	for curLength in result.shape:
		result = sum(result)

	return result

def Normalize2(psi):
	norm = real(InnerProduct2(psi, psi))
	psi.GetData()[:] /= sqrt(norm)

def GetEnergy2(prop):
	tempPsi = prop.GetTempPsi()
	tempPsi.GetData()[:] = 0
	Normalize2(prop.psi)
	prop.MultiplyHamiltonian(tempPsi)

	return InnerProduct2(prop.psi, tempPsi)

def GetBasisMatrix(**args):
	prop = SetupProblem(**args)

	A = pyprop.GetBasisExpansionMatrix(prop)
	E, V = linalg.eig(A)

	print "Eigenvalues = %s" % sort(E.real)
	return A

def GetGroundstateConvergence(**args):
	args['silent'] = True

	n = r_[1:80:5]
	groundstateEnergy = []
	
	for i in xrange(len(n)):
		args['gridSize'] = n[i]
		prop = SetupProblem(**args)
		A = pyprop.GetBasisExpansionMatrix(prop)
		E, V = linalg.eig(A)
		groundstateEnergy += [sort(E.real)[0]]

	return groundstateEnergy
	

def Propagate():
	prop = SetupProblem()
	prop.psi.Normalize()
	x = prop.psi.GetRepresentation().GetLocalGrid(0)
	for t in prop.Advance(50):
		print "t = ", t, ", Norm = ", prop.psi.GetNorm()
		data = abs(prop.psi.GetData())**2
		pylab.plot(x, data)

def SetupProblem(**args):
	conf = SetupConfig(**args)
	prop = pyprop.Problem(conf)
	prop.SetupStep()	
	return prop

def SetupConfig(**args):
	conf = pyprop.Load("config.ini")

	if 'scaling' in args:
		scaling = args['scaling']
		conf.SubRepresentation.scaling = scaling

	if 'silent' in args:
		silent = args['silent']
		conf.Propagation.silent = silent

	if 'gridSize' in args:
		gridSize = args['gridSize']
		conf.SubRepresentation.n = gridSize

	return conf

