import pylab
from numpy import * 

import sys
import os

#load pyprop
pyprop_path = "../../../"
sys.path.insert(1, os.path.abspath(pyprop_path))
import pyprop
pyprop = reload(pyprop)

#import local modules
from libhydrogen import *

def FindGroundstate(**args):
	#setup propagator
	args['imtime'] = True
	prop = SetupProblem(**args)
	
	#propagate with imaginary time to find groundstate
	for t in prop.Advance(10): 
		E = prop.GetEnergy()
		if pyprop.ProcId == 0:
			print "t = ", t, ", E = ", E

	return prop 

def CreateArpackSolver(**args):
	prop = SetupProblem(**args)
	solver = pyprop.ArpackSolver(prop)
	solver.Solve()
	return prop, solver

def FindPropagatedState(initialPsi, **args):
	#setup propagator
	args['imtime'] = False
	prop = SetupProblem(**args)

	#Load initial condition
	prop.psi.GetData()[:] = initialPsi.GetData()
	
	#propagate with real time to find advanced state
	for t in prop.Advance(10): 
		print "t = ", t

	return prop

def Measure(**args):
	
	#Find ground state
	groundstate = FindGroundstate(**args)

	#Use ground state to propagate in time
	propagation = FindPropagatedState(groundstate.psi, **args)

	#Find correlation on ground state (should be 1 if propagation was exact)
	corr = abs(groundstate.psi.InnerProduct(propagation.psi))**2

	#Find energy of the ground state. For soft == 0.0 this should be
	#E0 = -2.0 a.u.
	energy = groundstate.GetEnergy()

	print "dx = %(dx)s, dt = %(dt)s, soft = %(soft)s" % args
	print "E =", energy, ", corr =", corr

	return  energy, corr

def Plot2DSlice(prop, yIndex):
	xvector = prop.psi.GetRepresentation().GetLocalGrid(0)
	plot(xvector, abs(prop.psi.GetData()[:, yIndex])**2)

def Plot2DSliceFFT(prop, yIndex):
	range = prop.psi.GetRepresentation().GetRange(0)
	kmin = - pi / range.Dx
	kmax = -kmin
	dk = (kmax - kmin)/range.Count
	k = zeros(range.Count, dtype=double)
	k[:range.Count/2+1] = r_[0.0:kmax+dk:dk]
	k[range.Count/2:] = r_[kmin:0.0:dk]
	semilogy(k, abs(fft(prop.psi.GetData()[:, yIndex]))**2)

def run():
	softList = asarray([0.01, 0.02, 0.04, 0.08, 0.16])
	dxList = softList.copy()
	dtList = asarray([0.00115, 0.0025, 0.005, 0.01, 0.02, 0.04, 0.08, 0.16])

	energyArray = zeros((len(softList), len(dxList), len(dtList)), dtype=double)
	corrArray = energyArray.copy()

	for i in r_[0:len(softList)]:
		soft = softList[i]
		for j in r_[0:len(dxList)]:
			dx = dxList[j]
			for k in r_[0:len(dtList)]:
				dt = dtList[k]
				#mesure this point
				energy, corr = Measure(dx=dx, dt=dt, soft=soft)
				#update array
				energyArray[i,j,k] = energy
				corrArray[i,j,k] = corr

	print energyArray
	print corrArray
	pyprop.SavePickleArray("groundstate_energy.dat", energyArray)
	pyprop.SavePickleArray("groundstate_correlation.dat", corrArray)

def load_files():
	dxList = pyprop.LoadPickleArray('dxList.pickle')
	dtList = pyprop.LoadPickleArray('dtList.pickle')
	softList = pyprop.LoadPickleArray('softList.pickle')
	energy = pyprop.LoadPickleArray('groundstate_energy.pickle')
	correlation = pyprop.LoadPickleArray('groundstate_correlation.pickle')

	#correlation and energy is matrix of energy[dx, dt, soft]
	return dxList, dtList, softList, energy, correlation

def plot_corr_dt():
	for i in r_[0:len(softList)]:
		figure()
		title("correlation(dt) for dx=0.01...0.16, soft=" + str(softList[i]))
		for j in r_[0:len(dxList)]:
			plot(dtList, corr[j, :, i], label="dx = " + str(dxList[j]))
		legend()
		

def plot_energy_dt():
	for i in r_[0:len(softList)]:
		figure()
		title("correlation(dt) for dx=0.01...0.16, soft=" + str(softList[i]))
		for j in r_[0:len(dxList)]:
			plot(dtList, energy[j, :, i], label="dx = " + str(dxList[j]))
		legend()


def SetupConfig(**args):
	conf = pyprop.Load("config.ini")

	if 'silent' in args:
		silent = args['silent']
	else:
		silent = pyprop.ProcId != 0
	conf.Propagation.silent=silent


	if 'dx' in args:	
		dx = args['dx']

		#Modify grid
		confSection = conf.Representation
		xmin = confSection.rank0[0]
		xmax = confSection.rank0[1]
		nx = int((xmax - xmin) / dx)

		for rank in r_[0:4]:
			range = [xmin, xmax, nx]
			confSection.Set('rank' + str(rank), range)

	if 'dt' in args:
		dt = args['dt']
		conf.Propagation.timestep = dt

	if 'soft' in args:
		soft = args['soft']
		conf.GridPotential.soft = soft

	if 'imtime' in args:
		imtime = args['imtime']
		if imtime:
			conf.Propagation.timestep = -1.0j * abs(conf.Propagation.timestep)
			conf.Propagation.renormalization = True
		else:
			conf.Propagation.timestep = abs(conf.Propagation.timestep)
			conf.Propagation.renormalization = True
	        	
	return conf	

def SetupProblem(**args):
	conf = SetupConfig(**args)
	prop = pyprop.Problem(conf)
	prop.SetupStep()
	return  prop


