from pylab import *
from numpy import * 

import sys
import os

#load pyprop
pyprop_path = "../../../"
sys.path.insert(1, os.path.abspath(pyprop_path))
import pyprop
pyprop = reload(pyprop)

from libhydrogen import *

conf = pyprop.Load("config.ini")

def ModifyConf(conf, dx=None, dt=None, soft=None):

	if dx != None:
		#Modify grid
		confSection = conf.Representation
		xmin = confSection.xmin
		xmax = confSection.xmax
		nx = int((xmax - xmin) / dx)
		#pyprop does not take dx as parameter, but rather
		#the number of grid points, from which the effective
		#dx is computed.
		effectiveDx = (xmax - xmin) / nx

		#in order to maximize the distance between the two points 
		#clostest to the origin, we translate the grid by delta.
		#x0 is the point closest to the origin on the negative side.
		x0 = xmin - effectiveDx * floor(xmin/effectiveDx)
		delta = x0 + effectiveDx / 2.0

		#x = r_[xmin+delta : xmax+delta : effectiveDx]
		#print "min(x) =", min(abs(x)/effectiveDx), "dx"
		
		for rank in r_[0:4]:
			range = [xmin + delta, xmax + delta, nx]
			confSection.Set('rank' + str(rank), range)

	if dt != None:
		conf.Propagation.timestep = dt

	if soft != None:
		conf.GridPotential.soft = soft
	
	return conf	

def FindGroundState(conf):
	#setup propagator
	conf.Propagation.timestep = -1.0j * abs(conf.Propagation.timestep)
	conf.Propagation.renormalization = True
	prop = pyprop.Problem(conf)
	prop.SetupStep()
	
	#propagate with imaginary time to find groundstate
	for t in prop.Advance(10): pass

	return prop 

def FindPropagatedState(conf, initialCondition):
	#setup propagator
	conf.Propagation.timestep = abs(conf.Propagation.timestep)
	conf.Propagation.renormalization = False
	prop = pyprop.Problem(conf)
	prop.SetupStep()

	#Load initial condition
	prop.psi.GetData()[:] = initialCondition.GetData()
	
	#propagate with real time to find advanced state
	for t in prop.Advance(10): pass

	return prop

def Measure(dx, dt, soft):
	conf = pyprop.Load("config.ini")
	conf = ModifyConf(conf, dx, dt, soft)

	#Find ground state
	groundstate = FindGroundState(conf)

	#Use ground state to propagate in time
	propagation = FindPropagatedState(conf, groundstate.psi)

	#Find correlation on ground state (should be 1 if propagation
	#was exact)
	corr = abs(groundstate.psi.InnerProduct(propagation.psi))**2

	#Find energy of the ground state. For soft == 0.0 this should be
	#E0 = -2.0 a.u.
	energy = groundstate.GetEnergy()

	print "dx =", dx, ", dt =", dt, ", soft =", soft
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
				energy, corr = Measure(dx, dt, soft)
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
