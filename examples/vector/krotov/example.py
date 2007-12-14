#import system modules
import sys
import os

#pytables
import tables

#Pyprop itself
import pyprop; 

#def GetMatrix():
#	return array([[0, c, d], [c, 0, e], [d,e,0]], dtype=complex)

def GetDiagonalElements():
	"""
	A funtion to provide diagonal (energy) matrix elements
	"""
	data = pylab.load('energies.dat')
	return array(data, dtype=complex)

def Setup():
	"""
	Setup Krotov problem
	"""
	conf = pyprop.Load('config.ini')
	prob = pyprop.Problem(conf)
	prob.SetupStep()
	krotov = pyprop.Krotov(prob)
	return krotov

def Run(krotov):
	krotov = Setup()
	krotov.Run()

	return krotov
