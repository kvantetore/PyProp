#import system modules
import sys
import os
from numpy import conj
import pylab
from numpy import fft

#pytables
import tables

#Pyprop itself
sys.path.append("pyprop")
import pyprop; 
pyprop = reload(pyprop)

def GetMatrixQdot4(psi, conf):
	"""
	Get couplings between the four states considered for
	the optimal control PRL (Saelen et al 2008).
	"""
	omega = conf.omega
	V12 = -1.0067/sqrt(omega)
	V13 = 0.0
	V14 = -0.0696/sqrt(omega)
	V23 = 1.0494/sqrt(omega)
	V24 = 0.0
	V34 = -3.0782/sqrt(omega)
	potential = array([ [0, V12, V13, V14], \
	                    [V12, 0, V23, V24], \
	                    [V13, V23, 0, V34], \
	                    [V14, V24, V34, 0]],\
	                    dtype=complex)
	
	return potential

def GetSparseMatrix(psi, config):
	matrix = pylab.load("d130_50stk-matel")
	row = array(matrix[:,0], dtype=int) - 1
	col = array(matrix[:,1], dtype=int) - 1
	matelem = array(matrix[:,2], dtype=complex)

	return row, col, matelem

def GetDenseMatrix(psi, config):
	row, col, matelem = GetSparseMatrix(psi, config)

	N = max(col)+1
	matrix = zeros((N, N), dtype=complex)
	for i in range(len(row)):
		r, c, v = row[i], col[i], matelem[i]
		matrix[r, c] = v
		matrix[c, r] = conj(v)

	return matrix

def GetDiagonalElementsQdot4(psi, config, potential):
	"""
	"""
	potential[:] = [0.241995, 0.324479, 0.380237, 0.382176]

def GetDiagonalElements(psi, config, potential):
	"""
	A funtion to provide diagonal (energy) matrix elements
	"""

	potential[:] = pylab.load(config.file_name)[:config.size] * config.scaling


def GetFieldMatrixElements(psi, config):
	"""
	"""

def Setup(**args):
	"""
	Setup Krotov problem
	"""
	conf = pyprop.Load('config.ini')

	if "timestep" in args:
		config.Propagation.timestep = args["timestep"]

	prob = pyprop.Problem(conf)
	prob.SetupStep()
	krotov = pyprop.Krotov(prob)
	return krotov

def Run():
	krotov = Setup()
	krotov.Run()

	return krotov

def GetPulseSpectrum(krotov):
	"""
	"""
	
	dw = 2 * pi * fft.fftshift(fft.fftfreq(len(krotov.ControlVector), krotov.TimeGridResolution))
	pulseSpectrum = fft.fftshift(fft.fft(krotov.ControlVector))

	return dw, pulseSpectrum


def Commutator(A,B):
	"""
	Computes the commutator [A,B] = AB-BA
	"""
	return dot(A,B) - dot(B, A)


def NestedCommutatorTypeI(A, B, depth):
	"""
	Computes the nested commutator [A,[A,...,[A,B]...]],
	that is, an inner commutator of C = [A,B] and then
	then recursively the commutator [A,C] until desired 
	depth is reached.
	"""
	C = Commutator(A,B)
	for d in range(depth):
		C = Commutator(A,C)
	
	return C
	

def TextToHDFDense(fileName, vectorSize, scaling):

	groupName = 'doubledot'
	datasetPath = '/' + groupName + '/matrixelements_50'
	data = pylab.load(fileName)

	fileh5 = tables.openFile(fileName + '.h5', 'w')
	
	try:
		group = fileh5.createGroup(fileh5.root, groupName)
		h5array = pyprop.serialization.CreateDataset(fileh5, datasetPath, (vectorSize,vectorSize))

		#Fill h5array with matrix element data
		for i in range(shape(data)[0]):
			row = int(data[i,0]) - 1
			col = int(data[i,1]) - 1
			matel = data[i,2] * scaling
			h5array[row, col] = matel 
			h5array[col,row] = matel

	finally:
		fileh5.close()
