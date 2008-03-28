import sys

sys.path.append("./pyprop")
import pyprop

from numpy import array
from numpy import complex
from numpy import zeros

execfile("TensorGenerator.py")

#------------------------------------------------------------------------------------
#                       Test functions
#------------------------------------------------------------------------------------


from libpotential import *

def test():
	conf = pyprop.Load("test.ini")
	gen = TensorPotentialGenerator(config=conf)
	A = gen.GeneratePotential(conf, conf.TestPotential)
	B = gen.GeneratePotential(conf, conf.KineticEnergy)
	Overlap = gen.GeneratePotential(conf, conf.OverlapMatrix)
	
	#A = MultiplyInverseOverlapMatrix(A, Overlap, 0)
	A = MultiplyInverseOverlapMatrix(A, Overlap, 1)
	#B = MultiplyInverseOverlapMatrix(B, Overlap, 0)
	B = MultiplyInverseOverlapMatrix(B, Overlap, 1)

	return A, B, Overlap

def test2():
	conf = pyprop.Load("config.ini")
	prop = pyprop.Problem(conf)
	prop.SetupStep()
	prop.AdvanceStep()
	return prop

import numpy.linalg

def MultiplyInverseOverlapMatrix(A, OverlapMatrix, bsplineRank):
	rank = len(A.shape)

	if rank == 1:
		S = OverlapMatrix
		s = UnpackMatrix(S)
		sInv = numpy.linalg.inv(S)
		a = UnpackMatrix
		b = dot(sInv, a)
		B = b.reshape((A.shape[0]))

	if rank == 2:
		B = A.copy()
		if bsplineRank == 1:
			for otherIndex in range(A.shape[0]):
				Aslice = A[otherIndex, :]
				a = UnpackMatrix(Aslice)
				S = OverlapMatrix[otherIndex, :]
				if norm(a) > 10e-5:
					s = UnpackMatrix(S)
					try:
						sInv = numpy.linalg.inv(s)
					except:
						print otherIndex
						return B
					
					b = dot(sInv, a)
					B[otherIndex, :] = b.reshape(A.shape[1])

		if bsplineRank == 0:
			for otherIndex in range(A.shape[1]):
				Aslice = A[:, otherIndex]
				a = UnpackMatrix(Aslice)
				S = OverlapMatrix[:, otherIndex]
				if norm(a) > 10e-10:
					s = UnpackMatrix(S)
					sInv = numpy.linalg.inv(s)
					
					b = dot(sInv, a)
					B[:,otherIndex] = b.reshape(A.shape[0])

	return B
	
def UnpackMatrix(A):
	rank = len(A.shape)
	if rank == 1:
		N = sqrt(len(A))
		a = A.reshape((N,N))

	if rank == 2:
		N0 = sqrt(A.shape[0])
		N1 = sqrt(A.shape[1])

		a = zeros((N0*N1, N0*N1))
		for i0 in range(N0):
			for j0 in range(N0):
				for i1 in range(N1):
					for j1 in range(N1):
						a[i0*N1 + i1, j0*N1 + j1] = A[i0 + j0*N0, i1 + j1*N1]

	return a


def MatrixMultiply(matrix, inVector, outVector, geometryList):
	rank = len(matrix.shape)

	source = inVector
	dest = outVector
	dest[:] = 0

	if rank == 1:
		pairs0 = geometryList[0].GetBasisPairs()
		for i in xrange(pairs0.shape[0]):
			r0, c0 = pairs0[i, 0], pairs0[i,1]
			dest[r0] += matrix[i] * source[c1]

	if rank == 2:
		pairs0 = geometryList[0].GetBasisPairs()
		pairs1 = geometryList[1].GetBasisPairs()
		
		for i in xrange(pairs0.shape[0]):
			r0, c0 = pairs0[i, 0], pairs0[i,1]
			for j in xrange(pairs1.shape[0]):
				r1, c1 = pairs1[j, 0], pairs1[j,1]
		
				dest[r0, r1] += matrix[i, j] * soure[c0, c1]



