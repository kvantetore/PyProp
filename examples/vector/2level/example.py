#import system modules
import sys
import os

#pytables
import tables

#numerical modules
from numpy import *
from numpy import max as nmax
import pylab

#home grown modules
pyprop_path = "../../../"
sys.path.insert(1, os.path.abspath(pyprop_path))
import pyprop
pyprop = reload(pyprop)

c = 0.1
d = 0.005
e = 0.15

def GetSparseMatrix(psi, config):
	#matrix = pylab.load("/home/raymond/sci/dev/Krotov/matrixElements/d130_50stk-matel")
	matrix = pylab.load("/home/raymond/sci/dev/qdot4d/CreateMatrixElements/lene/d130_50stk_newRay-matel")
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

def GetDiagonalElements(psi, config, potential):
	potential[:] = pylab.load(config.file_name) * config.scaling

def Propagate():
	conf = pyprop.Load("config.ini")
	prop = pyprop.Problem(conf)
	prop.SetupStep()

	init = prop.psi.Copy()
	corr = []
	times = []
	corr.append(abs(prop.psi.GetData()[:])**2)
	times += [0]

	for t in prop.Advance(20):
		#corr += [abs(prop.psi.InnerProduct(init))**2]
		corr.append(abs(prop.psi.GetData()[:])**2)
		times += [t]
		print "Time = ", t, ", initial state correlation = ", corr[-1][0]

	corr.append(abs(prop.psi.GetData()[:])**2)
	times += [prop.PropagatedTime]
	print "Time = ", prop.PropagatedTime, ", initial state correlation = ", corr[-1][0]
	prop.corr = asarray(corr)

	return prop


def CompareFortran(**args):
	conf = pyprop.Load("config_compare_fortran.ini")
	prop = pyprop.Problem(conf)
	prop.SetupStep()

	init = prop.psi.Copy()

	for t in prop.Advance(5):
		corr = abs(prop.psi.InnerProduct(init))**2
		print "Time = %f, initial state correlation = %f" % (t, corr)

	corr = abs(prop.psi.InnerProduct(init))**2
	t = prop.PropagatedTime
	print "Time = %f, initial state correlation = %f" % (t, corr)

	#Load fortran data and compare
	fdata = pylab.load("fortran_propagation.dat")
	print "Max difference pyprop/fortran: %e" % nmax(abs(prop.psi.GetData())**2 - fdata[1:])

	return prop



#class SparseMatrix(tables.IsDescription):
#	RowIndex = tables.Int32Col(pos=1)
#	ColIndex = tables.Int32Col(pos=2)
#	MatrixElement = tables.ComplexCol(itemsize=16, pos=3)

#def CreateHDF5FileFromText(fileName,vectorSize):
#	RowIndex = tables.Int32Col(pos=1)
#	ColIndex = tables.Int32Col(pos=2)
#	MatrixElement = tables.ComplexCol(itemsize=16, pos=3)
#
#	data = pylab.load(fileName)
#
#	fileh5 = tables.openFile(fileName + '.h5', 'w')
#	# Create a new group
#	group = fileh.createGroup(fileh.root, "doubledot")
#
#	# Create a new table in group groupName
#	tableName = "matrixElements"
#	tableTitle = "Double quantum dot matrix elements"
#	table = fileh.createTable(group, tableName, SparseMatrix, tableTitle, Filters(1))
#	sparseMatrixRow = table.row
#
## Fill the table with 10 particles
#for i in xrange(10):
#    # First, assign the values to the Particle record
#    particle['name']  = 'Particle: %6d' % (i)
#    particle['lati'] = i
#    particle['longi'] = 10 - i
#    particle['pressure'] = float(i*i)
#    particle['temperature'] = float(i**2)
#    # This injects the row values.
#    particle.append()


def TextToHDFDense(fileName, vectorSize):

	groupName = 'doubledot'
	datasetPath = '/' + groupName + '/matrixelements'
	data = pylab.load(fileName)

	fileh5 = tables.openFile(fileName + '.h5', 'w')
	
	try:
		group = fileh5.createGroup(fileh5.root, groupName)
		h5array = pyprop.serialization.CreateDataset(fileh5, datasetPath, (vectorSize,vectorSize))

		#Fill h5array with matrix element data
		for i in range(shape(data)[0]):
			row = int(data[i,0]) - 1
			col = int(data[i,1]) - 1
			matel = data[i,2]
			h5array[row, col] = matel 
			h5array[col,row] = matel

	finally:
		fileh5.close()



