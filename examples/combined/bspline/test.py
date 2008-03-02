import pyprop
pyprop = reload(pyprop)

from libpotential import *

def Setup():
	conf = pyprop.Load('config.ini')
	bsplinerepr = pyprop.core.BSplineRepresentation()
	bsplinerepr.ApplyConfigSection(conf.BSplineRepresentation)

	bsplinegridrepr = pyprop.core.BSplineGridRepresentation()
	bsplinegridrepr.SetupRepresentation(bsplinerepr.GetBSplineObject())

	prop = pyprop.Problem(conf)
	prop.SetupStep()

	return prop

def Run():
	prop = Setup();
	for t in prop.Advance(10): print "t = %f, E = %f" % (t, prop.GetEnergy())

	return prop

def GetHamiltonMatrix(prop):
	size = prop.psi.GetData().size
	matrix = zeros((size, size), dtype=complex)
	tempPsi = prop.GetTempPsi()

	for i in range(size):
		prop.psi.GetData()[:] = 0
		prop.psi.GetData()[i] = 1

		tempPsi.GetData()[:] = 0
		prop.MultiplyHamiltonian(tempPsi)
		
		matrix[:, i] = tempPsi.GetData()
		
	return matrix
	
