

import tables
	
def SaveTensorPotential(filename, pot):
	f = tables.openFile(filename, "w")
	try:
		#Standard parameters
		argList = [pot.PotentialData, 1.0, pot.psi.GetData(), pot.psi.GetData()]
		#Parameters for each geometry info
		for i, geom in enumerate(pot.GeometryList):
			argList += geom.GetMultiplyArguments(pot.psi)

		for i, arg in enumerate(argList):
			if asarray(arg).size != 0:
				f.createArray(f.root, "argument_%i" % i, arg)

	finally:
		f.close()
