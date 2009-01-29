

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

def GenerateBenchmarkPotential(outputFile=None, potentialList=["ElectronicCouplingPotential"], **args):
	args["useDefaultPotentials"] = False
	args["additionalPotentials"] = potentialList
	args["config"] = "config_benchmark.ini"

	if outputFile == None:
		outputFile = "benchmark/potential_%i.h5" % pyprop.ProcId

	prop = SetupProblem(**args)

	potList = prop.Propagator.BasePropagator.PotentialList
	if len(potList) != 1:
		raise Exception("More than one potential generated, only specify potentials that can be consolidated")

	SaveTensorPotential(outputFile, potList[0])
