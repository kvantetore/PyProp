
#------------------------------------------------------------------------------------
#                       Serialization Functions
#------------------------------------------------------------------------------------

def SetupProblemFromFile(file, nodeName=None):
	"""
	Set up problem object and load wavefunction from file.
	"""
	prop = None
	cfgObj = pyprop.serialization.GetConfigFromHDF5(file)
	cfgObj.set("InitialCondition", "type", "None")
	conf = pyprop.Config(cfgObj)
	prop = pyprop.Problem(conf)
	prop.SetupStep()

	GetWavefunctionFromFile(file, prop.psi, nodeName=nodeName)
	
	return prop


def GetWavefunctionFromFile(file, psi, nodeName=None):
	h5file = tables.openFile(file, "r")
	try:
		if nodeName == None:
			for node in h5file.walkNodes():
				if node._v_name == "wavefunction":
					psi.GetData()[:] = node[:]
		else:
			psi.GetData()[:] = h5file.getNode(nodeName)[:]
	finally:
		h5file.close()


def GetArrayFromFile(file, nodeName):
	h5file = tables.openFile(file, "r")
	try:
		for node in h5file.walkNodes():
			if node._v_name == nodeName:
				dataArray = node[:]
	finally:
		h5file.close()
	
	return dataArray


def StoreTensorPotentialMTX(prop, whichPotentials, outFileName, eps = 1e-14):
	fh = open(outFileName, "w")
	for potNum in whichPotentials:
		potential = prop.Propagator.BasePropagator.PotentialList[potNum]
		print "    Processing potential %i: %s" % (potNum, potential.Name)

		basisPairs0 = potential.BasisPairs[0]
		basisPairs1 = potential.BasisPairs[1]
		basisPairs2 = potential.BasisPairs[2]
		
		Count0 = prop.psi.GetData().shape[0]
		Count1 = prop.psi.GetData().shape[1]
		Count2 = prop.psi.GetData().shape[2]

		for i, (x0,x0p) in enumerate(basisPairs0):
			for j, (x1,x1p) in enumerate(basisPairs1):
				for k, (x2,x2p) in enumerate(basisPairs2):
					indexLeft = x2 + (x1 * Count2) + (x0 * Count1 * Count2) 
					indexRight = x2p + (x1p * Count2) + (x0p * Count1 * Count2) 

					#Skip current element if less than eps
					if abs(potential.PotentialData[i, j, k]) < eps:
						continue

					#Write data line to file
					potReal = potential.PotentialData[i, j, k].real
					potImag = potential.PotentialData[i, j, k].imag
					outStr = "%i %i %1.16f %1.16f\n" % (indexLeft, indexRight, potReal, potImag)
					fh.write(outStr)

	fh.close()


def SetupStoredPotentials(potentialNames, **args):
	"""
	Sets up a tensor potential for a list of potentials, 
	and save them to file

	the potentials are saved to the files
	output/potential/<grid_postfix>/<potentialName>.h5
	"""

	#Setup a problem without potentials
	args["useStoredPotentials"] = False
	args["useDefaultPotentials"] = False
	args["useDefaultPreconditionPotentials"] = False
	prop = SetupProblem(**args)
	distr = prop.psi.GetRepresentation().GetDistributedModel()
	conf = prop.Config

	postfix = "_".join(GetRadialGridPostfix(config=conf) + GetAngularGridPostfix(config=conf))
	folder = "output/potentials/%s" % (postfix)
	if pyprop.ProcId == 0:
		if not os.path.exists(folder):
			os.makedirs(folder)

	propagator = prop.Propagator.BasePropagator
	for potName in potentialNames:
		#Generate Potential
		configSection = prop.Config.GetSection(potName)
		potential = propagator.GeneratePotential(configSection)
		#Save to disk
		filename = os.path.join(folder, "%s.h5" % potName)
		pyprop.serialization.SaveTensorPotential(filename, "/potential", potential, distr, prop.Config)
		potential = None


def SetupAllStoredPotentials(**args):

	#Check whether this is for a single L-shell or not
	conf = SetupConfig(**args)
	cfg = conf.AngularRepresentation
	Llist = unique([L for l1, l2, L, M in cfg.index_iterator])
	singleLShell = len(Llist) == 1
	
	potentialNames = [\
		"RadialKineticEnergy1", \
		"RadialKineticEnergy2", \
		"AngularKineticEnergy", \
		"CoulombPotential", \
		"ElectronicCouplingPotential", \
		"ElectronicCouplingPotentialMonopoleTerm", \
		"OverlapPotential", \
		"Absorber", \
		]
	if not singleLShell:
		potentialNames += [\
			"LaserPotentialVelocityDerivativeR1", \
			"LaserPotentialVelocityDerivativeR2", \
			"LaserPotentialVelocity", \
			"DipolePotentialLength", \
			]
	SetupStoredPotentials(potentialNames, **args)
