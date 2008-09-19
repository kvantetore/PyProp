

def CommonSetupConfig(**args):
	"""
	Loads a config file and sets a some common
	configuration variables based on args

	conf.SetValue is used to enable serialization
	of the configuration file
	"""

	configFile = args.get("configFile", "config.ini")
	conf = pyprop.Load(configFile)

	if "silent" in args:
		silent = args["silent"]
		conf.Propagation.silent = silent

	if "initialResidual" in args:
		conf.SetValue("Arpack", "krylov_use_random_start", False)

	if "dt" in args:
		dt = args["dt"]
		conf.SetValue("Propagation", "timestep", dt)

	if "imtime" in args:
		imtime = args["imtime"]
		if imtime:
			conf.SetValue("Propagation", "timestep", -1.0j * abs(conf.Propagation.timestep))
			conf.SetValue("Propagation", "renormalization", True)
		else:
			conf.SetValue("Propagation", "timestep", abs(conf.Propagation.timestep))
			conf.SetValue("Propagation", "renormalization", False)

	if "potentialList" in args:
		potentialList = args["potentialList"]
		newPotentialList = conf.Propagation.potential_evaluation + potentialList
		conf.SetValue("Propagation", "potential_evaluation", newPotentialList)

	if "nuclearSeparation" in args:
		nuclearSeparation = args["nuclearSeparation"]
		conf.SetValue("Potential", "nuclear_separation", nuclearSeparation)
	
	if "duration" in args:
		duration = args["duration"]
		conf.SetValue("Propagation", "duration", duration)

	if "innerGridCount" in args:
		conf.ElectronRepresentation.inner_count = args["innerGridCount"]

	return conf


def SetupProblem(**args):
	"""
	Creates a config object from args, uses it to
	create and set up a Problem object
	"""
	conf = SetupConfig(**args)
	prop = pyprop.Problem(conf)
	prop.SetupStep()
	return prop



