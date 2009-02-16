
#------------------------------------------------------------------------------------
#                       Job Submit Functions
#------------------------------------------------------------------------------------


def Submit(executable=None, writeScript=False, installation="hexagon", **args):
	"""
	Set up job scripts and other necessary stuff to run ionization rate
	cycle scan experiment.
	"""

	#Create jobscript 
	numProcs = args.get("numProcs", 1)

	jscript = submitpbs.SubmitScript()
	jscript.jobname = args.get("jobname", "pyprop")
	jscript.walltime = timedelta(hours=args.get("runHours",1), minutes=0, seconds=0)
	jscript.ppn = args.get("ppn", 4)
	jscript.proc_memory = args.get("proc_memory", "1000mb")
	jscript.nodes = int(ceil(numProcs / float(jscript.ppn)))
	jscript.interconnect = args.get("interconnect", "ib")
	if installation == "stallo":
		jscript.workingdir = args.get("workingDir", "/home/nepstad/proj/argon/")
		jscript.executable = "mpirun -n %s " % (jscript.ppn*jscript.nodes)
		jscript.executable += "python %s" % executable
	elif installation == "hexagon":
		jscript.workingdir = args.get("workingDir")
		jscript.executable = "aprun -n %s " % numProcs
		jscript.executable += "./python-exec %s" % executable
	jscript.parameters = commands.mkarg(repr(args))
	jscript.account = args.get("account", "fysisk")

	#Submit this job
	if writeScript:
		print "\n".join(jscript.CreateScript())
	else:
		jscript.Submit()



def RunSubmit(function, procCount=1, procPerNode=4, *arglist, **argdict):
	"""
	Runs a function on the compute nodes.
	
	function is either a string (name of the function) or a function declaration, 
	which is defined in a namespace available from exec'ing example.py

	function is run by creating a job script for run-function.py, and passing the
	function name and all arguments to the job script.
	"""

	if isinstance(function, str):
		arg1 = function
	else:
		arg1 = function.func_name

	arg2 = commands.mkarg(repr(arglist))
	arg3 = commands.mkarg(repr(argdict))

	jobId = None
	if INSTALLATION == "hexagon":
		submit = submitpbs.SubmitScript()
		submit.procs = procCount
		submit.ppn = min(procPerNode, procCount)
		submit.executable = "./python-exec run-function.py"
		submit.parameters = arg1 + arg2 + arg3
		submit.WriteScript("test.job")
		jobId = submit.Submit()

	elif INSTALLATION == "stallo":
		raise Exception("please to implement")
	
	elif INSTALLATION == "local":
		raise Exception("please to implement")
	
	else:
		raise Exception("Unknown installation '%s'" % INSTALLATION)

	return jobId


