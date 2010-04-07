from Configure import conftest
import TaskGen
from TaskGen import feature, extension, before, after
import Task
import Utils
import Build
import Options 

import fabricate as fab

import os.path
import os


curdir = os.path.abspath(os.path.curdir)
fab.setup(dirs=[curdir], )

#------------------------------ Tasks ---------------------------------------

def compile_cpp(self):
	"""
	Task function for compiling a cpp file
	"""
	fab.run(self.env.MPICXX, 
		"-c", self.inputs[0].srcpath(self.env), 
		"-o", self.outputs[0].bldpath(self.env),
		["-I%s" % i for i in getattr(self, "additional_includes", [])],
		["-I%s" %i for i in self.env.INCLUDE],
		["-D%s" %i for i in self.env.DEFINES],
		self.env.CXXFLAGS,
	)

cls = Task.task_type_from_func('cpp', compile_cpp, color='GREEN', ext_in='.cpp', ext_out=".o")
Task.always_run(cls)

def compile_f90(self):
	"""
	Task function for compiling a cpp file
	"""
	fab.run(self.env.MPIF90, 
		"-c", self.inputs[0].srcpath(self.env), 
		"-o", self.outputs[0].bldpath(self.env),
		["-I%s" % i for i in getattr(self, "additional_includes", [])],
		["-D%s" %i for i in self.env.DEFINES],
		self.env.F90FLAGS,
	)

cls = Task.task_type_from_func('f90', compile_f90, color='GREEN', ext_in='.f90', ext_out=".o")
Task.always_run(cls)



def run_link(self):
	"""
	Task function for linking pyprop
	"""

	#Get the local libraries
	localLibs = []
	localLibPaths = []
	for lib in self.local_lib_nodes:
		#get the path where the lib is built
		libPath = lib.parent.abspath(self.env)
		localLibPaths.append(libPath)

		#get the name without the postfixing .so
		baseLibName = os.path.splitext(lib.name)[0]

		#get the name with out the prefixing lib
		if not baseLibName.startswith("lib"):
			raise Exception("Shared library %s does not start with lib!" % lib) 
		localLibs.append(baseLibName[3:])


	fab.run(self.env.MPICXX,
		"-shared",
		"-o", self.outputs[0].bldpath(self.env),
		[i.srcpath(self.env) for i in self.inputs],
		["-L%s" %i for i in localLibPaths],
		["-L%s" %i for i in self.env.LIBPATH],
		["-l%s" %i for i in localLibs],
		["-l%s" %i for i in self.env.LIB],
		["-Wl,-rpath=%s" % i for i in localLibPaths]
	)

cls = Task.task_type_from_func('link', run_link, color='BLUE', ext_in='.o')
Task.always_run(cls)


def pyste_multiple(self):
	"""
	Task function for generating wrappers from several
	pyste files
	"""
	output_folder = self.inputs[0].bld_dir(self.env) + "/pysteoutput"
	if not os.path.exists(output_folder):
		os.makedirs(output_folder)

	includes = [self.path.abspath()]
	includes += self.env.INCLUDE
	includes += self.env.PYSTE_INCLUDE
	includes += getattr(self, "additional_includes", [])

	#run pyste
	fab.run(self.env.PYSTE,
		["-I%s" %i for i in includes],
		["-D%s" %i for i in self.env.DEFINES],
		"--multiple",
		"--out=%s" % output_folder,
		"--module=%s" % self.pyste_module,
		[s.srcpath(self.env) for s in self.inputs],
	)

	#copy output files back into src folder for keeping in git
	for o in self.outputs:
		src = o.abspath(self.env)
		dst = os.path.join(self.path.srcpath(), o.relpath_gen(self.path))
		dst_folder, dst_name = os.path.split(dst)
		if not os.path.exists(dst_folder):
			os.makedirs(dst_folder)
		fab.shell("cp", src, dst)


cls = Task.task_type_from_func('pyste_multiple', pyste_multiple, color='ORANGE', ext_in='.pyste', ext_out=".cpp")
Task.always_run(cls)


def pyste_multiple_generate_main(self):
	"""
	Task function for generating the main wrapper file for
	multiple pyste files
	"""
	output_folder = self.inputs[0].bld_dir(self.env) + "/pysteoutput"
	if not os.path.exists(output_folder):
		os.makedirs(output_folder)

	inputs = [s.srcpath(self.env) for s in self.inputs]
	fab.shell(self.env.PYSTE,
		"--multiple",
		"--out=%s" % output_folder,
		"--module=%s" % self.pyste_module,
		"--generate-main",
		*inputs)

	#copy output files back into src folder for keeping in git
	for o in self.outputs:
		src = o.abspath(self.env)
		dst = os.path.join(self.path.srcpath(), o.relpath_gen(self.path))
		dst_folder, dst_name = os.path.split(dst)
		if not os.path.exists(dst_folder):
			os.makedirs(dst_folder)
		fab.shell("cp", src, dst)



cls = Task.task_type_from_func('pyste_multiple_generate_main', pyste_multiple_generate_main, color='ORANGE', ext_in='.pyste')
Task.always_run(cls)


def pyste_placeholder(self):
	"""
	Placeholder task function for pyste
	"""

	#make sure pysteoutput folder exists
	output_folder = self.inputs[0].bld_dir(self.env) + "/pysteoutput"
	if not os.path.exists(output_folder):
		os.makedirs(output_folder)

	#copy src-stored versions of pyste files to dst
	for o in self.outputs:
		src = os.path.join(self.path.srcpath(), o.relpath_gen(self.path))
		dst = o.abspath(self.env)
		fab.shell("cp", src, dst)

cls = Task.task_type_from_func("pyste_multiple_placeholder", pyste_placeholder, ext_in=".pyste", ext_out=".cpp")
cls = Task.task_type_from_func("pyste_multiple_generate_main_placeholder", pyste_placeholder, ext_in=".pyste", ext_out=".cpp")


def run_python_generate(task):
	"""
	Task for running a python file which generates another
	source file (cpp or fortran)
	"""

	if len(task.inputs) != 1:
		raise Exception("python_generate task requires one python file as input")
	if len(task.outputs) != 1:
		raise Exception("python_generate task requires one file as output")
	
	python = task.env.PYTHON_EXEC
	src = task.inputs[0].srcpath(task.env)
	dst = task.outputs[0].bldpath(task.env)
	params = Utils.to_list(getattr(task, "params", []))
	data = fab.shell(python, src, params)
	f = file(dst, "w")
	f.write(data)
	f.close()


cls = Task.task_type_from_func("python_generate_cpp", run_python_generate, ext_in=".py", ext_out=".cpp")
cls = Task.task_type_from_func("python_generate_fortran", run_python_generate, ext_in=".py", ext_out=".f90")

#------------------------------ Features ---------------------------------------

@feature("compile")
@before("init_compile")
def init_compile_variables(self):
	self.additional_includes = []


@feature("pyprop")
@after("init_compile_variables")
@before("apply_core")
def init_pyprop(self):
	"""
	Feature for pyprop modules depending on core.
	Uses the following variables from the taskgen constructor
		pyprop_path    -   relative path from wscript to pyprop base dir

	"""
	pypropAbsPath = os.path.abspath(os.path.join(self.path.abspath(), self.pyprop_path, "pyprop"))
	self.additional_includes += [pypropAbsPath]


@feature('pyste')
@after("init_compile_variables")
@before('apply_core')
def init_pyste(self):
	"""
	The pyste feature creates tasks to use pyste to create 
	boost::python wrappers

	Furthermore, the generated boost::python wrappers are added
	to the list of files to be compiled. As this function is run before
	apply_core, we can add the generated files to self.allnodes, and they
	will be picked up by apply_core and sent to the corresponding extension
	"""

	pysteFiles = Utils.to_list(self.pyste_files)
	if len(pysteFiles) == 0:
		return

	#Set up nodes for the generated boost::python wrappers 
	modulePath = self.path.abspath(self.env)
	wrapperNodes = []
	for f in Utils.to_list(self.pyste_files):
		pysteNode = self.path.find_resource(f)
		if not pysteNode:
			raise Exception("pyste file %s not found" % f)

		#pyste output path is the path to the current node + pysteoutput
		relativeNodePath = pysteNode.dir(self.env)[len(modulePath)+1:]
		pysteOutput = os.path.join(relativeNodePath, "pysteoutput")
		if not os.path.exists(os.path.join(self.path.abspath(), pysteOutput)):
			print "MAKING %s " % os.path.join(self.path.abspath(), pysteOutput)
			os.makedirs(os.path.join(self.path.abspath(), pysteOutput))

		wrapperFilename = "_" + pysteNode.change_ext(".cpp").file()
		wrapperPath = os.path.join(pysteOutput, wrapperFilename)
		wrapperNode = self.path.find_or_declare(wrapperPath)
		if wrapperNode == None:
			raise Exception("None! %s" % wrapperPath)

		wrapperNodes.append(wrapperNode)

	#main wrapper file
	mainFile = os.path.join(pysteOutput, "_main.cpp") 
	wrapperNodes.append(self.path.find_or_declare([mainFile]))

	#add to all nodes
	self.allnodes += wrapperNodes

	#check whether pyste is enabled, and create the corresponding tasks
	import Options
	pysteForce = getattr(self, "pyste_force", False)
	if Options.options.PysteEnabled or pysteForce:
		tskGenerateWrapper = self.create_task("pyste_multiple")
		tskGenerateWrapper.additional_includes = self.additional_includes
		tskGenerateMain = self.create_task("pyste_multiple_generate_main")
	else:
		#if pyste is not enabled, we create placeholder tasks, 
		#which checks that the generated boost::python wrapper files exist
		tskGenerateWrapper = self.create_task("pyste_multiple_placeholder")
		tskGenerateMain = self.create_task("pyste_multiple_generate_main_placeholder")

	#the wrapper task creates all wrapper files except the main wrapper
	tskGenerateWrapper.set_outputs(wrapperNodes[:-1])
	#the generate main task creates the _main.cpp file
	tskGenerateMain.set_outputs(wrapperNodes[-1])

	#both tasks have all pyste files as input files
	for tsk in [tskGenerateWrapper, tskGenerateMain]:
		tsk.pyste_module = self.pyste_module
		tsk.path = self.path
		for f in Utils.to_list(self.pyste_files):
			n = self.path.find_resource(f)
			tsk.set_inputs(n)



@feature('compile')
@before('apply_core')
def init_compile(self):
	"""
	The compile feature is performed before apply_core, which takes all
	self.allnodes and checks for suitable extension mappings. The extension
	mapping creates compilation tasks, and should add their object output
	files to self.compiled_files, which will be processed by the linker
	after apply_core
	"""
	self.compiled_files = []

@feature('link')
@after('apply_core')
def init_link(self):
	"""
	The link feature takes all compiled files (by readinng self.compiled_files)
	and generates a linking task (static or shared depending on the configuration)
	"""
	tsk = self.create_task("link")
	tsk.set_inputs(self.compiled_files)
	tsk.set_outputs(self.path.find_or_declare("%s.so" % self.pyste_module))
	self.link_task = tsk


@feature('link')
@after("init_link")
def apply_local_libs(self):
	"""
	Applies local library dependencies. This is done after
	init_link to allow all link-tasks to be created, and thus all
	library nodes to be declared. They can then be found by find_resource
	"""
	#get the pyprop path so we know where to look for the libs
	pyprop_path = os.path.join(getattr(self, "pyprop_path", "./"), "pyprop")

	#get task
	tsk = self.link_task
	tsk.local_lib_nodes = []

	#Get the local libs
	libs = getattr(self, "local_libraries", "")
	for lib in Utils.to_list(libs):
		#find the node of the library path
		libpath = os.path.join(pyprop_path, "%s.so" % lib)
		libnode = self.path.find_resource(libpath)
		if libnode == None:
			raise Exception("Could not find library %s (%s)" % (lib, libnode))

		#add dependency
		tsk.deps_nodes.append(libnode)
		tsk.local_lib_nodes.append(libnode)


#------------------------ Extension Mapping----------------------------------


@extension(".cpp")
def extension_cpp(self, node):
	#Create task to compile .cpp -> .o
	tsk = self.create_task('cpp')
	tsk.set_inputs(node)
	outputNode = node.change_ext(".o")
	tsk.set_outputs(outputNode)
	tsk.always = True
	additional_includes = list(self.additional_includes)
	additional_includes += [self.path.srcpath(self.env), self.path.bldpath(self.env)]
	tsk.additional_includes = additional_includes
	
	#Add .o to list of compiled files for linking
	self.compiled_files.append(outputNode)

	#return task, in case someone wants to add options
	return tsk

@extension(".f90")
def extension_f90(self, node):
	#Create task to compile .f90 -> .o
	tsk = self.create_task('f90')
	tsk.set_inputs(node)
	outputNode = node.change_ext(".o")
	tsk.set_outputs(outputNode)
	tsk.always = True
	additional_includes = list(self.additional_includes)
	additional_includes += [self.path.srcpath(self.env), self.path.bldpath(self.env)]
	tsk.additional_includes = additional_includes
	
	#Add .o to list of compiled files for linking
	self.compiled_files.append(outputNode)

	#return task, in case someone wants to add options
	return tsk



#----------------------------- Configure ------------------------------------


@conftest
def detect(conf):
	print "  Detecting PyProp"

	def consolidate_conf(prefix):
		conf.env.INCLUDE += getattr(conf.env, prefix + "_INC")
		conf.env.LIB += getattr(conf.env, prefix + "_LIB")
		conf.env.LIBPATH += getattr(conf.env, prefix + "_LIBPATH")
		conf.env.DEFINES += getattr(conf.env, prefix + "_DEFINES")
		conf.env.CXXFLAGS += getattr(conf.env, prefix + "_CXXFLAGS")
		conf.env.F90FLAGS += getattr(conf.env, prefix + "_F90FLAGS")
		conf.env.LDFLAGS += getattr(conf.env, prefix + "_LDFLAGS")

	check_mpi(conf)
	check_fortran(conf)
	consolidate_conf("FORTRAN")

	check_python(conf)
	consolidate_conf("PYTHON")

	check_boost_python(conf)
	consolidate_conf("BOOST_PYTHON")

	check_blitz(conf)
	consolidate_conf("BLITZ")

	check_trilinos(conf)
	consolidate_conf("TRILINOS")

	check_blas_lapack(conf)
	consolidate_conf("BLAS_LAPACK")

	check_fftw(conf)
	consolidate_conf("FFTW")
	
	check_gsl(conf)
	consolidate_conf("GSL")
	
	check_pyste(conf)

@conftest
def check_fortran(conf):
	print "  - Detecting Fortran90"
	conf.env.FORTRAN_F90FLAGS = ["-Wall",]
	#free form
	conf.env.FORTRAN_F90FLAGS += ["-ffree-form", "-fimplicit-none", "-ffree-line-length-none",]
	#shared module on 64 bit systems
	conf.env.FORTRAN_F90FLAGS += ["-fPIC",]
	#debug
	#conf.env.FORTRAN_F90FLAGS += ["-fbounds-check"]
	#optimization
	#conf.env.FORTRAN_F90FLAGS += ["-O2", "-g"]
	conf.env.FORTRAN_LIB   = ["mpi_f90", "mpi_f77", "gfortran",]

@conftest
def check_blas_lapack(conf):
	print "  - Detecting BLAS/LAPACK"
	conf.env.BLAS_LAPACK_LIB = ["mkl", "guide", "mkl_core", "mkl_lapack", "pthread",]
	conf.env.BLAS_LAPACK_LIBPATH = ["/opt/intel/mkl/10.0.1.014/lib/em64t",]
	conf.env.BLAS_LAPACK_INC = ["/opt/intel/mkl/10.0.1.014/include",]
	conf.env.BLAS_LAPACK_DEFINES = ["PYPROP_USE_BLAS"]

	#TODO add compile tests


@conftest
def check_python(conf):
	print "  - Detecting python"
	import distutils.sysconfig
	import sys
	conf.env.PYTHON_EXEC = sys.executable
	conf.env.PYTHON_LIBPATH = [distutils.sysconfig.get_python_lib()]
	conf.env.PYTHON_INC = [distutils.sysconfig.get_python_inc()]
	conf.env.PYTHON_LIB = []
	conf.env.PYTHON_CXXFLAGS = ["-fPIC"]
	
	#TODO add compile tests


@conftest
def check_boost_python(conf):
	print "  - Detecting boost::python"
	conf.env.BOOST_PYTHON_LIB = ["boost_python-mt"]
	conf.env.BOOST_PYTHON_LIBPATH = ["/opt/boost/lib"]
	conf.env.BOOST_PYTHON_INC = ["/opt/boost/include"]
	
	#TODO add compile tests


@conftest
def check_blitz(conf):
	print "  - Detecting blitz"
	blitz_path = "/home/torebi/prog/pyprop/extern/blitz/build"
	conf.env.BLITZ_LIB = ["blitz"]
	conf.env.BLITZ_LIBPATH = [os.path.join(blitz_path, "lib")]
	conf.env.BLITZ_INC = [os.path.join(blitz_path, "include")]



@conftest
def check_trilinos(conf):
	print "  - Detecting Trilinos"
	trilinos_path = "/opt/trilinos"
	conf.env.TRILINOS_LIB = ["epetra", "teuchos", "ifpack", "anasazi"]
	conf.env.TRILINOS_LIBPATH = [os.path.join(trilinos_path, "lib")]
	conf.env.TRILINOS_INC = [os.path.join(trilinos_path, "include")]
	conf.env.TRILINOS_DEFINES = ["PYPROP_USE_TRILINOS"]
	
	#TODO add compile tests


@conftest
def check_fftw(conf):
	print "  - Detecting fftw"
	conf.env.FFTW_LIB = ["fftw3"]
	conf.env.FFTW_LIBPATH = []
	conf.env.FFTW_INC = []
	conf.env.FFTW_DEFINES = []
	
	#TODO add compile tests



@conftest
def check_mpi(conf):
	print "  - Detecting MPI compilers"
	conf.env.MPICXX = "mpicxx"
	conf.env.MPIF90 = "mpif90"
	conf.env.MPIRUN = "mpirun"
	conf.env.MPIRUN_FLAG_NUMPROCS = "-n %i"
	#TODO add compile tests


@conftest
def check_pyste(conf):
	print "  - Detecting pyste"
	conf.env.GCCXML = conf.find_program(["gccxml"], mandatory=True)
	conf.env.GCCXML_FLAGS = ""
	curpath = os.path.abspath(os.path.curdir)
	conf.env.PYSTE = conf.find_program([curpath + "/extern/pyste/pyste.py"], mandatory=True)
	conf.env.PYSTE_INCLUDE = ["/usr/lib/openmpi/include"]

@conftest
def check_gsl(conf):
	print "  - Detecting gsl"
	conf.env.GSL_LIB = ["gsl"]
	conf.env.GSL_LIBPATH = []
	conf.env.GSL_INC = []
	conf.env.GSL_DEFINES = ["PYPROP_USE_TRILINOS"]


def is_module_enabled(self, moduleName):
	"""
	Checks whether a module is enabled or disabled from the command line
	
	"""
	enabled = True
	
	# If enable-modules is specified, only the explicitly enabled modules 
	# are enabled
	if Options.options.EnableModules != "":
		enabled = moduleName in Options.options.EnableModules.split()
		
	# If the module is mentioned in disable-modules, it is disabled regardless
	if moduleName in Options.options.DisableModules.split():
		enabled = False
		
	return enabled 


import Build
Build.BuildContext.is_module_enabled = is_module_enabled 
