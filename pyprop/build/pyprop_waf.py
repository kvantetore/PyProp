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

def set_options(opt):
	options_mpi(opt)
	options_fortran(opt)
	options_python(opt)
	options_boost_python(opt)
	options_blitz(opt)
	options_trilinos(opt)
	options_blas_lapack(opt)
	options_fftw(opt)
	options_gsl(opt)
	options_pyste(opt)

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


def add_option_flags_fortran(opt, optName, variableName, **args):
	opt.add_option(optName, action="store", dest=variableName, default=None, **args)
	opt.add_option(optName + "-disable-default", action="store_true", dest=variableName + "DisableDefault", default=False)

def set_option_flags(conf, variableName, flagName):
	optVar = getattr(Options.options, variableName)
	disableDefault = getattr(Options.options, variableName + "DisableDefault")

	if not disableDefault:
		flagValue = getattr(conf.env, flagName)
	else:
		flagValue = []

	if not optVar is None:
		flagValue += optVar.split(" ")

	setattr(conf.env, flagName, flagValue)

def add_options_lib(opt, libname, defaultLib, defaultPath="", defaultDefines=""):
	l1 = libname.lower().replace("_","-")
	l2 = libname.upper()

	opt.add_option("--%s-lib" % l1, action="store", dest="%s_LIB" % l2, default=defaultLib, help="List of Libraries from %s to link against, default=%s" % (l1, defaultLib))
	opt.add_option("--%s-path" % l1, action="store", dest="%s_PATH" % l2, default=defaultPath, help="Path to %s, default=%s" % (l1, defaultPath))
	opt.add_option("--%s-libpath" % l1, action="store", dest="%s_LIBPATH" % l2, default="lib", help="Path to the libraries of %s, default=lib" % l1)
	opt.add_option("--%s-inc" % l1, action="store", dest="%s_INC" % l2, default="include", help="Path to the include files of %s, default=include" % l1)
	opt.add_option("--%s-defines" % l1, action="store", dest="%s_DEFINES" % l2, default=defaultDefines, help="Extra defines required by %s, default=%s" % (l1, defaultDefines))

def split_names(namestring):
	return namestring.replace(",", " ").replace(";", " ").split()


def set_options_lib(conf, libname):
	path = getattr(Options.options, "%s_PATH" % libname.upper())

	libpath = getattr(Options.options, "%s_LIBPATH" % libname.upper())
	if not os.path.isabs(libpath):
		if path == "":
			libpath = []
		else:
			libpath = [os.path.join(path, libpath)]
	else:
		libpath = [libpath]
	setattr(conf.env, "%s_LIBPATH" % libname.upper(), libpath)

	inc = getattr(Options.options, "%s_INC" % libname.upper())
	if not os.path.isabs(inc):
		if path == "":
			inc = []
		else:
			inc = [os.path.join(path, inc)]
	else:
		inc = [inc]
	setattr(conf.env, "%s_INC" % libname.upper(), inc)

	lib = getattr(Options.options, "%s_LIB" % libname.upper())
	setattr(conf.env, "%s_LIB" % libname.upper(), split_names(lib))

	defines = getattr(Options.options, "%s_DEFINES" % libname.upper())
	setattr(conf.env, "%s_DEFINES" % libname.upper(), split_names(defines))


def options_fortran(opt):
	opt = opt.parser.add_option_group("Fortran options")
	opt.add_option("--fortran-compiler", action="store", dest="FortranCompiler", default="mpif90", help="Fortran compiler to use")
	add_option_flags_fortran(opt, "--fortran-flags", "FortranFlags", help="Flags to pass to the fortran compiler")

@conftest
def check_fortran(conf):
	print "  - Detecting Fortran90"

	conf.env.FORTRAN_F90FLAGS = ["-Wall",]
	#free form
	conf.env.FORTRAN_F90FLAGS += ["-ffree-form", "-fimplicit-none", "-ffree-line-length-none",]
	#shared module on 64 bit systems
	conf.env.FORTRAN_F90FLAGS += ["-fPIC",]
	set_option_flags(conf, "FortranFlags", "FORTRAN_F90FLAGS")
	#debug
	#conf.env.FORTRAN_F90FLAGS += ["-fbounds-check"]
	#optimization
	#conf.env.FORTRAN_F90FLAGS += ["-O2", "-g"]
	conf.env.FORTRAN_LIB   = ["mpi_f90", "mpi_f77", "gfortran",]

def options_blas_lapack(opt):
	opt = opt.parser.add_option_group("BLAS/LAPACK options")
	add_options_lib(opt, "blas_lapack", "blas lapack")
	opt.add_option("--enable-blas-mkl", action="store_true",
			dest="PypropUseMKL", default=False,
			help="Enable Intel MKL for BLAS/LAPACK.")

@conftest
def check_blas_lapack(conf):
	if Options.options.PypropUseMKL:
		print "  - Detecting Intel MKL BLAS/LAPACK"
		Options.options.BLAS_LAPACK_LIB = "mkl, guide, mkl_core, mkl_lapack, pthread"
		set_options_lib(conf, "blas_lapack")
		conf.env.BLAS_LAPACK_DEFINES += ["PYPROP_USE_BLAS_MKL"]
	else:
		print "  - Detecting BLAS/LAPACK"
		set_options_lib(conf, "blas_lapack")

	#TODO add compile tests


def options_python(opt):
	opt = opt.parser.add_option_group("python options")
	import distutils.sysconfig
	import sys
	opt.add_option("--python-exec", action="store", dest="PYTHON_EXEC", default=sys.executable)
	opt.add_option("--python-libpath", action="store", dest="PYTHON_LIBPATH", default=distutils.sysconfig.get_python_lib())
	opt.add_option("--python-inc", action="store", dest="PYTHON_INC", default=distutils.sysconfig.get_python_inc())
	opt.add_option("--python-lib", action="store", dest="PYTHON_LIB", default="")
	opt.add_option("--python-cxxflags", action="store", dest="PYTHON_CXXFLAGS", default="-fPIC")


@conftest
def check_python(conf):
	print "  - Detecting python"
	conf.env.PYTHON_EXEC = Options.options.PYTHON_EXEC
	conf.env.PYTHON_LIBPATH = split_names(Options.options.PYTHON_LIBPATH)
	conf.env.PYTHON_INC = split_names(Options.options.PYTHON_INC)
	conf.env.PYTHON_LIB =  split_names(Options.options.PYTHON_LIB)
	conf.env.PYTHON_CXXFLAGS = split_names(Options.options.PYTHON_CXXFLAGS)

	#TODO add compile tests

def options_boost_python(opt):
	opt = opt.parser.add_option_group("boost::python options")
	add_options_lib(opt, "boost_python", "boost_python-mt", "/opt/boost")

@conftest
def check_boost_python(conf):
	print "  - Detecting boost::python"
	set_options_lib(conf, "boost_python")

	#TODO add compile tests


def options_blitz(opt):
	opt = opt.parser.add_option_group("blitz options")
	add_options_lib(opt, "blitz", "blitz", os.path.abspath("./extern/blitz/build"))

@conftest
def check_blitz(conf):
	print "  - Detecting blitz "
	set_options_lib(conf, "blitz")


def options_trilinos(opt):
	opt = opt.parser.add_option_group("trilinos options")
	add_options_lib(opt, "trilinos", "epetra, teuchos, ifpack, anasazi", "/opt/trilinos")
	opt.add_option("--disable-trilinos", action="store_false", dest="PypropUseTrilinos", default=True, help="Disables integration with trilinos completely.")
	opt.add_option("--disable-trilinos-tpetra", action="store_false", dest="PypropUseTrilinosTpetra", default=True, help="Disables trilinos/tpetra integration. Only applies if trilinos is enabled")

@conftest
def check_trilinos(conf):
	if Options.options.PypropUseTrilinos:
		print "  - Detecting Trilinos"
		set_options_lib(conf, "trilinos")
		conf.env.TRILINOS_DEFINES += ["PYPROP_USE_TRILINOS"]

		if Options.options.PypropUseTrilinosTpetra:
			conf.env.TRILINOS_DEFINES += ["PYPROP_USE_TRILINOS_TPETRA"]

	#TODO add compile tests


def options_fftw(opt):
	opt = opt.parser.add_option_group("fftw options")
	add_options_lib(opt, "fftw", "fftw3")

@conftest
def check_fftw(conf):
	print "  - Detecting fftw"
	set_options_lib(conf, "fftw")

	#TODO add compile tests


def options_mpi(opt):
	opt = opt.parser.add_option_group("mpi options")
	opt.add_option("--fortran-compiler", action="store", dest="MPIF90", default="mpif90", help="Fortran compiler to use")
	opt.add_option("--cxx-compiler", action="store", dest="MPICXX", default="mpicxx", help="C++ compiler to use")
	opt.add_option("--mpirun", action="store", dest="MPIRUN", default="mpirun", help="Command to start parallel jobs")
	opt.add_option("--mpirun-flag-numprocs", action="store", dest="MPIRUN_FLAG_NUMPROCS", default="-n %i", help="Parameters to mpirun")


@conftest
def check_mpi(conf):
	print "  - Detecting MPI compilers"

	testprogram = "int main(int argc, char** argv) { return 0; }"

	conf.env.MPICXX = conf.find_program(Options.options.MPICXX, mandatory=True)
	conf.env.MPIF90 = conf.find_program(Options.options.MPIF90, mandatory=True)
	conf.env.MPIRUN = conf.find_program(Options.options.MPIRUN, mandatory=False)
	conf.env.MPIRUN_FLAG_NUMPROCS = Options.options.MPIRUN_FLAG_NUMPROCS
	#TODO add compile tests


def options_pyste(opt):
	opt = opt.parser.add_option_group("pyste options")
	pyste = os.path.abspath("./extern/pyste/pyste.py")
	opt.add_option("--gccxml", action="store", dest="GCCXML", default="gccxml", help="Command to run gccxml")
	opt.add_option("--gccxml-flags", action="store", dest="GCCXML_FLAGS", default="", help="Flags for gccxml")
	opt.add_option("--pyste",  action="store", dest="PYSTE", default=pyste, help="Command to run pyste")
	opt.add_option("--pyste-include", action="store", dest="PYSTE_INCLUDE", default="/usr/lib/openmpi/include", help="Extra include paths for pyste")
	opt.add_option("--enable-pyste", action="store_true", dest="PysteEnabled", default=False, help="Pyste is only run if this option is supplied to the build command")

@conftest
def check_pyste(conf):
	print "  - Detecting pyste"
	conf.env.GCCXML = conf.find_program(Options.options.GCCXML, mandatory=True)
	conf.env.GCCXML_FLAGS = Options.options.GCCXML_FLAGS
	curpath = os.path.abspath(os.path.curdir)
	conf.env.PYSTE = conf.find_program(Options.options.PYSTE, mandatory=True)
	conf.env.PYSTE_INCLUDE = split_names(Options.options.PYSTE_INCLUDE)

def options_gsl(opt):
	add_options_lib(opt, "gsl", "gsl")

@conftest
def check_gsl(conf):
	print "  - Detecting gsl"
	set_options_lib(conf, "gsl")


def is_module_enabled(self, moduleName):
	"""
	Checks whether a module is enabled or disabled from the command line

	"""
	enabled = True

	#defaults to the options set at config time
	enableModules = self.env.EnableModules
	disableModules = self.env.DisableModules

	#overridden now?
	if Options.options.EnableModules != "":
		enableModules = Options.options.EnableModules
	if Options.options.DisableModules != "":
		disableModules = Options.options.DisableModules


	# If enable-modules is specified, only the explicitly enabled modules
	# are enabled
	if enableModules != "":
		enabled = moduleName in enableModules.split()

	# If the module is mentioned in disable-modules, it is disabled regardless
	if moduleName in disableModules.split():
		enabled = False

	return enabled


Build.BuildContext.is_module_enabled = is_module_enabled
