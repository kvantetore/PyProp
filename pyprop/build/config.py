import os.path

class Config(object):
	"""
	Configuration object populated by the various check-functions
	depending on the actual configuration
	"""
	mpi_enabled = True
	compiler_cxx = "mpicxx"
	compiler_ftn = "mpif90"
	python_module_mode = "shared"

	library = []
	library_path = []
	include_path = [".", "core"]
	defines = []

	


def check_blas_lapack(conf):
	print "Deteting BLAS/LAPACK"
	conf.library += ["mkl", "guide", "mkl_core", "mkl_lapack", "pthread",]
	conf.library_path += ["/opt/intel/mkl/10.0.1.014/lib",]
	conf.include_path += ["/opt/intel/mkl/10.0.1.014/include",]
	#TODO add compile tests


def check_python(conf):
	print "Detecting python"
	import distutils.sysconfig;
	conf.library += []
	conf.library_path += [distutils.sysconfig.get_python_lib()]
	conf.include_path += [distutils.sysconfig.get_python_inc()]
	#TODO add compile tests


def check_boost_python(conf):
	print "Detecting boost::python"
	conf.library += ["boost_python-mt"]
	conf.library_path += ["/opt/boost/lib"]
	conf.include_path += ["/opt/boost/include"]
	#TODO add compile tests


def check_blitz(conf):
	print "Detecting blitz"
	blitz_path = "/home/torebi/prog/pyprop/extern/blitz/build"
	conf.library += ["blitz"]
	conf.library_path += [os.path.join(blitz_path, "lib")]
	conf.include_path += [os.path.join(blitz_path, "include")]


def check_trilinos(conf):
	print "Detecting Trilinos"
	trilinos_path = "/opt/trilinos"
	conf.library += ["epetra", "teuchos", "ifpack", "anasazi"]
	conf.library_path += [os.path.join(trilinos_path, "lib")]
	conf.include_path += [os.path.join(trilinos_path, "include")]
	conf.defines += ["PYPROP_USE_TRILINOS"]


def setup(configFile=None):
	conf = Config()
	check_blas_lapack(conf)
	check_python(conf)
	check_boost_python(conf)
	check_blitz(conf)
	check_trilinos(conf)
	return conf

conf = setup()
