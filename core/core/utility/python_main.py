import sys
import os
import os.path

"""
Generates python_loader.cpp, which initializes the pyprop python modules
"""

fileContent = \
"""
#include "%(utility_path)sboostpythonhack.h"
#include <iostream>

using namespace std;

%(declare_modules)s

void LoadPypropModules()
{
  %(load_modules)s
}

int main(int argc, char** argv)
{
	LoadPypropModules();
	Py_Main(argc, argv);
}
								
"""


def is_library(name):
	"""
	Checks if a filename is a name to a library
	"""
	return name.startswith("lib") and name.endswith(".a")

path, name = os.path.split(__file__)
pypropPath = os.path.abspath(path + "/../..")
libPath = os.path.join(pypropPath, "core/lib")

if len(sys.argv) < 2:
	raise Exception("Please specify outputfile, i.e.\n    python python_main.py python_main.cpp ")
outputFile = sys.argv[1]
	
#Get base name (without ext) of all libs in 
moduleList = [os.path.splitext(name)[0] for name in os.listdir(libPath) if is_library(name)]
#add all libs specified on the command line
moduleList += sys.argv[2:]

load_modules = "\n".join(['PyImport_AppendInittab("%s", init%s);' % (name, name) for name in moduleList])
declare_modules = "\n".join(["PyMODINIT_FUNC init%s();" % name for name in moduleList])
utility_path = path + "/"

f = file(outputFile, "w")
try:
	f.write(fileContent % locals())
finally:
	f.close()


