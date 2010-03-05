import sys
import os

#Try to load trilinos modules with ctypes
#for some reason they must be loaded with RTLD_GLOBAL
#for pyprop with trilinos to work.
try:
	import ctypes
	ctypes.CDLL("libepetra.so", ctypes.RTLD_GLOBAL)
	ctypes.CDLL("libteuchos.so", ctypes.RTLD_GLOBAL)
	ctypes.CDLL("libifpack.so", ctypes.RTLD_GLOBAL)
	LOAD_TRILINOS_OK = True
except:
	LOAD_TRILINOS_FALSE = True


#find out if we're an install path or source path
import pyprop
if pyprop.IsRunningFromSource:
	sys.path.append(os.path.join(pyprop.BuildPath, "core"))
else:
	sys.path.append(os.path.join(__path__[0], "core"))

#import core c++ module
try:
	from libcore import *
	LOAD_CORE_OK = True
except:
	LOAD_CORE_OK = False


def EnumerateRankClasses(baseName):
	"""
	List all classes that starts with baseName + "_"

	This can be used to find all rank-instances of a certain
	function or class
	"""
	for className, classObject in dict(globals()).iteritems():
		if className.startswith(baseName + '_'):
			rank = int(className[len(baseName)+1:])
			yield classObject
			
#Representation
def RepresentationHash(self):
	"""
	Hack to get equality checks to work on representations
	"""
	return self.GetId()

def RepresentationEq(self,other):
	"""
	Hack to get equality checks to work on representations
	"""
	if hasattr(other, "GetId"):
		return self.GetId() == other.GetId()
	else:
		return False

if LOAD_CORE_OK:
	#Hack to get equality checks to work on representations
	for classObject in EnumerateRankClasses('Representation'):
		classObject.__hash__ = RepresentationHash
		classObject.__eq__ = RepresentationEq
		
##CartesianRange
#def CartesianRange__str__(self):
#	mystr = "(min, max, count, step, istranslated) = ("
#	mystr += str(self.Min) + ", "
#	mystr += str(self.Max) + ", "
#	mystr += str(self.Count) + ", "
#	mystr += str(self.Dx) + ", "
#	mystr += str(self.TranslatedGrid) + ")"
#	return mystr
#	
#CartesianRange.__str__ = CartesianRange__str__
