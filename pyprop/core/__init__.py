import sys
from libcore import *
from libredirect import *

try:
	from libexpokit import *
except:
	#print "Warning: could not load EXPOKIT wrapper (%s)" % sys.exc_info()[1]
	pass
try:
	from libarpack import *
except:
	#print "Warning: could not load ARPACK wrapper (%s)" % sys.exc_info()[1]
	pass
try:
	from libpiram import *
except:
	print "Warning: could not load pIRAM wrapper (%s)" % sys.exc_info()[1]
try:
	from libode import *
except:
	print "Warning: could not load ODE wrapper (%s)" % sys.exc_info()[1]


def EnumerateRankClasses(baseName):
	for className, classObject in dict(globals()).iteritems():
		if className.startswith(baseName + '_'):
			rank = int(className[len(baseName)+1:])
			yield classObject
			
#Representation
def RepresentationHash(self):
	return self.GetId()

def RepresentationEq(self,other):
	if hasattr(other, "GetId"):
		return self.GetId() == other.GetId()
	else:
		return False
	
for classObject in EnumerateRankClasses('Representation'):
	classObject.__hash__ = RepresentationHash
	classObject.__eq__ = RepresentationEq
	
#CartesianRange
def CartesianRange__str__(self):
	mystr = "(min, max, count, step, istranslated) = ("
	mystr += str(self.Min) + ", "
	mystr += str(self.Max) + ", "
	mystr += str(self.Count) + ", "
	mystr += str(self.Dx) + ", "
	mystr += str(self.TranslatedGrid) + ")"
	return mystr
	
CartesianRange.__str__ = CartesianRange__str__
