from libcore import *
from libredirect import *

def EnumerateRankClasses(baseName):
	for className, classObject in dict(globals()).iteritems():
		if className.startswith(baseName + '_'):
			rank = int(className[len(baseName)+1:])
			yield classObject
			
#Representation
def RepresentationHash(self):
	return self.GetId()

def RepresentationEq(self,other):
	return self.GetId() == other.GetId()
	
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
