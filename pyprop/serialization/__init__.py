try:
	import tables
except:
	print "Warning: Could not load module tables (pytables). Serialization to HDF files will not be available"

execfile(__path__[0] + "/Pickle.py")
execfile(__path__[0] + "/CommonHDF.py")
execfile(__path__[0] + "/TensorPotentialHDF.py")
execfile(__path__[0] + "/WavefunctionHDF.py")

