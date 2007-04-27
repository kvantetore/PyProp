import numpy

def SavePickleArray(fb, array):
	autoClose = False
	if fb.__class__.__name__ == 'str':
		fb = file(fb, 'w')
		autoClose = True
		
	pickle.dump(array, fb, 2)
	
	if (autoClose):
		fb.close()

def LoadPickleArray(fb):
	autoClose = False
	if fb.__class__.__name__ == 'str':
		fb = file(fb, 'r')
		autoClose = True
		
	array = pickle.load(fb)
	
	if (autoClose):
		fb.close()

	return array
	
