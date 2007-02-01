import numpy

def SaveMatlabArray(fb, array):
	array = asarray(array)
	if len(array.shape) > 2:
		raise "Cannot save matlab arrays with rank > 2."
		
	autoClose = False
	if fb.__class__.__name__ == 'str':
		fb = file(fb, 'w')
		autoClose = True
		
	if len(array.shape) == 1:
		for index, val in enumerate(array):
			fb.write(str(val))
			fb.write(" ")
		fb.write("\n")
		
	if len(array.shape) == 2:
		mrange = range(array.shape[0])
		for m in mrange:
			for index, val in enumerate(array[m,:]):
				fb.write(str(val))
				fb.write(" ")
			fb.write("\n")
	
	if (autoClose):
		fb.close()

def AppendMatlabArray(fileName, array):
	fb = file(fileName, 'a')
	SaveMatlabArray(fb, array)
	fb.close()

def LoadMatlabArray(fb):
	autoClose = False
	if fb.__class__.__name__ == 'str':
		fb = file(fb, 'r')
		autoClose = True
		
	initialRows = 100
	arr = None
	curRow = 0
	while True:
		line = fb.readline()
		if bool(line) == False:
			break
		
		row = numpy.fromstring(line.strip(), sep=' ', dtype=double)
			
		if arr == None:
			arr = empty( ((initialRows), len(row)), dtype=double )
			
		if curRow >= arr.shape[0]:
			arr.resize((arr.shape[0]+initialRows, arr.shape[1]))
		
		arr[curRow,:] = row
		curRow += 1
		
	arr.resize((curRow, arr.shape[1]))
	
	if (autoClose):
		fb.close()
		
	return arr

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
	
