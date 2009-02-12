from pyprop.utilities.folderwatch import FolderWatch, CHANGE_NEW, CHANGE_MODIFIED
import os
import tables

from pyprop.serialization import RemoveExistingDataset

def FindEigenvaluesNearShift(shift, **args):
	#Setup output files
	messageFile = args["messageFile"]
	eigenvectorFile = args["eigenvectorFile"]
	pyprop.PrintOut("Message File = %s" % messageFile)

	#Setup Problem
	prop = SetupProblem(silent=True, eigenvalueShift=shift, **args)

	#Setup shift invert solver in order to perform inverse iterations
	shiftInvertSolver = pyprop.GMRESShiftInvertSolver(prop)
	prop.Config.Arpack.inverse_iterations = True
	prop.Config.Arpack.matrix_vector_func = shiftInvertSolver.InverseIterations

	#Setup eiganvalue solver
	solver = pyprop.PiramSolver(prop)
	solver.Solve()

	#save eigenvectors
	SaveEigenvalueSolverShiftInvert(eigenvectorFile, solver, shift)

	#Get the converged eigenvalues
	ev = solver.Solver.GetEigenvalues().copy()
	estimates = solver.Solver.GetConvergenceEstimates().copy()
	idx = where(estimates < 0)[0]
	eigenvalues = ev[idx]

	print "SHIFT = ", prop.Config.GMRES.shift

	#convert from shift inverted eigenvalues to "actual" eigenvalues
	eigenvalues = 1.0 / eigenvalues + shift

	shiftInvertSolver.PrintStatistics()

	print "GMRESERRORS = ", shiftInvertSolver.Solver.GetErrorEstimateList()
	PrintOut("EIGENVALUES = %s" % eigenvalues)

	sortIdx = argsort(real(eigenvalues))
	eigenvalues = eigenvalues[sortIdx]

	#Save message file
	if pyprop.ProcId == 0:
		f = tables.openFile(messageFile, "w")
		f.createArray(f.root, "eigenvalues", eigenvalues)
		f.close()

	return solver, shiftInvertSolver

def SaveEigenvalueSolverShiftInvert(filename, solver, shift):

	#Get eigenvalue error estimate
	errorEstimatesPIRAM = solver.Solver.GetErrorEstimates()
	convergenceEstimatesEig = solver.Solver.GetConvergenceEstimates()

	#Get eigenvalues
	prop = solver.BaseProblem
	E = 1.0 / array(solver.GetEigenvalues()) + shift

	#Store eigenvalues and eigenvectors
	PrintOut("Now storing eigenvectors...")
	for i in range(len(E)):
		solver.SetEigenvector(prop.psi, i)
		prop.SaveWavefunctionHDF(filename, "/Eig/Eigenvector%03i" % i)

	if pyprop.ProcId == 0:
		RemoveExistingDataset(filename, "/Eig/Eigenvalues")
		RemoveExistingDataset(filename, "/Eig/ErrorEstimateListGMRES")
		RemoveExistingDataset(filename, "/Eig/ErrorEstimateListPIRAM")
		RemoveExistingDataset(filename, "/Eig/ConvergenceEstimateEig")
		h5file = tables.openFile(filename, "r+")
		try:
			#myGroup = h5file.createGroup("/", "Eig")
			myGroup = h5file.getNode("/Eig")
			h5file.createArray(myGroup, "Eigenvalues", E)
			#h5file.createArray(myGroup, "ErrorEstimateListGMRES", errorEstimatesGMRES)
			h5file.createArray(myGroup, "ErrorEstimateListPIRAM", errorEstimatesPIRAM)
			h5file.createArray(myGroup, "ConvergenceEstimateEig", convergenceEstimatesEig)

			#Store config
			myGroup._v_attrs.configObject = prop.Config.cfgObj
			
			#PIRAM stats
			myGroup._v_attrs.opCount = solver.Solver.GetOperatorCount()
			myGroup._v_attrs.restartCount = solver.Solver.GetRestartCount()
			myGroup._v_attrs.orthCount = solver.Solver.GetOrthogonalizationCount()
		finally:
			h5file.close()


def RunSubmit(function, procCount=1, procPerNode=4, *arglist, **argdict):
	if isinstance(function, str):
		arg1 = function
	else:
		arg1 = function.func_name

	arg2 = commands.mkarg(repr(arglist))
	arg3 = commands.mkarg(repr(argdict))

	jobId = None
	if INSTALLATION == "hexagon":
		submit = submitpbs.SubmitScript()
		submit.procs = procCount
		submit.ppn = min(procPerNode, procCount)
		submit.executable = "./python-exec run-function.py"
		submit.parameters = arg1 + arg2 + arg3
		submit.WriteScript("test.job")
		jobId = submit.Submit()

	elif INSTALLATION == "stallo":
		raise Exception("please to implement")
	
	elif INSTALLATION == "local":
		raise Exception("please to implement")
	
	else:
		raise Exception("Unknown installation '%s'" % INSTALLATION)

	return jobId

class SpectrumFinder(object):
	"""
	Class for finding all eigenvalues/eigenvectors in a 
	range of eigenvalues. This is done by running a series
	of FindEigenvaluesNearShift. Each call returns a number
	of eigenvalues near a given shift. The intervals between 
	startEigenvalue and endEigenvalue are subdivided until 
	the eigenvalues from two nearby shifts overlap. 

	This algorithm has one serious flaw: we might not find all
	eigenvalues, if the edge-eiggenvalues between two shifts are
	degenerate.

	Another inefficiency is that as we do not try to guess the
	density of states, the final subdivision may cause the eigenvalues
	of two nearby shifts to overlap severly.
	
	"""

	class Message(object):
		def __init__(self, shift, interval, messageId, startTime, msgFolder, dataFolder):
			self.MessagesFolder = msgFolder
			self.DataFolder = dataFolder
			self.Shift = shift
			self.MessageId = messageId
			self.StartTime = startTime
			self.Interval = interval

		def GetMessageFile(self):
			return os.path.join(self.MessagesFolder, "message_%s.h5" % (self.MessageId))
			
		def GetDataFile(self):
			return os.path.join(self.DataFolder, "eigenvectors_%s.h5" % (self.MessageId))
		
	def __init__(self, startEigenvalue, endEigenvalue, folderPrefix, **args):
		self.StartEigenvalue = startEigenvalue
		self.EndEigenvalue = endEigenvalue
		self.Arguments = args

		#Setup folders
		self.MessagesFolder = os.path.join(folderPrefix, "messages")
		self.DataFolder = os.path.join(folderPrefix, "data")
		self.StatusFolder = folderPrefix
		if not os.path.exists(self.MessagesFolder):
			os.makedirs(self.MessagesFolder)
		if not os.path.exists(self.DataFolder):
			os.makedirs(self.DataFolder)
		if not os.path.exists(self.StatusFolder):
			os.makedirs(self.StatusFolder)

		procCount = args.get("procCount", None)
		#if no proccount is specified, use one for each angular basis function
		if procCount == None:
			conf = SetupConfig(**args)
			lcount = len([1 for l in conf.AngularRepresentation.index_iterator])
			self.ProcCount = lcount
		else:
			self.ProcCount = procCount
		self.ProcPerNode = 4

		#Setup maps to keep track of which eigenvalues we have found
		self.NextMessageId = 0
		self.ActiveMessages = []  #Messages we are waiting for to arrive
		self.CompletedMessages = []     #List over completed messages 
		self.ShiftMap = {}		  #Maps shifts to a list of eigenvalues
		self.ActiveIntervals = []
		self.ErrorMessages = []
	
	def GetNextMessageId(self):
		msg = self.NextMessageId
		self.NextMessageId += 1
		return msg

	def StartMessage(self, shift, interval, messageId):
		#Create message
		message = self.Message(shift, interval, messageId, time.time(), self.MessagesFolder, self.DataFolder)
		messageFile = message.GetMessageFile()
		eigenvectorFile = message.GetDataFile()

		#Remove files if it exists
		if os.path.exists(messageFile):
			os.remove(messageFile)
		if os.path.exists(eigenvectorFile):
			os.remove(eigenvectorFile)

		#Start job
		jobId = RunSubmit(FindEigenvaluesNearShift, self.ProcCount, self.ProcPerNode, shift=shift, \
			messageFile=messageFile, eigenvectorFile=eigenvectorFile, **self.Arguments)
		message.JobId = jobId
		self.ActiveMessages.append(message)

	def WaitForCompletedMessages(self):
		timeout = False
		pyprop.PrintOut("Waiting for %i messages to complete" % len(self.ActiveMessages))
		while len(self.ActiveMessages) > 0 and not timeout:
			completedMessages = []
			errorMessages = []

			#check for finished messages
			for msg in list(self.ActiveMessages):
				if submitpbs.CheckJobCompleted(msg.JobId):
					msgPath = os.path.abspath(msg.GetMessageFile())
					if os.path.exists(msgPath):
						completedMessages.append(msg)
					else:
						print "Job Status for error message: "
						print submitpbs.GetJobStatus(msg.JobId)
						self.ErrorMessages.append(msg)
						self.ActiveMessages.remove(msg)

			#process completed messages
			map(self.MessageCompleted, completedMessages)
		
			if len(completedMessages) == 0:
				time.sleep(10)

	def MessageCompleted(self, msg):
		#get eigenvalues
		f = tables.openFile(msg.GetMessageFile())
		try:
			eigenvalues = f.root.eigenvalues[:]
		finally:
			f.close()

		#Update message lists
		self.ActiveMessages.remove(msg)
		self.CompletedMessages.append(msg)

		pyprop.PrintOut("Found %i eigenvalues around %f" % (len(eigenvalues), msg.Shift))
		self.ShiftMap[msg.Shift] = eigenvalues

		#Check if this is one of the first two eigenvalues
		if msg.Interval == None:
			#if this is message the last of the boundary messages
			if len(self.ShiftMap) == 2:
				self.ActiveIntervals.append((self.StartEigenvalue, self.EndEigenvalue))

		else:
			#Otherwise subdivide the interval and mark it active
			start, end = msg.Interval
			shift = msg.Shift
			self.ActiveIntervals.append((start, shift))
			self.ActiveIntervals.append((shift, end))

	def RunIteration(self):
		#if this is the first iteration:
		if len(self.ShiftMap) == 0:
			self.StartMessage(self.StartEigenvalue, None, self.GetNextMessageId())
			self.StartMessage(self.EndEigenvalue, None, self.GetNextMessageId())

		for start, end in self.ActiveIntervals:
			evStart = self.ShiftMap[start]
			evEnd = self.ShiftMap[end]

			#Check whether there can be a gab between the last eigenvalue of start 
			#and the smallest eigenvalue of end
			if evStart[-1] < evEnd[0]:
				shift = (end + start) / 2.
				print "Starting message at shift %s" % shift
				self.StartMessage(shift, (start, end), self.GetNextMessageId())
		self.ActiveIntervals = []

		self.WaitForCompletedMessages()
	
	def FilterActiveIntervals(self):
		def intervalNeedsSubdivision(interval):
			start, end = interval
			evStart = self.ShiftMap[start]
			evEnd = self.ShiftMap[end]
		
			#Check whether there can be a gab between the last eigenvalue of start 
			#and the smallest eigenvalue of end
			return evStart[-1] < evEnd[0]

		return filter(intervalNeedsSubdivision, self.ActiveIntervals)

	def RunToEnd(self):
		self.RunIteration()
		while len(self.ActiveIntervals) > 0:
			self.RunIteration()

	def LoadAvailableMessages(self):
		dataFiles = os.listdir(self.DataFolder)
		dataFiles = filter(lambda x: x.startswith("eigenvectors_"), dataFiles)
		dataFiles = map(lambda x: os.path.join(self.DataFolder, x), dataFiles)

		def getMessage(dataFile):
			messageId = int(os.path.splitext(dataFile)[0].split("_")[1])
			startTime = os.path.getctime(dataFile)
			f = tables.openFile(dataFile, "r")
			try:
				shift = float(f.root.Eig._v_attrs.configObject.get("GMRES", "shift"))
				eigenvalues = sorted(f.root.Eig.Eigenvalues[:])
			finally:
				f.close()

			message = self.Message(shift, (-Inf, Inf), messageId, startTime, self.MessagesFolder, self.DataFolder)
			return message, eigenvalues

		#Get messages and eigenvalues from datafiles
		first = lambda x: x[0]
		messageList, eigenvalueList = zip(*map(lambda x: getMessage(x), dataFiles))

		#Setup Completed Messages list
		self.CompletedMessages += messageList
		#Setup ShiftMap
		for msg, ev in zip(messageList, eigenvalueList):
			self.ShiftMap[msg.Shift] = ev
		#Setup ActiveIntervals
		shiftList = sorted(self.ShiftMap.keys())
		self.ActiveIntervals += [(start, end) for start, end in zip(shiftList[0:-1], shiftList[1:])]
		#next message id
		self.NextMessageId = 1 + max(map(lambda x: x.MessageId, self.CompletedMessages))

	def CancelMessage(self, message):
		pyprop.PrintOut("Canceling message for shift = %s" % message.Shift)
		if message not in self.ActiveMessages:
			raise Exception("Message not active")
		#Put the interval back into active intervals
		self.ActiveIntervals.append(message.Interval)
		#Remove message from active messages
		self.ActiveMessages.remove(message)

	def CancelAllMessages(self):
		for msg in list(self.ActiveMessages):
			self.CancelMessage(msg)
	
	def GetSortedEigenvalues(self):
		eigenvalueList = []
		overlappingEigenvalues = 0
		for shift in sorted(self.ShiftMap.keys()):
			for curE in self.ShiftMap[shift]:
				if len(eigenvalueList) == 0 or curE > eigenvalueList[-1]:
					eigenvalueList.append(curE)
				else:
					overlappingEigenvalues += 1
	
		print "Eiganvalues found       = %i" % (len(eigenvalueList))
		print "Overlapping eigenvalues = %i" % (overlappingEigenvalues)
	
		return eigenvalueList

	def GetSortedEigenvectors(self):
		eigenvalues = self.GetSortedEigenvalues()
		eigenvectorFile = os.path.join(self.StatusFolder, "eigenvectors.h5")

		shiftToMessage = dict([(msg.Shift, msg) for msg in self.CompletedMessages])

		f = tables.openFile(eigenvectorFile, "w")
		try:
			eigGroup = f.createGroup(f.root, "Eig")
		
			#Store metadata from this spetrum finder
			attrs = eigGroup._v_attrs
			attrs.StartEigenvalue = self.StartEigenvalue
			attrs.EndEigenvalue = self.EndEigenvalue
			attrs.Arguments = self.Arguments
			attrs.FolderPrefix = self.StatusFolder
		
			index = 0
			eigenvalueList = []
			for shift in sorted(self.ShiftMap.keys()):
				msg = shiftToMessage[shift]
				curFile = tables.openFile(msg.GetDataFile(), "r")
				try:
					curEigenvalues = curFile.Eig.Eigenvalues[:]
					eigIndex = argsort(curEigenvalues)
					for curIndex, curE in zip(eigIndex, curEigenvalues[eigIndex]):
						if index == 0 or curE > eigenvalueList[-1]:
							infoStr =  "Progress: %3i%% (%i/%i)" % ((index * 100)/ len(eigenvalues))
							sys.stdout.write("\b"*30 + infoStr)
							sys.stdout.flush()

							curNode = curFile.getNode("/Eig/Eigenvector%03i" % curIndex)[:]
							newNode = f.createArray(eigGroup, "Eigenvector_%i" % index, curData)
							newNode._v_attrs.configObject = curFile.Eig._v_attrs.configObject

							index += 1

				finally:
					curFile.close()
						
		finally:
			f.close()



