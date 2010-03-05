
class InterruptHandlerClass:
	"""
	KeyboardInterrupt handler, which is used in Problem.Advance, in order to not
	to put the problem in an invalid state when using ctrl+c during propagation
	"""
	def __init__(self):
		self.__interrupt = False
		self.__signum = -1
		self.__frame = -1
		self.__origHandler = None

	def Register(self):
		if self.__origHandler != None:
			raise "Already registered"

		self.__origHandler = signal.signal(signal.SIGINT, self.Handler)

	def UnRegister(self):
		if self.__origHandler == None:
			return False
			
		signal.signal(signal.SIGINT, signal.default_int_handler)
		self.__origHandler = None
		self.__signum = -1
		self.__frame = -1
		self.__interrupt = False
	
	def Handler(self, signum, frame):
		print "Got keyboard interrupt, will terminate next timestep"
		self.__signum = signum
		self.__frame = frame
		self.__interrupt = True

	def ProcessInterrupt(self):
		if not self.IsInterrupted():
			print "HM? Not interrupted"
			return

		hndlr = self.__origHandler
		signum = self.__signum
		frame = self.__frame
		self.UnRegister()
		hndlr(signum, frame)

	def IsInterrupted(self):
		return self.__interrupt
	
InterruptHandler = InterruptHandlerClass()
