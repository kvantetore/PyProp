import time

class Timer(object):
	def __init__(self):
		self.Duration = 0
		self.IsRunning = False

	def Start(self):
		if self.IsRunning:
			raise Exception("Timer already running")
		self.Duration -= time.time()
		self.IsRunning = True

	def Stop(self):
		if not self.IsRunning:
			raise Exception("Timer not running")

		self.Duration += time.time()
		self.IsRunning = False

		return self.Duration

	def GetDuration(self):
		return self.Duration

	def __str__(self):
		duration = int(self.Duration)
		h = (duration / 3600)
		m = (duration / 60) % 60
		s = self.Duration - h * 3600 - m * 60
		
		str = []
		if h>0: str += ["%ih" % h]
		if m>0: str += ["%im" % m]
		if s>0: str += ["%fs" % s]
		
		return " ".join(str)



class Timers(dict):
	"""
	Dictonary of Timer objects. Automatically instantiates new
	Timer objects on first use of a key

	timers = Timers()
	timers["Total"].Start()
	#do stuff...
	Timers[Total"].Stop()

	print Timers
	"""
	def __getitem__(self, ident):
		if not self.has_key(ident):
			self[ident] = Timer()
		return dict.__getitem__(self, ident)

	def __setitem__(self, ident, value):
		if not isinstance(value, Timer):
			raise ValueError("Only Timer-objects are allowed as value in timers")
		dict.__setitem__(self, ident, value)

	def __str__(self):
		str =  "Timers:\n"
		for key, timer in self.iteritems():
			str += "    %s = %s\n" % (key, timer)
		return str
	
class Counters(dict):
	"""
	Dictionary of counters
	"""
	def __getitem__(self, ident):
		if not self.has_key(ident):
			self[ident] = 0
		return dict.__getitem__(self, ident)

	def __setitem__(self, ident, value):
		if not isinstance(value, int):
			raise ValueError("Only integer-objects are allowed as value in counters")
		dict.__setitem__(self, ident, value)

	def __str__(self):
		str =  "Counters:\n"
		for key, counter in self.iteritems():
			str += "    %s = %s\n" % (key, counter)
		return str
	
		
