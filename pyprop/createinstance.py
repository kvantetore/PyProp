import sys

def FindObjectStack(name, startdepth=0, maxdepth=100):
	"""
	Finds an object in locals or globals of one of the preceding
	stack frames. The closest stack frames are searched for a matching
	object first.
	"""
	for depth in range(startdepth+1,maxdepth):
		try:
			f = sys._getframe(depth)
		except ValueError:
			stacktrace = ""
			for i in range(depth-1, 0, -1):
				stacktrace += str(sys._getframe(i).f_code) + "\n"
			raise Exception("Could not find object ""%s"" in globals of the following frames: \n%s" % (name, stacktrace))

		try:
			obj = eval(name, f.f_globals, f.f_locals)
			return obj
		except Exception:
			pass

	raise Exception("Exceeded maxdepth (%s) while traversing stack" % maxdepth)


def CreateInstanceRank(className, rank):
	try:
		name = "%s_%i" % (className, rank)
		cls = FindObjectStack(name, startdepth=1)
		instance = cls()
		return instance
	except Exception:
		print "Could not create instance of class ", className, "with rank ", rank
		raise
	
