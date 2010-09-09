"""
Logger
======

Logging utilities for use in any Pyprop project/module

"""
import logging
import inspect

def GetClassLogger(obj):
	"""Return a logger for 'obj'

	Using attribute __module__ and __name__ of 'obj.__class__', create
	a logger with full path name.
	"""
	cls = obj.__class__
	return logging.getLogger("%s.%s" % (cls.__module__, cls.__name__))


def GetFunctionLogger():
	"""Return a logger for calling function
	
	Using the 'inspect' module, determine information about the
	calling function, and create the appropriate logger.
	"""
	callerFrame = inspect.stack()[1]
	moduleName = inspect.getmodule(callerFrame[0]).__name__
	functionName = callerFrame[3]
	return logging.getLogger("%s.%s" % (moduleName, functionName))
