import sys

import core
from distribution import ProcCount, ProcId, pypar


DEBUG_PRINT_MEMORY_USAGE = False

if hasattr(core, "PapiSetup"):
	core.PapiSetup()

def GetMemoryUsage():
	if hasattr(core, "PapiGetDynamicMemoryInfo"):
		mem = core.PapiGetDynamicMemoryInfo()
	return mem.size, mem.peak


def SerialPrint(str, proc=-1):
	if ProcCount == 1:
		print str
	else:
		if proc==-1: procList = range(ProcCount)
		else: procList = [proc]
		for i in procList:
			if i == ProcId:
				print "Proc %4i: %s" % (ProcId, str,)
			sys.stdout.flush()
			pypar.barrier()

def PrintMemoryUsage(header, proc=-1):
	if DEBUG_PRINT_MEMORY_USAGE:
		mem, peak = GetMemoryUsage()
		SerialPrint(header, proc=0)
		SerialPrint("Memory Usage = %10iKB, Peak = %10iKB" % (mem, peak), proc=proc)



