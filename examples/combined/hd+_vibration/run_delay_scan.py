#!/usr/bin/env python

execfile("example.py")

argstring = " ".join(sys.argv[1:])
args = eval(argstring)
delaySlice = args["delayList"]
args["delayList"] = r_[delaySlice]
RunDelayScan(**args)

