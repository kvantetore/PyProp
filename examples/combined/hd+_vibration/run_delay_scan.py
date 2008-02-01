#!/usr/bin/env python

execfile("example.py")

argstring = " ".join(sys.argv[1:])
args = eval(argstring)
RunDelayScan(**args)

