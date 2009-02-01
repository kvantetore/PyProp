import sys
import pypar
execfile("example.py")

pyprop.PrintOut("function = %s" % sys.argv[1])
pyprop.PrintOut("args = %s" % sys.argv[2])
pyprop.PrintOut("     + %s" % sys.argv[3])

function = eval(sys.argv[1])
arglist = eval(sys.argv[2])
argdict = eval(sys.argv[3])

result = function(*arglist, **argdict)

print result

pypar.barrier()
pypar.finalize()
