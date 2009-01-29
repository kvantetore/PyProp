import sys

execfile("example.py")

pyprop.PrintOut("function = %s" % sys.argv[1])
pyprop.PrintOut("args = %s" sys.argv[2])
pyprop.PrintOut("     + %s" sys.argv[3])


