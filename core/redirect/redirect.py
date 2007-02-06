import libredirect as redir

cout_redirect = None
stdout_redirect = None

class RedirectedOutput:
	def __init__(self, stdout):
		self.stdout = stdout

	def write(self, *data):
		self.stdout.write("Redirected cout: ")
		self.stdout.write(*data)

import sys

stdout_redirect = RedirectedOutput(sys.stdout)
cout_redirect = redir.redirect_cout()
sys.stdout = myout

#print with redirected stdout
print "Output from Python"
redir.test()

redir.restore_cout(cout_redirect)
sys.stdout = sys.__stdout__

#print with original stdout
print "Output from Python"
redir.test()


