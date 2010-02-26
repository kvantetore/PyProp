import os
import sys
import commands

from datetime import timedelta
from math import ceil

"""
Utilities for submitting jobs to the PBS scheduling system used
on stallo.uit.no.
"""


class SubmitScript:
	#resources
	walltime = timedelta(hours=0, minutes=30, seconds=0)
	nodes = 1
	ppn = 1
	proc_memory = None 
	
	account = None 
	jobname = "pyprop"

	stdout = None
	stdin = None
	stderr = None

	executable = "a.out"
	parameters = ""
	workingdir = None

	interconnect = "ib" # "ib" or None
	

	def CreateScript(self):
		script = list()

		#preamble
		script.append("#! /bin/bash -")
		script.append("#PBS -S /bin/bash")

		#Resources
		#time
		hours = self.walltime.days * 24 + (self.walltime.seconds / 3600)
		minutes = (self.walltime.seconds / 60) % 60
		seconds = self.walltime.seconds % 60
		script.append("#PBS -l walltime=%i:%i:%i" % (hours, minutes, seconds))

		if hasattr(self, "procs"):
			procCount = self.procs
			self.nodes = int(ceil(procCount / float(self.ppn)))
		else:
			procCount = self.nodes*self.ppn

		#procs
		if self.interconnect == "ib":
			interconnect = ":ib"
		else:
			interconnect = ":gige"
		script.append("#PBS -l nodes=%i:ppn=%i%s" % (self.nodes, self.ppn, interconnect))

		#mem
		if self.proc_memory != None:
			script.append("#PBS -l pmem=" + str(self.proc_memory))

		#Administrative
		script.append('#PBS -N "%s"' % (self.jobname,))
		if self.account != None:
			script.append("#PBS -A %s" % (self.account,))

		#IO redirection
		if self.stdout != None:
			script.append('#PBS -o "%s"' % (self.stdout))
		if self.stderr != None:
			script.append('#PBS -e "%s"' % (self.stderr))

		#Working dir
		if self.workingdir == None:
			self.workingdir = os.path.abspath(os.curdir)
		script.append("cd %s" % (self.workingdir))

		#Check if we're redirecting stdin
		instr = ""
		if self.stdin != None:
			instr = "< " + str(self.stdin)
	
		#Create script line
		script.append("mpirun -np %s " % procCount + str(self.executable) + " " + str(self.parameters) + instr)
		
		#exit stuff
		script.append("exit $?")

		return script


	def WriteScript(self, filename):
		script = self.CreateScript()
		f = file(filename, 'w')
		for line in script:
			f.write(line)
			f.write("\n")
		f.close()
		
	def Submit(self):
		#create a temporary file for the script
		tempName = os.tempnam(".", "scrpt")
		self.WriteScript(tempName)

		#submit script
		jobName = commands.getoutput("qsub " + tempName)

		#delete temporary script file
		os.remove(tempName)

		print jobName + " submitted"

		return jobName
	
