import os
import sys
import commands

from datetime import timedelta

"""
Utilities for submitting jobs to the PBS scheduling system used
on hexagon.bccs.uib.no.
"""


class SubmitScript:
	#resources
	walltime = timedelta(hours=0, minutes=30, seconds=0)
	nodes = 1
	ppn = 1
	proc_memory = None
	
	account = None
	jobname = "myjob"

	stdout = None
	stdin = None
	stderr = None

	executable = "a.out"
	parameters = ""
	workingdir = None

	def CreateScript(self):
		script = list()

		#preamble
		script.append("#! /bin/bash -")
		script.append("#PBS -S /bin/bash")

		#Resources
		hours = self.walltime.days * 24 + (self.walltime.seconds / 3600)
		minutes = (self.walltime.seconds / 60) % 60
		seconds = self.walltime.seconds % 60
		script.append("#PBS -l walltime=" + str(hours) + ":" + str(minutes) + ":" + str(seconds))
		script.append("#PBS -l mppwidth=" + str(self.nodes*self.ppn))
		if self.proc_memory != None:
			script.append("#PBS -l pmem=" + str(self.proc_memory))

		#Administrative
		if self.jobname != None:
			script.append('#PBS -N "' + str(self.jobname) + '"')
		if self.account != None:	
			script.append("#PBS -A " + self.account)

		#IO redirection
		if self.stdout != None:
			script.append('#PBS -o "' + str(self.stdout) + '"')
		if self.stderr != None:
			script.append('#PBS -e "' + str(self.stderr) + '"')

		#Working dir
		if self.workingdir == None:
			self.workingdir = os.path.abspath(os.curdir)
		script.append("cd " + str(self.workingdir))

		#Check if we're redirecting stdin
		instr = ""
		if self.stdin != None:
			instr = "< " + str(self.stdin)
	
		#Create script line
		script.append(str(self.executable) + " " + str(self.parameters) + instr)
		
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
	
