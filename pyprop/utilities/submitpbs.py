import os
import sys
import commands

from datetime import timedelta

"""
Utilities for submitting jobs to the PBS scheduling system used
on fimm.bccs.uib.no.
"""


class SubmitScript:
	#resources
	walltime = timedelta(hours=0, minutes=30, seconds=0)
	nodes = 1
	ppn = 1
	proc_memory = "500mb"
	
	account = "matematisk"
	jobname = "myjob"

	stdout = None
	stdin = None
	stderr = None

	executable = "a.out"
	parameters = ""
	workingdir = None

	compiler = "gnu"
	mpi = "openmpi"
	mpirun_enabled = True

	compilerscript = "~/use_compiler"
	mpiscript = "~/use_mpi"

	#Environment variables to be copied from current environment
	env_copy = list()
	#Extra environment variables to add to 
	env_extra = dict()


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
		script.append("#PBS -l nodes=" + str(self.nodes) + ":ppn=" + str(self.ppn))
		script.append("#PBS -l pmem=" + str(self.proc_memory))

		#Administrative
		script.append('#PBS -N "' + str(self.jobname) + '"')
		script.append("#PBS -A " + self.account)

		#Environment variables
		for env in self.env_copy:
			script.append("#PBS -v " + str(env))

		for env in self.env_extra:
			script.append("#PBS -v " + str(env) + str(self.env_extra[env]))

		#IO redirection
		if self.stdout != None:
			script.append('#PBS -o "' + str(self.stdout) + '"')
		if self.stderr != None:
			script.append('#PBS -e "' + str(self.stderr) + '"')

		#Working dir
		if self.workingdir == None:
			self.workingdir = os.path.abspath(os.curdir)
		script.append("cd " + str(self.workingdir))

		#Load the correct paths for the given mpi and compiler
		script.append("source " + str(self.compilerscript)  + " " + str(self.compiler))
		script.append("source " + str(self.mpiscript) + " " + str(self.mpi))
		
		#Check if we're redirecting stdin
		instr = ""
		if self.stdin != None:
			instr = "< " + str(self.stdin)
	
		#HACK: if we're using openmpi, we must specify number of procs
		procstr = ""
		if self.mpi == "openmpi":
			procstr = "-np " + str(self.nodes * self.ppn)

		#Check if we're using mpirun or mpiexec to start our program
		mpirun = ""
		if self.mpirun_enabled:
			mpirun = "$MPIRUN " + procstr + " "

		#Create script line
		script.append(mpirun + str(self.executable) + " " + str(self.parameters) + instr)
		
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
	
