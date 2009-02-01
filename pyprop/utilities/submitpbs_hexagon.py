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
	jobname = "pyprop"

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

		if hasattr(self, "procs"):
			procCount = self.procs
		else:
			procCount = self.nodes*self.ppn

		#Resources
		hours = self.walltime.days * 24 + (self.walltime.seconds / 3600)
		minutes = (self.walltime.seconds / 60) % 60
		seconds = self.walltime.seconds % 60
		script.append("#PBS -l walltime=" + str(hours) + ":" + str(minutes) + ":" + str(seconds))
		script.append("#PBS -l mppwidth=" + str(procCount))
		script.append("#PBS -l mppnppn=" + str(self.ppn))
		if self.proc_memory != None:
			script.append("#PBS -l mppmem=" + str(self.proc_memory))

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

		#check if user supplied aprun 
		if not self.executable.lower().startswith("aprun "):
			mem = ""
			if self.proc_memory != None:
				memstr = "-m %s" % self.proc_memory
			self.executable = "aprun -n %i -N %i %s %s" % (procCount, self.ppn, mem, self.executable)
	
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


def GetJobStatus(jobId):
	"""
	Returns a dict containing the info from qstat -f jobid
	if the job is not found None is returned
	"""
	status, output = commands.getstatusoutput("qstat -f %s" % jobName)
	if status != 0:
		return None

	statusList = [s.strip() for s in output.split("\n") if not s.startswith("\t")]
	statusDict = {"job_id": jobId}
	for curStatus in infoList:
		info = curStatus.split("=")
		if len(info) == 2:
			infoDict[info[0].strip().lower()] = info[1].strip()

	return statusDict

STATE_COMPLETED = "C"
STATE_EXITING   = "E"
STATE_HELD      = "H"
STATE_QUEUED    = "Q"
STATE_RUNNING   = "R"
STATE_MOVED     = "M"
STATE_WAITING   = "W"
STATE_SUSPENDED = "S"

def CheckJobCompleted(jobId):
	"""	
	Check if a job is completed. If the job 
	does not exist, it is considered to be completed
	"""
	status = GetJobStatus(jobId)
	if status == None:
		return True

	state = status["job_state"]
	if state == STATE_COMPLETED or state == STATE_EXITING:
		return False

	return True
