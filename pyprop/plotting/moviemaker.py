import os
import sys
import pylab
import tables
import numpy

"""
Example config sections:

[Movie]
encoder = "ffmpeg"
framerate = 25
total_frames = 200
frame_dpi = 100
frame_size = 600
bitrate = 8000
tmpdir = "movie"
movie_name = "output.avi"

[Rebuilder]
rank0 = [-10, 10, 100]
rank1 = [-10, 10, 100]
lmax = 1
"""

class MovieMaker:

	def __init__(self, wavefunctionRebuilder=None):

		self.Rebuilder = wavefunctionRebuilder
		if wavefunctionRebuilder != None:
			self.Propagator = self.Rebuilder.Propagator

	def ApplyConfigSection(self, configSection):
		self.Encoder = configSection.encoder
		self.FrameRate = configSection.framerate
		self.TotalFrames = configSection.total_frames
		self.FrameDPI = configSection.frame_dpi
		self.FrameSize = configSection.frame_size
		self.BitRate = configSection.bitrate
		self.TmpDir = configSection.tmpdir
		self.MovieName = configSection.movie_name
		self.InputFileName = None
		if hasattr(configSection, "input_name"):
			self.InputFileName = configSection.input_name

	def PlotFrame(self, curIter, t):
		X, Y, psi = self.Rebuilder.BuildWavefunction()
		pylab.imshow(numpy.abs(psi)**2, extent=(self.Rebuilder.GridExtent), interpolation="gaussian", vmin=0.0)


	def CreateFrames(self, plotFunction = PlotFrame):
		"""
		Propagate the wavefunction, storing movie frames for the duration.
		"""

		#Plot non-interactively
		interactive = pylab.rcParams['interactive']
		pylab.rcParams['interactive'] = False
		figSizeInch = self.FrameSize / float(self.FrameDPI)
		pylab.figure(figsize = (figSizeInch, figSizeInch), dpi = self.FrameDPI)

		#Create tmp directory to hold movie and image frames
		try:
			os.mkdir(self.TmpDir, 0755)
		except:
			print "Could not make tmp directory (%s)" % self.TmpDir
			raise Exception

		#Create movie frames
		print "Creating movie frames "
		for i, t in enumerate(self.Propagator.Advance(self.TotalFrames)):
			#Write progress info
			sys.stdout.write(10 * "\b")
			sys.stdout.write("%3i / %3i" % (i, self.TotalFrames))
			sys.stdout.flush()
			#X, Y, psi = self.Rebuilder.BuildWavefunction()
			#pylab.imshow(abs(psi)**2, extent=(self.Rebuilder.GridExtent))
			plotFunction(i, t)
			pylab.savefig("%s/frame%03i.png" % (self.TmpDir, i), dpi = self.FrameDPI)
			pylab.clf()

		print "Done creating movie frames!"
		pylab.rcParams['interactive'] = interactive


	def PlotFramesFromHDF5(self, **args):
		"""
		Make movie frames from wavefunctions stored in HDF5 files.
		"""

		wavefunctionFile = args["wavefunctionFile"]
		wavefunctionDataset = args["wavefunctionDataset"]

		sys.stdout.write("Creating frames ")

		#Create tmp directory to hold movie and image frames
		if not os.path.isdir(self.TmpDir):
			try:
				os.mkdir(self.TmpDir, 0755)
			except:
				print "Could not make tmp directory (%s)" % self.TmpDir
				raise Exception

		#Plot non-interactively
		interactive = pylab.rcParams['interactive']
		pylab.rcParams['interactive'] = False
		figSizeInch = self.FrameSize / float(self.FrameDPI)
		pylab.figure(figsize = (figSizeInch, figSizeInch), dpi = self.FrameDPI)
	

		psi = self.Rebuilder.Propagator.psi
		sys.stdout.write("%04i/%04i" % (0, self.TotalFrames)) 
		for i in range(self.TotalFrames):
			
			#Print progress info
			infoString = "\b"*9 + "%04i/%04i" % (i, self.TotalFrames)
			sys.stdout.write(infoString)
			sys.stdout.flush()

			#Load wavefunction
			h5file = tables.openFile("%s%04i.h5" % (wavefunctionFile, i), "r")
			psi.GetData()[:] = h5file.getNode(wavefunctionDataset)[:]
			h5file.close()

			#Create movie frame
			self.PlotFrame(i, 0)
			pylab.savefig("%s/frame%03i.png" % (self.TmpDir, i), dpi = self.FrameDPI)
			pylab.clf()

		print "Done creating movie frames!"
		pylab.rcParams['interactive'] = interactive


	def CreateMovie(self):
		"""
		Make movie from pre-made movie frames.
		"""
		print "Now creating movie..."
		encoderCmd = self.Encoder
		if self.Encoder == "ffmpeg":
			if self.FrameRate < 20:
				print "WARNING: ffmpeg with mpeg-codec does not support fps < 20!"
	
			framerate = max(self.FrameRate, 20)
			if self.InputFileName != None:
				encoderCmd += " -i %s" % self.InputFileName
			else:
				encoderCmd += " -i %s%s " % (self.TmpDir, "/frame%3d.png")
			encoderCmd += " -r %i" % framerate
			encoderCmd += " -b %ik -bufsize %ik -f mpeg2video -s %ix%i" % (self.BitRate, self.BitRate, self.FrameSize, self.FrameSize)
			encoderCmd += "  %s" % self.MovieName
	
		elif self.Encoder == "mencoder":
			encoderCmd += " -ovc lavc -fps %i -lavcopts" % (self.FrameRate)
			encoderCmd += " vcodec=mpeg4:mbd=2:cbp:trell:vbitrate=%i:autoaspect" % self.BitRate
			encoderCmd += " -vf scale=%i:%i -ffourcc DX50" % (self.FrameSize, self.FrameSize)
			encoderCmd += " mf://%s/*png" % self.TmpDir
			encoderCmd += " -o %s" % self.MovieName

		else:
			raise Exception("Unknown encoder '%s'" % encoderCmd)

		print "Running Encoder: %s" % encoderCmd
		os.system(encoderCmd)

