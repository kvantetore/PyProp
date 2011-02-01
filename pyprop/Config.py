import os
import ConfigParser
import pyproplogging

class Section:

	def __init__(self, name, cfg=None):
		self.Logger = pyproplogging.GetClassLogger(self)
		self.name = name
		if cfg != None:
			self.LoadConfig(name, cfg)

	def LoadConfig(self, name, cfg):
		for optionName in cfg.options(name):
			if optionName == "base":
				importName = eval(cfg.get(name, optionName))
				self.LoadConfig(importName, cfg)
			else:
				self.SetOptionString(optionName, cfg.get(name, optionName))
			
	def SetOptionString(self, optionName, optionString):
		glob = dict(ProjectNamespace)
		glob.update(globals())
		try:
			optionValue = eval(optionString, glob, self.__dict__)
			self.Set(optionName, optionValue)
		except:
			self.Logger.error("Missing from namespace: %s" % optionString)
			self.Logger.debug("Locals: %s" % self.__dict__)
			#raise 
				
	def Set(self, optionName, optionValue):
		if optionName in ["name", "SetOptionString", "Set", "Get"]:
			raise Exception, "Invalid optionName %s" % optionName
		self.__dict__[optionName] = optionValue
	
	def Get(self, optionName):
		return self.__dict__[optionName]

	def Exists(self, optionName):
		return hasattr(self, optionName)
		
	def Apply(self, other):
		if hasattr(other,"ApplyConfigSection"):
			#print "Applying config section to ", other.__class__.__name__
			other.ApplyConfigSection(self)

class Config:

	def __init__(self, cfg):	
		self.Logger = pyproplogging.GetClassLogger(self)
		for sectionName in cfg.sections():
			if sectionName in ["GetSection"]:
				raise Exception, "Invalid sectionName %s" % sectionName
			section = Section(sectionName, cfg)
			section.Config = self
			self.__dict__[sectionName] = section
			self.cfgObj = cfg

	def GetSection(self, sectionName):
		return self.__dict__[sectionName]
		
	def Apply(self, other):
		if hasattr(other, "ApplyConfig"):
			#print "Applying config to ", other.__class__.__name__
			other.ApplyConfig(self)

	def SetValue(self, section, optionName, value, infoLevel=0, ignoreSerializationError=False):
		"""
		Update 'optionName' in 'section' with 'value'. Both the Config
		object and the original ConfigParser object is updated.
		If an un-repr-able is passed, exception is raised unless 'ignoreSerializationError'
		is set to True.
		"""
		try:
			eval(repr(value))
		except:
			if not ignoreSerializationError:
				raise Exception("Object cannot be serialized: %s", repr(value))

		self.cfgObj.set(section, optionName, repr(value))
		self.__dict__[section].Set(optionName, value)

		self.Logger.debug("Changing %s.%s to %s" % (section, optionName, value))

def Load(fileName, silent=True):
	#Find all imported files in the config files
	#hierarchy
	if fileName is None:
		raise Exception("Filename can not be none")
	configFiles = []
	newFiles = [(".", fileName)]
	while len(newFiles) > 0:
		curDir, curFile = newFiles.pop()
		localFile = curDir.rstrip("/") + "/" + curFile
		absFile = os.path.abspath(os.path.expanduser(os.path.expandvars(localFile)))
		absDir = os.path.dirname(absFile)
		if not silent:
			PrintOut("Using config file %s" % absFile)
		configFiles.insert(0, absFile)
		
		#find Import statements from curFile
		cfg = ConfigParser.ConfigParser()
		cfg.readfp(open(absFile))
		if cfg.has_section("Import"):
			if cfg.has_option("Import", "files"):
				for f in eval(cfg.get("Import", "files")):
					newFiles.insert(0, (absDir, f))
		del cfg

	#Create a fresh config parser object
	cfg = ConfigParser.ConfigParser()
	#read all config files
	for file in configFiles:
		cfg.readfp(open(file))
	
	#parse the configparser object
	parsedConfig = Config(cfg)

	return parsedConfig


def UpdateConfig(conf, updateParams, addMissingSections=False):
	"""
	Update Pyprop config object with values given in updateParams.

	Input
	-----
	conf: pyprop config object
	updateParams: list of tuples, [('section', 'param', 'value'), ...]
	
	
	Returns
	-------
	Updated pyprop config object
	
	
	Note: Section references (i.e. through 'base') will be updated as well.

	"""
	tmpConf = Config(conf.cfgObj)
	logger = GetFunctionLogger()
	
	#Update config
	for section, param, val in updateParams:
		if not hasattr(tmpConf, section) and addMissingSections:
			logging.info("Config object did not contain section: %s, creating it now." % section)
			cfg = tmpConf.cfgObj
			cfg._dict = dict
			cfg.add_section(section)
			tmpConf = Config(cfg)
		logger.debug("Updating config: %s(%s): %s" % (section, param, val))
		tmpConf.SetValue(section, param, val)

	#Update config object from possible changed ConfigParser object
	newConf = Config(tmpConf.cfgObj)
	
	return newConf
