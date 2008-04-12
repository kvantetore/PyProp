import os
import ConfigParser

class Section:
	
	def __init__(self, name, cfg=None):
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
		except:
			print "Locals:", self.__dict__
			print "Option:", optionString
			raise 
				
		self.Set(optionName, optionValue)
		
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

def Load(fileName, silent=True):
	#Find all imported files in the config files
	#hierarchy
	configFiles = []
	newFiles = [(".", fileName)]
	while len(newFiles) > 0:
		curDir, curFile = newFiles.pop()
		localFile = curDir.rstrip("/") + "/" + curFile
		absFile = os.path.abspath(os.path.expanduser(os.path.expandvars(localFile)))
		absDir = os.path.dirname(absFile)
		if not silent:
			print "Using config file", absFile
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
