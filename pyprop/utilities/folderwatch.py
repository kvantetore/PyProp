import sys
import os
import os.path

CHANGE_NEW = 1
CHANGE_DELETED = 2
CHANGE_MODIFIED = 3

class FolderWatch(object):
	
	def __init__(self, folder, changeThreshold=60):
		self.Folder = os.path.abspath(folder)
		self.CurrentFiles = []
		self.ChangeThreshold = changeThreshold
	
	def GetUpdatedFiles(self):
		files = filter(self.__FileFilter, os.listdir(self.Folder))
		updatedFiles = [(f, self.__LastModified(f)) for f in files]
		return updatedFiles
	
	def GetChanges(self):
		currentFiles = list(self.CurrentFiles)
		updatedFiles = self.GetUpdatedFiles()
		newCurrentFiles = []
		changes = []

		#Loop through all updated files, and see whether they are 
		#new or modified. 
		for filename, lastModified in updatedFiles:
			change = None
			existingFile = self.__CheckExistingFile(filename)

			#if modification date is newer then threshold seconds, 
			#ignore (unless it is delete) the change until it is older
			ignoreChanges = (time.time() - lastModified) < self.ChangeThreshold

			#If the file is in the current list, remove it
			#so we can figure out which files are missing
			if existingFile:
				currentFiles.remove(existingFile) 

			#check if file is new, changed, or not modified
			if existingFile:
				f, prevLastModified = existingFile
				if ignoreChanges:
					lastModified = prevLastModified
				else:
					if lastModified > prevLastModified:
						change = CHANGE_MODIFIED
			else:
				if ignoreChanges:
					filename = None
				else:
					change = CHANGE_NEW

			#Update current files
			if filename != None:
				newCurrentFiles.append((filename, lastModified))

			#record change
			if change:
				changes.append( (filename, change) )

		#All items left in current files are missing, i.e. deleted
		for filename, lastModified in currentFiles:
			changes.append( (filename, CHANGE_DELETED) )

		#replace the list of current files
		self.CurrentFiles = newCurrentFiles

		return changes

	def __FileFilter(self, filename):
		if not self.__IsFile(filename):
			return False
		return True
	
	def __IsFile(self, filename):
		return stat.S_ISREG(os.stat(filename).st_mode)
	
	def __LastModified(self, filename):
		return os.stat(filename)[stat.ST_MTIME]

	def __CheckExistingFile(self, filename):
		fileList = self.CurrentFiles
		for f, modified in fileList:
			if f == filename:
				return (f, modified)
		return None


def TestFolderWatch():
	watch = FolderWatch(".", changeThreshold=3)
	while True:
		#look for changes
		changes = watch.GetChanges()
		for filename, change in changes:
			if change == CHANGE_NEW:
				print "New file %s" % (filename)
			elif change == CHANGE_MODIFIED:
				print "Modified file %s" % (filename)
			elif change == CHANGE_DELETED:
				print "Deleted file %s" % (filename)
			else:
				print "what?"

		#sleep a bit
		sys.stdout.flush()
		time.sleep(0.5)

