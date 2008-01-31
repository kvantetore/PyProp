
def RemoveNodeIfExists(f, path, node):
	"""
	Removes a node from a pytables file f if it exists.

	f: an open pytables file
	path: a group or a path to a group where the node is
	node: the name of the node in the group specified by path
	"""
	try:
		f.removeNode(path, node)
	except tables.NoSuchNodeError: pass

def SaveArray(f, path, node, data):
	RemoveNodeIfExists(f, path, node)
	tab = f.createArray(path, node, data, createparents=True)
	tab.close()

