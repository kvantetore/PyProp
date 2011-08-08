import os
import os.path as path

#Find out whether we running from a source or installed pyprop
IsRunningFromSource = not path.exists(path.join(__path__[0], "core", "libcore.so"))

if IsRunningFromSource:
	#find out which build kind we are
	BuildKind = os.environ.get("PYPROP_BUILD_KIND", "default")
	BuildPath = path.abspath(path.join(__path__[0], "..", "build", BuildKind, "pyprop"))

import core
import config
from distribution import ProcId, ProcCount
