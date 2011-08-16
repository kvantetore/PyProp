import sys
import os

import pyprop
if pyprop.IsRunningFromSource:
	sys.path.append(os.path.join(pyprop.BuildPath, "modules",
		"discretizations", "finitedifference"))
import libfinitedifference

from libfinitedifference import *
from basis import BasisfunctionFiniteDifference
