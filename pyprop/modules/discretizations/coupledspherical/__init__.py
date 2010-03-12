import sys
import os

import pyprop
if pyprop.IsRunningFromSource:
	sys.path.append(os.path.join(pyprop.BuildPath, "modules", "discretizations", "bspline"))
import libcoupledspherical

from libcoupledspherical import *