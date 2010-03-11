import sys
import os

import pyprop
if pyprop.IsRunningFromSource:
	sys.path.append(os.path.join(pyprop.BuildPath, "modules", "discretizations", "bspline"))
import libbspline

from bspline import InitBSpline
from libbspline import *