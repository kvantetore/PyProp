import sys
import os

import pyprop
if pyprop.IsRunningFromSource:
	sys.path.append(os.path.join(pyprop.BuildPath, "modules",
		"discretizations", "sphericalbasis"))
import libsphericalbasis

from libsphericalbasis import *
from .basis import BasisfunctionSphericalHarmonic
