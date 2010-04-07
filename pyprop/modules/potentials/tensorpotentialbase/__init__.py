import sys
import os

import pyprop
if pyprop.IsRunningFromSource:
    sys.path.append(os.path.join(pyprop.BuildPath, "modules", "potentials", "tensorpotentialbase"))
import libtensorpotentialbase

from libtensorpotentialbase import *