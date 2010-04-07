import sys
import os

import pyprop
if pyprop.IsRunningFromSource:
    sys.path.append(os.path.join(pyprop.BuildPath, "modules", "potentials", "oscillator"))
import liboscillatorpotential

from liboscillatorpotential import *