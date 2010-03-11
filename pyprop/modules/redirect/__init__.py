import sys
import os.path

import pyprop
if pyprop.IsRunningFromSource:
	sys.path.append(os.path.join(pyprop.BuildPath, "modules", "redirect"))
import libredirect


from redirect import Redirect

