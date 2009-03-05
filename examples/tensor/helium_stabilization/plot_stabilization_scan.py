
import tables
from pylab import *
from numpy import *

f = tables.openFile("stabilization_scan.h5")
I = f.root.intensity

plot(I, f.root.totalIon, "-", label="Total Ionization")
plot(I, f.root.singleIon, "--", label="Single Ionization" )
plot(I, f.root.doubleIon, ":", label="Double Ionization")
legend(loc="lower right")
xlabel("Intensity")

