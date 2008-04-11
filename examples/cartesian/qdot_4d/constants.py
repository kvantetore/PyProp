from numpy import arccos

#Load unit conversion
sys.path.append(os.path.abspath("./pyprop/pyprop/utilities"))
import units
from units import ElectricFieldAtomicFromIntensitySI as field_from_intensity
from units import AngularFrequencyAtomicFromWavelengthSI as freq_from_wavelength

#Unit conversion factors
femtosec_to_au = 1e-15 / units.constantsAU.time

#Pulse duration from intensity full with half maximum
#for a cos**2 pulse
fwhm_intensity = pi / arccos(0.5**0.25) / 2.0

#Duration of one pulse cycle
time_per_cycle = lambda w : 2 * pi / w


