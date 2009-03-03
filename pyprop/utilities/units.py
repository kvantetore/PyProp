from math import pi
from math import sqrt

"""
Unit conversion to and from atomic units

Original version, Raymond Nepstad, 24. Jan. 07
"""

UNIT_ANGSTROM = 'AA'
UNIT_CHARGE = 'C'
UNIT_ENERGY = 'J'
UNIT_FIELD_INTENSITY = 'W/cm2'
UNIT_FIELD_STRENGTH = "V/cm"
UNIT_FREQUENCY = "Hz"
UNIT_MASS = 'kg'
UNIT_NANOMETER = 'nm'
UNIT_TIME = 's'

def nice(toBeNiced):
  return '%.17e' % toBeNiced

#
# Unit system conversion
#
def AtomicToSI(number,what):
  "Convert from various atomic units to SI"
  if(what == UNIT_FIELD_STRENGTH):
    return number*(constantsAU.electric_field_strength/100)
  elif(what == UNIT_FREQUENCY):
    return number/(constantsAU.time)
  elif(what == UNIT_MASS):
    return number*constansAU.mass
  else:
    raise "Undefined type", what

def SItoAtomic(number,what):
  "Convert from various SI units to atomic units"
  if(what == UNIT_FIELD_INTENSITY):
    return number/(constantsAU.intensity/100.0**2)
  elif(what == UNIT_FIELD_STRENGTH):
    return number / (constantsAU.electric_field_strength/100.) 
  elif(what == UNIT_FREQUENCY):
    return number * constantsAU.time 
  else:
    raise "Undefined type", what

#
# Unit type conversion
#
def AngularFrequencyAtomicFromWavelengthSI(wavelength):
  """
  Wavelength [nm] -> angular freq. [a.u]
  """
  return SItoAtomic(2*pi*constantsSI.lightSpeed/(wavelength*1e-9),UNIT_FREQUENCY)

def ElectricFieldAtomicFromIntensitySI(intensity):
  """
  Intensity [W/cm**2] -> E-field strength [a.u.]

  Relation obtained from time-averaging over one cycle,

  E0 = sqrt( (2 <I>) / (eps0 * c) )
  """
  #return SItoAtomic(sqrt(intensity/0.0013),UNIT_FIELD_STRENGTH)
  return SItoAtomic(sqrt(2 * intensity / (constantsAU.electrostatic_constant / (4*pi) * constantsSI.lightSpeed)), UNIT_FIELD_STRENGTH)

def WavelengthSIFromAngularFrequencyAtomic(freq):
  """
  Angular frequency [a.u] -> Wavelength [nm]
  """
  return 1e9*2*pi*constantsSI.lightSpeed/(AtomicToSI(freq,UNIT_FREQUENCY))

def IntensitySIFromElectricFieldAtomic(efield):
  """
  Electric field [a.u.] -> Intensity [W/cm**2]
  """
  #return 0.0013*(AtomicToSI(efield,UNIT_FIELD_STRENGTH))**2 
  scaling = (constantsAU.electrostatic_constant / (4*pi) * constantsSI.lightSpeed / 2.0)
  return scaling * AtomicToSI(efield,UNIT_FIELD_STRENGTH)**2



def EnergyElectronVoltFromAU(energy):
	"""
	Energy [J] -> Energy [eV]
	"""
	return constantsAU.energy * energy / constantsAU.charge
#
# Constants
#
class constantsAU:
  mass = 9.10938e-31                    #  -> kg
  charge = 1.60218e-19                  #  -> C
  hbar = 1.05457e-34                    #  -> Js
  electrostatic_constant = 1.11265e-10  #  -> 1/Fm
  alpha = 1./137.03599911
  
  # Derived units
  lightSpeed = 1./alpha
  length = electrostatic_constant*hbar**2/(mass*charge**2)  # AA -> m
  energy = hbar**2/(mass*length**2)                         # hartree -> J
  time = hbar/energy                                        #  -> s
  velocity = length/time                                    #  -> m/s
  angular_frequency = velocity/(2*pi*length)                #  -> Hz
  electric_field_strength = charge/(electrostatic_constant*length**2) #  -> V/m
  intensity = energy**2 / (hbar * length**2)                #  -> W/m**2

class constantsSI:
  lightSpeed = 299792458
  alpha = 1./137.03599911

