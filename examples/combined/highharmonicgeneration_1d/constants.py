import os
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

#Converts wavelength in nm -> time of one cycle a.u.
cycletime_from_wavelength = lambda l: 2 * pi / freq_from_wavelength(l)

#Ponderomotive energy
ponderomotive_energy = lambda I, omega: I / (4.0 * omega**2)

#eV -> au
eV_to_au = 3.674932540e-2

def ponderomotive_energy_au_from_SI(I, wavelength):
	"""
		I: intensity in W/cm**2
		wavelength: wavelength in nm
		Returns ponderomotive energy in a.u.
	"""

	I_au = units.SItoAtomic(I, units.UNIT_FIELD_INTENSITY)
	omega_au = freq_from_wavelength(wavelength)
	Up = ponderomotive_energy(I_au, omega_au)
	
	return Up

def keldysh_parameter(I, wavelength, I0):
	"""
	I: Laser intensity [W/cm**2]
	wavelength: laser wavelength [nm]
	I0: Ionization potential [a.u.]
	"""
	I_au = units.SItoAtomic(I, units.UNIT_FIELD_INTENSITY)
	omega_au = freq_from_wavelength(wavelength)
	gamma = 2 * omega_au**2 * I0 / I_au

	return gamma
