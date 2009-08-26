import os
from numpy import arccos

#Load unit conversion
sys.path.append(os.path.abspath("./pyprop/pyprop/utilities"))
import units
from units import ElectricFieldAtomicFromIntensitySI as field_from_intensity
from units import AngularFrequencyAtomicFromWavelengthSI as freq_from_wavelength
from units import IntensitySIFromElectricFieldAtomic as intensity_from_field

#Unit conversion factors
femtosec_to_au = 1e-15 / units.constantsAU.time

#Pulse duration from intensity full with half maximum
#for a cos**2 pulse
fwhm_intensity = pi / arccos(0.5**0.25) / 2.0

#Converts wavelength in nm -> time of one cycle a.u.
cycletime_from_wavelength = lambda l: 2 * pi / freq_from_wavelength(l)

#Converts frequency in a.u. -> time of one cycle a.u.
cycletime_from_frequency = lambda f: 2 * pi / f

#Ponderomotive energy
#ponderomotive_energy = lambda I, omega: I / (4.0 * omega**2)
ponderomotive_energy = lambda E0, omega: E0**2 / (4.0 * omega**2)

#eV -> au
eV_to_au = 3.674932540e-2

def ponderomotive_energy_au_from_SI(I, wavelength):
	"""
		I: intensity in W/cm**2
		wavelength: wavelength in nm
		Returns ponderomotive energy in a.u.
	"""

	#I_au = units.SItoAtomic(I, units.UNIT_FIELD_INTENSITY)
	omega_au = freq_from_wavelength(wavelength)
	field_au = field_from_intensity(I)
	Up = ponderomotive_energy(field_au, omega_au)
	return Up

def keldysh_parameter(I, wavelength, I0):
	"""
	I: Laser intensity [W/cm**2]
	wavelength: laser wavelength [nm]
	I0: Ionization potential [a.u.]
	"""
	#I_au = units.SItoAtomic(I, units.UNIT_FIELD_INTENSITY)
	omega_au = freq_from_wavelength(wavelength)
	#gamma = 2 * omega_au**2 * I0 / I_au
	Up = ponderomotive_energy(field_from_intensity(I), freq_from_wavelength(wavelength))
	gamma = sqrt(I0 / (2 * Up))

	return gamma

def critical_ionization_intensity(Ip):
	return 4e9 * (Ip / eV_to_au)**4 / 1e15

def tpdi_crossection(freq, intensity, pIon, pulseDuration):
	tau = 35.0 / 128.0 * pulseDuration * units.constantsAU.time
	photonEnergy = units.constantsAU.energy * freq
	sigma = (photonEnergy / intensity)**2 / tau * pIon
	print
	print "Photon energy = %.1f eV" % (freq / eV_to_au)
	print "TPDI cross section = %1.2e cm**4 s" % sigma
	print
	return freq / eV_to_au, sigma
	

class Levels():
	_1s = {"n": 0, "l": 0, "name": "2p"}
	_2s = {"n": 1, "l": 0, "name": "2s"}
	_2p = {"n": 1, "l": 1, "name": "2p"}
	_3s = {"n": 2, "l": 0, "name": "3s"}
	_3p = {"n": 2, "l": 1, "name": "3p"}
	_3d = {"n": 2, "l": 2, "name": "3d"}
	_4s = {"n": 3, "l": 0, "name": "4s"}
	_4p = {"n": 3, "l": 1, "name": "4p"}
	_4d = {"n": 3, "l": 2, "name": "4d"}
	_4f = {"n": 3, "l": 3, "name": "4f"}
	_5s = {"n": 4, "l": 0, "name": "5s"}
	_5p = {"n": 4, "l": 1, "name": "5p"}
	_5d = {"n": 4, "l": 2, "name": "5d"}
	_5f = {"n": 4, "l": 3, "name": "5f"}
	_5g = {"n": 4, "l": 4, "name": "5g"}
