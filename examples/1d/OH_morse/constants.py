
#Physical constants
eV_to_au = 0.036764706
au_to_eV = 27.2114  
Ry_to_au = 0.5

inverse_cm_to_au = 4.5563276e-06
au_to_inverse_cm = 219475.

volt_per_metre_to_au = 1.944689e-12
femtosec_to_au = 41.36338

def ReducedMass(M1, M2):
	return M1 * M2 / (M1 + M2)

mass_hydrogen_au = 1837.3623595944798
mass_oxygen_au = 29166.217982728809
oh_reduced_mass_au = ReducedMass(mass_hydrogen_au, mass_oxygen_au)

def wavelength_to_frequency(nm):
	return inverse_cm_to_au/(nm*1e-7)


