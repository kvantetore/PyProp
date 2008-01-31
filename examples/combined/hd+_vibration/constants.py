
#Physical constants
eV_to_au = 0.036764706
au_to_eV = 27.2114  

inverse_cm_to_au = 4.5563276e-06
au_to_inverse_cm = 219475.

volt_per_metre_to_au = 1.944689e-12
femtosec_to_au = 41.36338

def ReducedMass(M1, M2):
	return M1 * M2 / (M1 + M2)

mass_proton_au = 1836.152666
mass_deuteron_au = 3670.482954
hdp_reduced_mass_au = ReducedMass(mass_proton_au, mass_deuteron_au)

def wavelength_to_frequency(nm):
	return inverse_cm_to_au/(nm*1e-7)


