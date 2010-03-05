
#Potential Type "enum"
##Matches core/representation/orthopol/orthopoltools.h PolynomialType
#if core.LOAD_CORE_OK:
#	PolynomialType = core.OrthoPolType
	
#Initial condition types
class IntegratorType:
	IntegratorRK2 = 0
	IntegratorRK4 = 1
	IntegratorRKF45 = 2
	IntegratorRKCK = 3
	IntegratorRK8PD = 4


