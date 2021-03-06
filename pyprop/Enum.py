
#Potential Type "enum"
class PotentialType:
	Static = 1
	Dynamic = 3
	FiniteDifference = 4
	CrankNicholson = 5
	Matrix = 6
	RankOne = 7

#Matches core/representation/orthopol/orthopoltools.h PolynomialType
PolynomialType = core.OrthoPolType
	
#Initial condition types
class InitialConditionType:
	Function = 1
	File = 2
	Class = 3
	Custom = 4

class WavefunctionFileFormat:
	Ascii  = 1
	Binary = 2
	HDF    = 3

class IntegratorType:
	IntegratorRK2 = 0
	IntegratorRK4 = 1
	IntegratorRKF45 = 2
	IntegratorRKCK = 3
	IntegratorRK8PD = 4


