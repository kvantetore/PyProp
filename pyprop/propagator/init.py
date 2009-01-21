import sys

execfile(__path__[0] + "/propagator/Propagator.py")
execfile(__path__[0] + "/propagator/CartesianPropagator.py")
execfile(__path__[0] + "/propagator/KrylovPropagator.py")
execfile(__path__[0] + "/propagator/PamPropagator.py")
execfile(__path__[0] + "/propagator/ArpackSolver.py")
execfile(__path__[0] + "/propagator/PiramSolver.py")
execfile(__path__[0] + "/propagator/CartesianMixedPropagator.py")
execfile(__path__[0] + "/propagator/CombinedPropagator.py")
execfile(__path__[0] + "/propagator/TransformedGridPropagator.py")
execfile(__path__[0] + "/propagator/ExponentialFiniteDifferencePropagator.py")
execfile(__path__[0] + "/propagator/OdePropagator.py")
execfile(__path__[0] + "/propagator/RungeKuttaPropagator.py")
execfile(__path__[0] + "/propagator/VectorPropagator.py")
execfile(__path__[0] + "/propagator/CayleyPropagator.py")
#execfile(__path__[0] + "/propagator/Krotov.py")

#init subpropagators
execfile(__path__[0] + "/propagator/subpropagator/init.py")

#init optimal control solvers
#try:
#execfile(__path__[0] + "/propagator/optimalcontrol/init.py")
oct_module_path = os.environ["OCT_PYPROP_MODULE_PATH"]
execfile(oct_module_path + "/init.py")
execfile(__path__[0] + "/propagator/Krotov.py")
#except:
#	print "Oops, could not load optimal control solvers!"
