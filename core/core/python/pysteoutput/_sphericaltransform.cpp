
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <transform/radialtransform.h>
#include <transform/spherical/shtools.h>
#include <transform/spherical/sphericaltransform.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_python_sphericaltransform()
{
    class_< SphericalTransformTensorGrid, boost::noncopyable >("SphericalTransformTensorGrid", no_init)
        .def("GenerateThetaQuadrature", &SphericalTransformTensorGrid::GenerateThetaQuadrature)
        .def("EvaluateAssociatedLegendrePolynomials", &SphericalTransformTensorGrid::EvaluateAssociatedLegendrePolynomials)
        .def("GetLMax", &SphericalTransformTensorGrid::GetLMax)
        .def("GetThetaGrid", &SphericalTransformTensorGrid::GetThetaGrid)
        .def("GetPhiGrid", &SphericalTransformTensorGrid::GetPhiGrid)
        .def("GetOmegaGrid", &SphericalTransformTensorGrid::GetOmegaGrid)
        .def("GetWeights", &SphericalTransformTensorGrid::GetWeights)
        .def("GetAssociatedLegendrePolynomial", &SphericalTransformTensorGrid::GetAssociatedLegendrePolynomial)
        .def("GetSphericalHarmonic", &SphericalTransformTensorGrid::GetSphericalHarmonic)
        .def("GetSphericalHarmonicDerivativeTheta", &SphericalTransformTensorGrid::GetSphericalHarmonicDerivativeTheta)
        .def("GetSphericalHarmonicDerivativePhi", &SphericalTransformTensorGrid::GetSphericalHarmonicDerivativePhi)
        .def("Initialize", &SphericalTransformTensorGrid::Initialize)
        .staticmethod("EvaluateAssociatedLegendrePolynomials")
        .staticmethod("GenerateThetaQuadrature")
    ;

    class_< SphericalTransform<2>, boost::noncopyable >("SphericalTransform_2", init<  >())
        .def_readwrite("transform", &SphericalTransform<2>::transform)
        .def("SetupStep", &SphericalTransform<2>::SetupStep)
        .def("ForwardTransform", &SphericalTransform<2>::ForwardTransform)
        .def("InverseTransform", &SphericalTransform<2>::InverseTransform)
        .def("CreateSphericalHarmonicRepr", &SphericalTransform<2>::CreateSphericalHarmonicRepr)
        .def("CreateAngularRepresentation", &SphericalTransform<2>::CreateAngularRepresentation)
    ;

    class_< SphericalTransform<3>, boost::noncopyable >("SphericalTransform_3", init<  >())
        .def_readwrite("transform", &SphericalTransform<3>::transform)
        .def("SetupStep", &SphericalTransform<3>::SetupStep)
        .def("ForwardTransform", &SphericalTransform<3>::ForwardTransform)
        .def("InverseTransform", &SphericalTransform<3>::InverseTransform)
        .def("CreateSphericalHarmonicRepr", &SphericalTransform<3>::CreateSphericalHarmonicRepr)
        .def("CreateAngularRepresentation", &SphericalTransform<3>::CreateAngularRepresentation)
    ;

    class_< RadialTransform<1>, boost::noncopyable >("RadialTransform_1", init<  >())
        .def("TransformRank", &RadialTransform<1>::TransformRank)
        .def("ForwardTransform", &RadialTransform<1>::ForwardTransform)
        .def("InverseTransform", &RadialTransform<1>::InverseTransform)
    ;

    class_< RadialTransform<2>, boost::noncopyable >("RadialTransform_2", init<  >())
        .def("TransformRank", &RadialTransform<2>::TransformRank)
        .def("ForwardTransform", &RadialTransform<2>::ForwardTransform)
        .def("InverseTransform", &RadialTransform<2>::InverseTransform)
    ;

    class_< RadialTransform<3>, boost::noncopyable >("RadialTransform_3", init<  >())
        .def("TransformRank", &RadialTransform<3>::TransformRank)
        .def("ForwardTransform", &RadialTransform<3>::ForwardTransform)
        .def("InverseTransform", &RadialTransform<3>::InverseTransform)
    ;

}

