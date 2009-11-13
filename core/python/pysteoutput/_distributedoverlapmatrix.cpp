
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <representation/distributedoverlapmatrix.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_python_distributedoverlapmatrix()
{
    class_< DistributedOverlapMatrix<1>, boost::noncopyable >("DistributedOverlapMatrix_1", init<  >())
        .def("SetupRank", &DistributedOverlapMatrix<1>::SetupRank)
        .def("SetupMultivector", &DistributedOverlapMatrix<1>::SetupMultivector)
        .def("MultiVectorToWavefunction", &DistributedOverlapMatrix<1>::MultiVectorToWavefunction)
        .def("MultiplyOverlapRank", &DistributedOverlapMatrix<1>::MultiplyOverlapRank)
        .def("SolveOverlapRank", &DistributedOverlapMatrix<1>::SolveOverlapRank)
    ;

    class_< DistributedOverlapMatrix<2>, boost::noncopyable >("DistributedOverlapMatrix_2", init<  >())
        .def("SetupRank", &DistributedOverlapMatrix<2>::SetupRank)
        .def("SetupMultivector", &DistributedOverlapMatrix<2>::SetupMultivector)
        .def("MultiVectorToWavefunction", &DistributedOverlapMatrix<2>::MultiVectorToWavefunction)
        .def("MultiplyOverlapRank", &DistributedOverlapMatrix<2>::MultiplyOverlapRank)
        .def("SolveOverlapRank", &DistributedOverlapMatrix<2>::SolveOverlapRank)
    ;

    class_< DistributedOverlapMatrix<3>, boost::noncopyable >("DistributedOverlapMatrix_3", init<  >())
        .def("SetupRank", &DistributedOverlapMatrix<3>::SetupRank)
        .def("SetupMultivector", &DistributedOverlapMatrix<3>::SetupMultivector)
        .def("MultiVectorToWavefunction", &DistributedOverlapMatrix<3>::MultiVectorToWavefunction)
        .def("MultiplyOverlapRank", &DistributedOverlapMatrix<3>::MultiplyOverlapRank)
        .def("SolveOverlapRank", &DistributedOverlapMatrix<3>::SolveOverlapRank)
    ;

    class_< DistributedOverlapMatrix<4>, boost::noncopyable >("DistributedOverlapMatrix_4", init<  >())
        .def("SetupRank", &DistributedOverlapMatrix<4>::SetupRank)
        .def("SetupMultivector", &DistributedOverlapMatrix<4>::SetupMultivector)
        .def("MultiVectorToWavefunction", &DistributedOverlapMatrix<4>::MultiVectorToWavefunction)
        .def("MultiplyOverlapRank", &DistributedOverlapMatrix<4>::MultiplyOverlapRank)
        .def("SolveOverlapRank", &DistributedOverlapMatrix<4>::SolveOverlapRank)
    ;

}

