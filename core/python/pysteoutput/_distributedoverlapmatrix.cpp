
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <representation/distributedoverlapmatrix.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

void (DistributedOverlapMatrix<1>::*DistributedOverlapMatrix_1___MultiplyOverlapRankWavefunction_1___Wavefunction_1___int_bool)(Wavefunction<1>&, Wavefunction<1>&, int, bool)  = &DistributedOverlapMatrix<1>::MultiplyOverlapRank;

void (DistributedOverlapMatrix<1>::*DistributedOverlapMatrix_1___MultiplyOverlapRankWavefunction_1___int_bool)(Wavefunction<1>&, int, bool)  = &DistributedOverlapMatrix<1>::MultiplyOverlapRank;

void (DistributedOverlapMatrix<1>::*DistributedOverlapMatrix_1___SolveOverlapRankWavefunction_1___Wavefunction_1___int_bool)(Wavefunction<1>&, Wavefunction<1>&, int, bool)  = &DistributedOverlapMatrix<1>::SolveOverlapRank;

void (DistributedOverlapMatrix<1>::*DistributedOverlapMatrix_1___SolveOverlapRankWavefunction_1___int_bool)(Wavefunction<1>&, int, bool)  = &DistributedOverlapMatrix<1>::SolveOverlapRank;

void (DistributedOverlapMatrix<2>::*DistributedOverlapMatrix_2___MultiplyOverlapRankWavefunction_2___Wavefunction_2___int_bool)(Wavefunction<2>&, Wavefunction<2>&, int, bool)  = &DistributedOverlapMatrix<2>::MultiplyOverlapRank;

void (DistributedOverlapMatrix<2>::*DistributedOverlapMatrix_2___MultiplyOverlapRankWavefunction_2___int_bool)(Wavefunction<2>&, int, bool)  = &DistributedOverlapMatrix<2>::MultiplyOverlapRank;

void (DistributedOverlapMatrix<2>::*DistributedOverlapMatrix_2___SolveOverlapRankWavefunction_2___Wavefunction_2___int_bool)(Wavefunction<2>&, Wavefunction<2>&, int, bool)  = &DistributedOverlapMatrix<2>::SolveOverlapRank;

void (DistributedOverlapMatrix<2>::*DistributedOverlapMatrix_2___SolveOverlapRankWavefunction_2___int_bool)(Wavefunction<2>&, int, bool)  = &DistributedOverlapMatrix<2>::SolveOverlapRank;

void (DistributedOverlapMatrix<3>::*DistributedOverlapMatrix_3___MultiplyOverlapRankWavefunction_3___Wavefunction_3___int_bool)(Wavefunction<3>&, Wavefunction<3>&, int, bool)  = &DistributedOverlapMatrix<3>::MultiplyOverlapRank;

void (DistributedOverlapMatrix<3>::*DistributedOverlapMatrix_3___MultiplyOverlapRankWavefunction_3___int_bool)(Wavefunction<3>&, int, bool)  = &DistributedOverlapMatrix<3>::MultiplyOverlapRank;

void (DistributedOverlapMatrix<3>::*DistributedOverlapMatrix_3___SolveOverlapRankWavefunction_3___Wavefunction_3___int_bool)(Wavefunction<3>&, Wavefunction<3>&, int, bool)  = &DistributedOverlapMatrix<3>::SolveOverlapRank;

void (DistributedOverlapMatrix<3>::*DistributedOverlapMatrix_3___SolveOverlapRankWavefunction_3___int_bool)(Wavefunction<3>&, int, bool)  = &DistributedOverlapMatrix<3>::SolveOverlapRank;

void (DistributedOverlapMatrix<4>::*DistributedOverlapMatrix_4___MultiplyOverlapRankWavefunction_4___Wavefunction_4___int_bool)(Wavefunction<4>&, Wavefunction<4>&, int, bool)  = &DistributedOverlapMatrix<4>::MultiplyOverlapRank;

void (DistributedOverlapMatrix<4>::*DistributedOverlapMatrix_4___MultiplyOverlapRankWavefunction_4___int_bool)(Wavefunction<4>&, int, bool)  = &DistributedOverlapMatrix<4>::MultiplyOverlapRank;

void (DistributedOverlapMatrix<4>::*DistributedOverlapMatrix_4___SolveOverlapRankWavefunction_4___Wavefunction_4___int_bool)(Wavefunction<4>&, Wavefunction<4>&, int, bool)  = &DistributedOverlapMatrix<4>::SolveOverlapRank;

void (DistributedOverlapMatrix<4>::*DistributedOverlapMatrix_4___SolveOverlapRankWavefunction_4___int_bool)(Wavefunction<4>&, int, bool)  = &DistributedOverlapMatrix<4>::SolveOverlapRank;


}// namespace 


// Module ======================================================================
void Export_python_distributedoverlapmatrix()
{
    class_< DistributedOverlapMatrix<1>, boost::noncopyable >("DistributedOverlapMatrix_1", init<  >())
        .def("SetupRank", &DistributedOverlapMatrix<1>::SetupRank)
        .def("SetupMultivector", &DistributedOverlapMatrix<1>::SetupMultivector)
        .def("MultiVectorToWavefunction", &DistributedOverlapMatrix<1>::MultiVectorToWavefunction)
        .def("WavefunctionToMultiVector", &DistributedOverlapMatrix<1>::WavefunctionToMultiVector)
        .def("MultiplyOverlapRank", DistributedOverlapMatrix_1___MultiplyOverlapRankWavefunction_1___Wavefunction_1___int_bool)
        .def("MultiplyOverlapRank", DistributedOverlapMatrix_1___MultiplyOverlapRankWavefunction_1___int_bool)
        .def("SolveOverlapRank", DistributedOverlapMatrix_1___SolveOverlapRankWavefunction_1___Wavefunction_1___int_bool)
        .def("SolveOverlapRank", DistributedOverlapMatrix_1___SolveOverlapRankWavefunction_1___int_bool)
    ;

    class_< DistributedOverlapMatrix<2>, boost::noncopyable >("DistributedOverlapMatrix_2", init<  >())
        .def("SetupRank", &DistributedOverlapMatrix<2>::SetupRank)
        .def("SetupMultivector", &DistributedOverlapMatrix<2>::SetupMultivector)
        .def("MultiVectorToWavefunction", &DistributedOverlapMatrix<2>::MultiVectorToWavefunction)
        .def("WavefunctionToMultiVector", &DistributedOverlapMatrix<2>::WavefunctionToMultiVector)
        .def("MultiplyOverlapRank", DistributedOverlapMatrix_2___MultiplyOverlapRankWavefunction_2___Wavefunction_2___int_bool)
        .def("MultiplyOverlapRank", DistributedOverlapMatrix_2___MultiplyOverlapRankWavefunction_2___int_bool)
        .def("SolveOverlapRank", DistributedOverlapMatrix_2___SolveOverlapRankWavefunction_2___Wavefunction_2___int_bool)
        .def("SolveOverlapRank", DistributedOverlapMatrix_2___SolveOverlapRankWavefunction_2___int_bool)
    ;

    class_< DistributedOverlapMatrix<3>, boost::noncopyable >("DistributedOverlapMatrix_3", init<  >())
        .def("SetupRank", &DistributedOverlapMatrix<3>::SetupRank)
        .def("SetupMultivector", &DistributedOverlapMatrix<3>::SetupMultivector)
        .def("MultiVectorToWavefunction", &DistributedOverlapMatrix<3>::MultiVectorToWavefunction)
        .def("WavefunctionToMultiVector", &DistributedOverlapMatrix<3>::WavefunctionToMultiVector)
        .def("MultiplyOverlapRank", DistributedOverlapMatrix_3___MultiplyOverlapRankWavefunction_3___Wavefunction_3___int_bool)
        .def("MultiplyOverlapRank", DistributedOverlapMatrix_3___MultiplyOverlapRankWavefunction_3___int_bool)
        .def("SolveOverlapRank", DistributedOverlapMatrix_3___SolveOverlapRankWavefunction_3___Wavefunction_3___int_bool)
        .def("SolveOverlapRank", DistributedOverlapMatrix_3___SolveOverlapRankWavefunction_3___int_bool)
    ;

    class_< DistributedOverlapMatrix<4>, boost::noncopyable >("DistributedOverlapMatrix_4", init<  >())
        .def("SetupRank", &DistributedOverlapMatrix<4>::SetupRank)
        .def("SetupMultivector", &DistributedOverlapMatrix<4>::SetupMultivector)
        .def("MultiVectorToWavefunction", &DistributedOverlapMatrix<4>::MultiVectorToWavefunction)
        .def("WavefunctionToMultiVector", &DistributedOverlapMatrix<4>::WavefunctionToMultiVector)
        .def("MultiplyOverlapRank", DistributedOverlapMatrix_4___MultiplyOverlapRankWavefunction_4___Wavefunction_4___int_bool)
        .def("MultiplyOverlapRank", DistributedOverlapMatrix_4___MultiplyOverlapRankWavefunction_4___int_bool)
        .def("SolveOverlapRank", DistributedOverlapMatrix_4___SolveOverlapRankWavefunction_4___Wavefunction_4___int_bool)
        .def("SolveOverlapRank", DistributedOverlapMatrix_4___SolveOverlapRankWavefunction_4___int_bool)
    ;

}

