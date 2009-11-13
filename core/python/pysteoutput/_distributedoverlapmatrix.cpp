
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
        .def("MultiplyOverlapRank", &DistributedOverlapMatrix<1>::MultiplyOverlapRank)
    ;

    class_< DistributedOverlapMatrix<2>, boost::noncopyable >("DistributedOverlapMatrix_2", init<  >())
        .def("MultiplyOverlapRank", &DistributedOverlapMatrix<2>::MultiplyOverlapRank)
    ;

    class_< DistributedOverlapMatrix<3>, boost::noncopyable >("DistributedOverlapMatrix_3", init<  >())
        .def("MultiplyOverlapRank", &DistributedOverlapMatrix<3>::MultiplyOverlapRank)
    ;

    class_< DistributedOverlapMatrix<4>, boost::noncopyable >("DistributedOverlapMatrix_4", init<  >())
        .def("MultiplyOverlapRank", &DistributedOverlapMatrix<4>::MultiplyOverlapRank)
    ;

}

