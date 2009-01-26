
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <utility/databuffer.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_python_databuffer()
{
    scope* DataBuffer_1_scope = new scope(
    class_< DataBuffer<1> >("DataBuffer_1", init<  >())
        .def(init< blitz::TinyVector<int,1>& >())
        .def(init< const DataBuffer<1>& >())
        .def("IsAvailable", &DataBuffer<1>::IsAvailable)
        .def("Lock", &DataBuffer<1>::Lock)
        .def("UnLock", &DataBuffer<1>::UnLock)
        .def("GetArray", &DataBuffer<1>::GetArray2)
        .def("ResizeArray", &DataBuffer<1>::ResizeArray)
        .def("FreeArray", &DataBuffer<1>::FreeArray)
        .def( self == other< blitz::TinyVector<int,1> >() )
    );
    register_ptr_to_python< boost::shared_ptr< DataBuffer<1> > >();
    delete DataBuffer_1_scope;

    scope* DataBuffer_2_scope = new scope(
    class_< DataBuffer<2> >("DataBuffer_2", init<  >())
        .def(init< blitz::TinyVector<int,2>& >())
        .def(init< const DataBuffer<2>& >())
        .def("IsAvailable", &DataBuffer<2>::IsAvailable)
        .def("Lock", &DataBuffer<2>::Lock)
        .def("UnLock", &DataBuffer<2>::UnLock)
        .def("GetArray", &DataBuffer<2>::GetArray2)
        .def("ResizeArray", &DataBuffer<2>::ResizeArray)
        .def("FreeArray", &DataBuffer<2>::FreeArray)
        .def( self == other< blitz::TinyVector<int,2> >() )
    );
    register_ptr_to_python< boost::shared_ptr< DataBuffer<2> > >();
    delete DataBuffer_2_scope;

    scope* DataBuffer_3_scope = new scope(
    class_< DataBuffer<3> >("DataBuffer_3", init<  >())
        .def(init< blitz::TinyVector<int,3>& >())
        .def(init< const DataBuffer<3>& >())
        .def("IsAvailable", &DataBuffer<3>::IsAvailable)
        .def("Lock", &DataBuffer<3>::Lock)
        .def("UnLock", &DataBuffer<3>::UnLock)
        .def("GetArray", &DataBuffer<3>::GetArray2)
        .def("ResizeArray", &DataBuffer<3>::ResizeArray)
        .def("FreeArray", &DataBuffer<3>::FreeArray)
        .def( self == other< blitz::TinyVector<int,3> >() )
    );
    register_ptr_to_python< boost::shared_ptr< DataBuffer<3> > >();
    delete DataBuffer_3_scope;

    scope* DataBuffer_4_scope = new scope(
    class_< DataBuffer<4> >("DataBuffer_4", init<  >())
        .def(init< blitz::TinyVector<int,4>& >())
        .def(init< const DataBuffer<4>& >())
        .def("IsAvailable", &DataBuffer<4>::IsAvailable)
        .def("Lock", &DataBuffer<4>::Lock)
        .def("UnLock", &DataBuffer<4>::UnLock)
        .def("GetArray", &DataBuffer<4>::GetArray2)
        .def("ResizeArray", &DataBuffer<4>::ResizeArray)
        .def("FreeArray", &DataBuffer<4>::FreeArray)
        .def( self == other< blitz::TinyVector<int,4> >() )
    );
    register_ptr_to_python< boost::shared_ptr< DataBuffer<4> > >();
    delete DataBuffer_4_scope;

}

