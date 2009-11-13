
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <representation/representation.h>
#include <wavefunction.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

blitz::Array<std::complex<double>,1>& (Wavefunction<1>::*Wavefunction_1___GetData)()  = &Wavefunction<1>::GetData;

const blitz::Array<std::complex<double>,1>& (Wavefunction<1>::*Wavefunction_1___GetData_const)() const = &Wavefunction<1>::GetData;

int (Wavefunction<1>::*Wavefunction_1___AllocateDatablitz__TinyVector_int_1_)(blitz::TinyVector<int,1>)  = &Wavefunction<1>::AllocateData;

blitz::Array<std::complex<double>,1>& (Wavefunction<1>::*Wavefunction_1___GetDataint)(int)  = &Wavefunction<1>::GetData;

const blitz::Array<std::complex<double>,1>& (Wavefunction<1>::*Wavefunction_1___GetDataint_const)(int) const = &Wavefunction<1>::GetData;

void (Wavefunction<1>::*Wavefunction_1___AllocateData)()  = &Wavefunction<1>::AllocateData;

blitz::Array<std::complex<double>,2>& (Wavefunction<2>::*Wavefunction_2___GetData)()  = &Wavefunction<2>::GetData;

const blitz::Array<std::complex<double>,2>& (Wavefunction<2>::*Wavefunction_2___GetData_const)() const = &Wavefunction<2>::GetData;

int (Wavefunction<2>::*Wavefunction_2___AllocateDatablitz__TinyVector_int_2_)(blitz::TinyVector<int,2>)  = &Wavefunction<2>::AllocateData;

blitz::Array<std::complex<double>,2>& (Wavefunction<2>::*Wavefunction_2___GetDataint)(int)  = &Wavefunction<2>::GetData;

const blitz::Array<std::complex<double>,2>& (Wavefunction<2>::*Wavefunction_2___GetDataint_const)(int) const = &Wavefunction<2>::GetData;

void (Wavefunction<2>::*Wavefunction_2___AllocateData)()  = &Wavefunction<2>::AllocateData;

blitz::Array<std::complex<double>,3>& (Wavefunction<3>::*Wavefunction_3___GetData)()  = &Wavefunction<3>::GetData;

const blitz::Array<std::complex<double>,3>& (Wavefunction<3>::*Wavefunction_3___GetData_const)() const = &Wavefunction<3>::GetData;

int (Wavefunction<3>::*Wavefunction_3___AllocateDatablitz__TinyVector_int_3_)(blitz::TinyVector<int,3>)  = &Wavefunction<3>::AllocateData;

blitz::Array<std::complex<double>,3>& (Wavefunction<3>::*Wavefunction_3___GetDataint)(int)  = &Wavefunction<3>::GetData;

const blitz::Array<std::complex<double>,3>& (Wavefunction<3>::*Wavefunction_3___GetDataint_const)(int) const = &Wavefunction<3>::GetData;

void (Wavefunction<3>::*Wavefunction_3___AllocateData)()  = &Wavefunction<3>::AllocateData;

blitz::Array<std::complex<double>,4>& (Wavefunction<4>::*Wavefunction_4___GetData)()  = &Wavefunction<4>::GetData;

const blitz::Array<std::complex<double>,4>& (Wavefunction<4>::*Wavefunction_4___GetData_const)() const = &Wavefunction<4>::GetData;

int (Wavefunction<4>::*Wavefunction_4___AllocateDatablitz__TinyVector_int_4_)(blitz::TinyVector<int,4>)  = &Wavefunction<4>::AllocateData;

blitz::Array<std::complex<double>,4>& (Wavefunction<4>::*Wavefunction_4___GetDataint)(int)  = &Wavefunction<4>::GetData;

const blitz::Array<std::complex<double>,4>& (Wavefunction<4>::*Wavefunction_4___GetDataint_const)(int) const = &Wavefunction<4>::GetData;

void (Wavefunction<4>::*Wavefunction_4___AllocateData)()  = &Wavefunction<4>::AllocateData;


}// namespace 


// Module ======================================================================
void Export_python_wavefunction()
{
    scope* Wavefunction_1_scope = new scope(
    class_< Wavefunction<1>, boost::noncopyable >("Wavefunction_1", init<  >())
        .def_readwrite("Data", &Wavefunction<1>::Data)
        .def("SetRepresentation", &Wavefunction<1>::SetRepresentation)
        .def("GetRepresentation", &Wavefunction<1>::GetRepresentation)
        .def("GetData", Wavefunction_1___GetData, return_value_policy< return_by_value >())
        .def("GetData", Wavefunction_1___GetData_const, return_value_policy< return_by_value >())
        .def("GetMemoryFootprint", &Wavefunction<1>::GetMemoryFootprint)
        .def("AllocateData", Wavefunction_1___AllocateDatablitz__TinyVector_int_1_)
        .def("FreeData", &Wavefunction<1>::FreeData)
        .def("GetActiveBufferName", &Wavefunction<1>::GetActiveBufferName)
        .def("SetActiveBuffer", &Wavefunction<1>::SetActiveBuffer)
        .def("GetData", Wavefunction_1___GetDataint, return_value_policy< return_by_value >())
        .def("GetData", Wavefunction_1___GetDataint_const, return_value_policy< return_by_value >())
        .def("SetData", &Wavefunction<1>::SetData)
        .def("LockBuffer", &Wavefunction<1>::LockBuffer)
        .def("UnLockBuffer", &Wavefunction<1>::UnLockBuffer)
        .def("GetAvailableDataBufferName", &Wavefunction<1>::GetAvailableDataBufferName)
        .def("HasAvailableBuffer", &Wavefunction<1>::HasAvailableBuffer)
        .def("GetRank", &Wavefunction<1>::GetRank)
        .def("AllocateData", Wavefunction_1___AllocateData)
        .def("InnerProduct", &Wavefunction<1>::InnerProduct)
        .def("LocalInnerProduct", &Wavefunction<1>::LocalInnerProduct)
        .def("GetNorm", &Wavefunction<1>::GetNorm)
        .def("GetLocalNorm", &Wavefunction<1>::GetLocalNorm)
        .def("Normalize", &Wavefunction<1>::Normalize)
        .def("Clear", &Wavefunction<1>::Clear)
        .def("Copy", &Wavefunction<1>::Copy)
        .def("CopyDeep", &Wavefunction<1>::CopyDeep)
    );
    register_ptr_to_python< boost::shared_ptr< Wavefunction<1> > >();
    delete Wavefunction_1_scope;

    scope* Wavefunction_2_scope = new scope(
    class_< Wavefunction<2>, boost::noncopyable >("Wavefunction_2", init<  >())
        .def_readwrite("Data", &Wavefunction<2>::Data)
        .def("SetRepresentation", &Wavefunction<2>::SetRepresentation)
        .def("GetRepresentation", &Wavefunction<2>::GetRepresentation)
        .def("GetData", Wavefunction_2___GetData, return_value_policy< return_by_value >())
        .def("GetData", Wavefunction_2___GetData_const, return_value_policy< return_by_value >())
        .def("GetMemoryFootprint", &Wavefunction<2>::GetMemoryFootprint)
        .def("AllocateData", Wavefunction_2___AllocateDatablitz__TinyVector_int_2_)
        .def("FreeData", &Wavefunction<2>::FreeData)
        .def("GetActiveBufferName", &Wavefunction<2>::GetActiveBufferName)
        .def("SetActiveBuffer", &Wavefunction<2>::SetActiveBuffer)
        .def("GetData", Wavefunction_2___GetDataint, return_value_policy< return_by_value >())
        .def("GetData", Wavefunction_2___GetDataint_const, return_value_policy< return_by_value >())
        .def("SetData", &Wavefunction<2>::SetData)
        .def("LockBuffer", &Wavefunction<2>::LockBuffer)
        .def("UnLockBuffer", &Wavefunction<2>::UnLockBuffer)
        .def("GetAvailableDataBufferName", &Wavefunction<2>::GetAvailableDataBufferName)
        .def("HasAvailableBuffer", &Wavefunction<2>::HasAvailableBuffer)
        .def("GetRank", &Wavefunction<2>::GetRank)
        .def("AllocateData", Wavefunction_2___AllocateData)
        .def("InnerProduct", &Wavefunction<2>::InnerProduct)
        .def("LocalInnerProduct", &Wavefunction<2>::LocalInnerProduct)
        .def("GetNorm", &Wavefunction<2>::GetNorm)
        .def("GetLocalNorm", &Wavefunction<2>::GetLocalNorm)
        .def("Normalize", &Wavefunction<2>::Normalize)
        .def("Clear", &Wavefunction<2>::Clear)
        .def("Copy", &Wavefunction<2>::Copy)
        .def("CopyDeep", &Wavefunction<2>::CopyDeep)
    );
    register_ptr_to_python< boost::shared_ptr< Wavefunction<2> > >();
    delete Wavefunction_2_scope;

    scope* Wavefunction_3_scope = new scope(
    class_< Wavefunction<3>, boost::noncopyable >("Wavefunction_3", init<  >())
        .def_readwrite("Data", &Wavefunction<3>::Data)
        .def("SetRepresentation", &Wavefunction<3>::SetRepresentation)
        .def("GetRepresentation", &Wavefunction<3>::GetRepresentation)
        .def("GetData", Wavefunction_3___GetData, return_value_policy< return_by_value >())
        .def("GetData", Wavefunction_3___GetData_const, return_value_policy< return_by_value >())
        .def("GetMemoryFootprint", &Wavefunction<3>::GetMemoryFootprint)
        .def("AllocateData", Wavefunction_3___AllocateDatablitz__TinyVector_int_3_)
        .def("FreeData", &Wavefunction<3>::FreeData)
        .def("GetActiveBufferName", &Wavefunction<3>::GetActiveBufferName)
        .def("SetActiveBuffer", &Wavefunction<3>::SetActiveBuffer)
        .def("GetData", Wavefunction_3___GetDataint, return_value_policy< return_by_value >())
        .def("GetData", Wavefunction_3___GetDataint_const, return_value_policy< return_by_value >())
        .def("SetData", &Wavefunction<3>::SetData)
        .def("LockBuffer", &Wavefunction<3>::LockBuffer)
        .def("UnLockBuffer", &Wavefunction<3>::UnLockBuffer)
        .def("GetAvailableDataBufferName", &Wavefunction<3>::GetAvailableDataBufferName)
        .def("HasAvailableBuffer", &Wavefunction<3>::HasAvailableBuffer)
        .def("GetRank", &Wavefunction<3>::GetRank)
        .def("AllocateData", Wavefunction_3___AllocateData)
        .def("InnerProduct", &Wavefunction<3>::InnerProduct)
        .def("LocalInnerProduct", &Wavefunction<3>::LocalInnerProduct)
        .def("GetNorm", &Wavefunction<3>::GetNorm)
        .def("GetLocalNorm", &Wavefunction<3>::GetLocalNorm)
        .def("Normalize", &Wavefunction<3>::Normalize)
        .def("Clear", &Wavefunction<3>::Clear)
        .def("Copy", &Wavefunction<3>::Copy)
        .def("CopyDeep", &Wavefunction<3>::CopyDeep)
    );
    register_ptr_to_python< boost::shared_ptr< Wavefunction<3> > >();
    delete Wavefunction_3_scope;

    scope* Wavefunction_4_scope = new scope(
    class_< Wavefunction<4>, boost::noncopyable >("Wavefunction_4", init<  >())
        .def_readwrite("Data", &Wavefunction<4>::Data)
        .def("SetRepresentation", &Wavefunction<4>::SetRepresentation)
        .def("GetRepresentation", &Wavefunction<4>::GetRepresentation)
        .def("GetData", Wavefunction_4___GetData, return_value_policy< return_by_value >())
        .def("GetData", Wavefunction_4___GetData_const, return_value_policy< return_by_value >())
        .def("GetMemoryFootprint", &Wavefunction<4>::GetMemoryFootprint)
        .def("AllocateData", Wavefunction_4___AllocateDatablitz__TinyVector_int_4_)
        .def("FreeData", &Wavefunction<4>::FreeData)
        .def("GetActiveBufferName", &Wavefunction<4>::GetActiveBufferName)
        .def("SetActiveBuffer", &Wavefunction<4>::SetActiveBuffer)
        .def("GetData", Wavefunction_4___GetDataint, return_value_policy< return_by_value >())
        .def("GetData", Wavefunction_4___GetDataint_const, return_value_policy< return_by_value >())
        .def("SetData", &Wavefunction<4>::SetData)
        .def("LockBuffer", &Wavefunction<4>::LockBuffer)
        .def("UnLockBuffer", &Wavefunction<4>::UnLockBuffer)
        .def("GetAvailableDataBufferName", &Wavefunction<4>::GetAvailableDataBufferName)
        .def("HasAvailableBuffer", &Wavefunction<4>::HasAvailableBuffer)
        .def("GetRank", &Wavefunction<4>::GetRank)
        .def("AllocateData", Wavefunction_4___AllocateData)
        .def("InnerProduct", &Wavefunction<4>::InnerProduct)
        .def("LocalInnerProduct", &Wavefunction<4>::LocalInnerProduct)
        .def("GetNorm", &Wavefunction<4>::GetNorm)
        .def("GetLocalNorm", &Wavefunction<4>::GetLocalNorm)
        .def("Normalize", &Wavefunction<4>::Normalize)
        .def("Clear", &Wavefunction<4>::Clear)
        .def("Copy", &Wavefunction<4>::Copy)
        .def("CopyDeep", &Wavefunction<4>::CopyDeep)
    );
    register_ptr_to_python< boost::shared_ptr< Wavefunction<4> > >();
    delete Wavefunction_4_scope;

}

