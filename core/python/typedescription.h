#ifndef TYPEDESCRIPTION_H
#define TYPEDESCRIPTION_H

template<class T>
class PyArrayTraits
{
public:
	typedef T BasicType;
	static PyArray_Descr* GetTypeDescr();
};

template<> PyArray_Descr* PyArrayTraits<double>::GetTypeDescr()
{
	return PyArray_DescrFromType(PyArray_DOUBLE);
}

template<> PyArray_Descr* PyArrayTraits< std::complex<double> >::GetTypeDescr()
{
	return PyArray_DescrFromType(PyArray_CDOUBLE);
}


template<> PyArray_Descr* PyArrayTraits<int>::GetTypeDescr()
{
	return PyArray_DescrFromType(PyArray_INT);
}

template<> PyArray_Descr* PyArrayTraits<long>::GetTypeDescr()
{
	return PyArray_DescrFromType(PyArray_LONG);
}


#endif

