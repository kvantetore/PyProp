#ifndef TYPEDESCRIPTION_H
#define TYPEDESCRIPTION_H

/* 
Methods to convert from C type to NumPy type descriptor. If anyone knows
of any smart template tricks to avoid having T as an argument, please let me know
*/
template<class T> PyArray_Descr* type_to_descr(T)
{
	return 0;
}

template<> PyArray_Descr* type_to_descr(double)
{
	return PyArray_DescrFromType(PyArray_DOUBLE);
}

template<> PyArray_Descr* type_to_descr(float)
{
	return PyArray_DescrFromType(PyArray_FLOAT);
}

template<> PyArray_Descr* type_to_descr(int)
{
	return PyArray_DescrFromType(PyArray_LONG);
}

template<> PyArray_Descr* type_to_descr(std::complex<double>)
{
	return PyArray_DescrFromType(PyArray_CDOUBLE);
}

template<> PyArray_Descr* type_to_descr(std::complex<float>)
{
	return PyArray_DescrFromType(PyArray_CFLOAT);
}

#endif
