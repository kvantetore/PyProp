#ifndef H5TRAITS
#define H5TRAITS

#include <hdf5.h>

template<class T>
class H5Traits
{
public:
	hid_t GetTypeId();
};

template<class T>
class H5Traits< std::complex<T> >
{
public:
	hid_t GetTypeId();
};

template<> inline hid_t H5Traits<int>::GetTypeId()
{
	return H5Tcopy(H5T_NATIVE_INT);
}

template<> inline hid_t H5Traits<double>::GetTypeId()
{
	return H5Tcopy(H5T_NATIVE_DOUBLE);
}

template<> inline hid_t H5Traits<float>::GetTypeId()
{
	return H5Tcopy(H5T_NATIVE_FLOAT);
}

template<class T> inline hid_t H5Traits< std::complex<T> >::GetTypeId()
{
	typedef std::complex<T> cplx_t;
	hid_t baseTypeId = H5Traits<T>().GetTypeId();
	hid_t typeId = H5Tcreate(H5T_COMPOUND, sizeof(cplx_t));
	H5Tinsert(typeId, "r", 0, baseTypeId);
	H5Tinsert(typeId, "i", sizeof(T), baseTypeId);
	H5Tclose(baseTypeId);
	return typeId;
}

#endif

