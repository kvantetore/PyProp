#ifndef BLITZHDF
#define BLITZHDF

#include <stdexcept>

#include <blitz/array.h>
#include <hdf5.h>

#include "h5traits.h"

typedef std::complex<double> cplx;

template<int Rank> hid_t CreateDataSpace(const blitz::TinyVector<int, Rank> &shape)
{
	blitz::TinyVector<hsize_t, Rank> dims;
	for (int i=0; i<Rank; i++)
	{
		dims(i) = shape(i);
	}
	return H5Screate_simple(Rank, dims.data(), NULL);
}


template<class T, int Rank> void GetDataspaceDims(hid_t dataspaceId, blitz::TinyVector<T, Rank> &dims, blitz::TinyVector<T, Rank> &maxdims)
{
	hsize_t dimsPtr[Rank];
	hsize_t maxdimsPtr[Rank];
	int dataspaceRank = H5Sget_simple_extent_dims(dataspaceId, dimsPtr, maxdimsPtr);
	if (dataspaceRank != Rank)
	{
		throw std::runtime_error("Invalid number of dimensions");
	}
	for (int i=0; i<Rank; i++)
	{
		dims[i] = (T) dimsPtr[i];
		maxdims[i] = (T) maxdimsPtr[i];
	}
}

template<class T, int Rank> hid_t CreateDataset(hid_t fileId, const char* name, const blitz::TinyVector<int, Rank> &shape)
{
	hid_t dataSpaceId = CreateDataSpace(shape);
	hid_t dataTypeId = H5Traits<T>().GetTypeId();

	//Create dataset without compression
	hid_t datasetId = H5Dcreate(fileId, name, dataTypeId, dataSpaceId, H5P_DEFAULT);

	H5Tclose(dataTypeId);
	H5Sclose(dataSpaceId);

	return datasetId;
}

template<class T, int Rank> hid_t GetDataset(hid_t fileId, const char* name, const blitz::TinyVector<int, Rank> &shape)
{
	hid_t datasetId = H5Dopen(fileId, name);
	if (datasetId < 0)
	{
		datasetId = CreateDataset<T, Rank>(fileId, name, shape);
	}
	return datasetId;
}


template<class T, int Rank> bool CheckDatasetDims(hid_t datasetId, const blitz::Array<T, Rank> &array)
{
	//Get Dims of dataset
	blitz::TinyVector<int, Rank> dims, maxdims;
	hid_t dataspaceId = H5Dget_space(datasetId);
	GetDataspaceDims(dataspaceId, dims, maxdims);
	H5Sclose(dataspaceId);

	//Verify dims
	for (int i=0; i<Rank; i++)
	{
		if (dims[i] != array.shape()[i])
		{
			return false;
		}
	}
	return true;
}


template<class T, int Rank> herr_t WriteArray(hid_t datasetId, const blitz::Array<T, Rank> &array)
{
	if (!CheckDatasetDims(datasetId, array))
	{
		throw std::runtime_error("Cannot write to datataset, incorrect dims");
	}
	hid_t dataTypeId = H5Traits<T>().GetTypeId();
	herr_t status = H5Dwrite(datasetId, dataTypeId, H5S_ALL, H5S_ALL, H5P_DEFAULT, array.data());
	H5Tclose(dataTypeId);

	return status;
}

template<class T, int Rank> herr_t ReadArray(hid_t datasetId, blitz::Array<T, Rank> &array)
{
	if (!CheckDatasetDims(datasetId, array))
	{
		throw std::runtime_error("Cannot read entire datataset, incorrect dims");
	}
	hid_t dataTypeId = H5Traits<T>().GetTypeId();
	herr_t status = H5Dread(datasetId, dataTypeId, H5S_ALL, H5S_ALL, H5P_DEFAULT, array.data());
	H5Tclose(dataTypeId);

	return status;
}

template<class T, int Rank> blitz::Array<T, Rank> ReadArray(hid_t datasetId)
{
	//Get Dims of dataset
	blitz::TinyVector<int, Rank> dims, maxdims;
	hid_t dataspaceId = H5Dget_space(datasetId);
	GetDataspaceDims(dataspaceId, dims, maxdims);
	H5Sclose(dataspaceId);

	//Allocate and read data
	blitz::Array<T, Rank> array(maxdims);
	ReadArray(datasetId, array);

	return array;
}

template<class T> T ReadScalar(hid_t datasetId)
{
	T value;

	//Check that dataset is zero-dimensional
	hid_t dataspaceId = H5Dget_space(datasetId);
	int dataspaceRank = H5Sget_simple_extent_ndims(dataspaceId);
	H5Sclose(dataspaceId);
	if (dataspaceRank != 0)
	{
		throw std::runtime_error("Dataspace is not a scalar (not of dimension zero)");
	}

	//Read value
	hid_t dataTypeId = H5Traits<T>().GetTypeId();
	herr_t status = H5Dread(datasetId, dataTypeId, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
	H5Tclose(dataTypeId);

	return value;
}

template<int Rank> hid_t GetHyperslabSelection(hid_t datasetId, const blitz::TinyVector<int, Rank> &shape, const blitz::TinyVector<int, Rank> &offset)
{
	hid_t dataSpaceId = H5Dget_space(datasetId);

	blitz::TinyVector<hsize_t, Rank> slabShape = shape;
	blitz::TinyVector<hsize_t, Rank> slabOffset = offset;
	
	H5Sselect_hyperslab(dataSpaceId, H5S_SELECT_SET, slabOffset.data(), NULL, slabShape.data(), NULL);
	return dataSpaceId;
}

template<class T, int Rank> herr_t ReadArrayFromSlab(hid_t datasetId, blitz::Array<T, Rank> &array, const blitz::TinyVector<int, Rank> &offset)
{
	hid_t dataTypeId = H5Traits<T>().GetTypeId();

	//crate selections
	hid_t fileSelectionId = GetHyperslabSelection(datasetId, array.shape(), offset);
	hid_t memorySelectionId = CreateDataSpace(array.shape());

	//read data
	herr_t status = H5Dread(datasetId, dataTypeId, memorySelectionId, fileSelectionId, H5P_DEFAULT, array.data());

	//clean up
	H5Sclose(fileSelectionId);
	H5Sclose(memorySelectionId);
	H5Tclose(dataTypeId);

	return status;
}

template<class T, int Rank> herr_t WriteArrayToSlab(hid_t datasetId, const blitz::Array<T, Rank> &array, const blitz::TinyVector<int, Rank> &offset)
{
	hid_t dataTypeId = H5Traits<T>().GetTypeId();

	//Create selections
	hid_t fileSelectionId = GetHyperslabSelection(datasetId, array.shape(), offset);
	hid_t memorySelectionId = CreateDataSpace(array.shape());

	//Write data
	herr_t status = H5Dwrite(datasetId, dataTypeId, memorySelectionId, fileSelectionId, H5P_DEFAULT, array.data());

	//Clean up
	H5Sclose(fileSelectionId);
	H5Sclose(memorySelectionId);
	H5Tclose(dataTypeId);	

	return status;
}


inline hid_t GetWritableFile(const char* filename)
{
	hid_t fileId = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
	if (fileId < 0)
	{
		fileId = H5Fcreate(filename, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
	}
	return fileId;
}

#endif

