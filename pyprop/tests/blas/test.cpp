
#include <core/utility/blitzblas.cpp>
namespace blas
{
#include <core/utility/blitzblas_cblas.cpp>
}
namespace ref
{
#include <core/utility/blitzblas_ref.cpp>
}


typedef double T;

int main(int argc, char* argv[])
{
	cout << "Test" << endl;

	Array<T, 2> A(2,3);
	Array<T, 2> B(3,2);
	Array<T, 2> CRef(2,2), CBlas(2,2);
	Array<T, 1> v(3);
	Array<T, 1> wRef(2), wBlas(2);

	A = 1,2,3,4,5,6;
	B = 10,20,30,40,50,60;
	v = 5,6,7;

	blas::MatrixVectorMultiply(A, v, wBlas);
	cout << "w (blas) " << wBlas << endl;

	blas::MatrixMatrixMultiply(A, B, CBlas);
	cout << "C (blas) " << CBlas << endl;
	
	ref::MatrixVectorMultiply(A, v, wRef);
	cout << "w (ref): " << wRef << endl;

	ref::MatrixMatrixMultiply(A, B, CRef);
	cout << "C (ref) " << CRef << endl;

	if (any(wBlas != wRef)) 
	{
		cout << "Error in MatrixVectorMultiply" << endl;
	}

	if (any(CBlas != CRef))
	{
		cout << "Error in MatrixMatrixMultiply" << endl;
	}
}


