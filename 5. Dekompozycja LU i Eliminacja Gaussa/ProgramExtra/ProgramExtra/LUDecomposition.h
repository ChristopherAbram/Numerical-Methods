#ifndef __LUDECOMPOSITION
#define __LUDECOMPOSITION

#include<vector>
#include "Matrix.h"

#define BASIC_GAUSSIAN				1
#define PARTIAL_SELECTION			2
#define SCALABLE_CHOICE				3
#define FULL_SELECTION				4

/* LU - decomposes matrix U to LU representation.
* Matrix L has ones on its diagonal (L[i][i] = 1 for each i = 0, 1, ..., n).
* @params:
*	Matrix<T>& L - lower matrix,
*	Matrix<T>& U - upper matrix. Its also input matrix - which will be decomposed,
*	std::vector<long>& p - permutated set of indexes,
*   int method - special flag specifies algortihm for LU decomposition.
* @return: void
* @throws: -
*/
template<typename T>
void LU(Matrix<T>& L, Matrix<T>& U, std::vector<long>& p, int method = BASIC_GAUSSIAN) {
	if (U.empty())
		throw *(new MatrixException("Matrix is empty"));

	// Size matrix:
	long n = U.getRows();
	long m = U.getColumns();
	// Array of permutation (0, 1, 2, ..., n-1):
	for (long i = 0; i < n; ++i)
		p[i] = i;

	// First method is the most basic method, where base element is always U[k][k]:
	if (method == BASIC_GAUSSIAN) {
		for (long k = 0; k < n; ++k) {
			L[k][k] = T{ 1. };
			for (long i = k + 1; i < n; ++i) {
				T z = U[i][k] / U[k][k];
				L[i][k] = z;
				U[i][k] = T{};
				for (long j = k + 1; j < m; ++j)
					U[i][j] = U[i][j] - z * U[k][j];
			}
		}
	}

	// Second method is Gaussian elimination using partial selection of the base element,
	// which is such j >= k that |U[p[j]][k]| >= |U[p[i]][k]| for each i = k, k + 1, ..., n.
	else if (method == PARTIAL_SELECTION) {
		for (long k = 0; k < n; ++k) {
			// Selection of the base element:
			long max = k;
			for (long j = k + 1; j < n; ++j) {
				if (abs(U[p[j]][k]) > abs(U[p[max]][k]))
					max = j;
			}
			std::swap(p[k], p[max]);

			L[p[k]][k] = T{ 1. };
			for (long i = k + 1; i < n; ++i) {
				T z = U[p[i]][k] / U[p[k]][k];
				L[p[i]][k] = z;
				U[p[i]][k] = T{};
				for (long j = k + 1; j < m; ++j)
					U[p[i]][j] = U[p[i]][j] - z * U[p[k]][j];
			}
		}
	}
	else if (method == SCALABLE_CHOICE) {
		// TODO: implementation...
	}
	else if (method == FULL_SELECTION) {
		// TODO: implementation...
	}

	return;
}// end LU


 /* LU2 - decomposes matrix A to LU representation.
 * Matrix L and U are saved in matrix A.
 * @params:
 *	Matrix<T>& A - matrix,
 *	std::vector<long>& p - permutated set of indexes,
 *  int method - special flag specifies algortihm for LU decomposition.
 * @return: void
 * @throws: -
 */
template<typename T>
void LU2(Matrix<T>& A, std::vector<long>& p, int method = BASIC_GAUSSIAN) {
	if (A.empty())
		throw *(new MatrixException("Matrix is empty"));

	// Size of matrix:
	long n = A.getRows();
	long m = A.getColumns();
	// Array of permutation (0, 1, 2, ..., n-1):
	for (long i = 0; i < n; ++i)
		p[i] = i;

	// First method is the most basic method, where base element is always U[k][k]:
	if (method == BASIC_GAUSSIAN) {

	}

	// Second method is Gaussian elimination using partial selection of the base element,
	// which is such j >= k that |U[p[j]][k]| >= |U[p[i]][k]| for each i = k, k + 1, ..., n.
	else if (method == PARTIAL_SELECTION) {
		for (long k = 0; k < n; ++k) {
			// Selection of the base element:
			long max = k;
			for (long j = k + 1; j < n; ++j) {
				if (abs(A[p[j]][k]) > abs(A[p[max]][k]))
					max = j;
			}
			std::swap(p[k], p[max]);

			for (long i = k + 1; i < n; ++i) {
				T z = A[p[i]][k] / A[p[k]][k];
				A[p[i]][k] = z;
				for (long j = k + 1; j < m; ++j)
					A[p[i]][j] = A[p[i]][j] - z * A[p[k]][j];
			}
		}
	}
	else if (method == SCALABLE_CHOICE) {
		// TODO: implementation...
	}
	else if (method == FULL_SELECTION) {
		// TODO: implementation...
	}

	return;
}// end LU2

template <typename T>
Matrix<T> Kronecker(const std::vector<long>& i, const std::vector<long>& j) {
	Matrix<T> delta(i.size(), j.size(), T{});
	for (long k = 0; k < i.size(); ++k)
		for (long l = 0; l < j.size(); ++l)
			if (i[k] == j[l])
				delta[k][l] = T{ 1 };
	return delta;
}// end Kronecker

#endif /* __LUDECOMPOSITION */