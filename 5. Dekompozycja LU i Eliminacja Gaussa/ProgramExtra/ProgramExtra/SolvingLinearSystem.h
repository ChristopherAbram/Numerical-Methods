#ifndef __SOLVINGLINEARSYSTEM
#define __SOLVINGLINEARSYSTEM

#include<vector>
#include<cmath>  // nan()

#include "Matrix.h"
#include "LUDecomposition.h"

class SLSException : public std::runtime_error {
	public:
		explicit SLSException(const std::string& what_str): std::runtime_error(what_str) {}
};

/* Solving Linear System */
template<typename T>
class SLS {
		using value_type = T;
		using vector = std::vector<value_type>;
		using vector_reference = vector&;
		using matrix = Matrix<value_type>;
		using matrix_reference = matrix&;

		matrix __tmp;
		vector __vtmp;

		// Main A matrix:
		matrix_reference __A;

		// b matrix:
		vector_reference __b;

		// Matrix A as vector - only for diagonal A:
		vector_reference __Ad;

		// 3 vectors for Thomas' algorithm:
		vector_reference __u;
		vector_reference __d;
		vector_reference __l;

		// Aliases:
		using self = SLS;
		using self_reference = SLS&;
		using const_self_reference = const SLS&;


		/** __LUThomas - make LU decomposition for tridiagonal matrix (band matrix).
		* @return	-
		* @throws	-
		*/
		void __LUThomas() {
			long n = __d.size();
			for (long i = 1; i < n; ++i)
				__d[i] = __d[i] - (__l[i - 1] / __d[i - 1]) * __u[i - 1];
			return;
		}// end __LUThomas

		 /** __bidiagonalThomas - solves system using Thomas algorithm.
		 * @return	-
		 * @throws	-
		 */
		void __bidiagonalThomas() {
			long n = __d.size();
			//std::vector<value_type> x(n);

			for (long i = 1; i < n; ++i)
				__b[i] = __b[i] - (__l[i - 1] / __d[i - 1]) * __b[i - 1];

			__b[n - 1] = __b[n - 1] / __d[n - 1];
			for (long i = n - 2; i >= 0; --i)
				__b[i] = (__b[i] - __u[i] * __b[i + 1]) / __d[i];
			
			return;
		}// end __bidiagonalThomas

	public:


		explicit SLS(matrix_reference A, vector_reference b) : __A{ A }, __b{b}, __l{ __vtmp }, __d{ __vtmp }, __u{ __vtmp }, __Ad{ __vtmp } {
			if(!is_square(A))
				throw *(new SLSException("SLS::SLS(): Matrix is not square or is empty"));
			if(A.getRows() != b.size())
				throw *(new SLSException("SLS::SLS(): Matrix A must have the same amount of rows as matrix b"));
		}


		explicit SLS(vector_reference A, vector_reference b) : __Ad{ A }, __b{b}, __A{ __tmp }, __l{ __vtmp }, __d{ __vtmp }, __u{ __vtmp } {
			if (A.size() != b.size())
				throw *(new SLSException("SLS::SLS(): Matrix A must have the same amount of rows as matrix b"));
		}


		explicit SLS(vector_reference l, vector_reference d, vector_reference u) : __l{ l }, __d{ d }, __u{ u }, __A{ __tmp }, __b{ __vtmp }, __Ad{__vtmp} { }


		/** solveDiagonal - solves linear system when matrix A has non-zero values only on its diagonal,
		* and zero values above and under diagonal.
		* @return	SLS& - reference to self (fluent interface)
		* @throws	-
		*/
		self_reference solveDiagonal(void) {
			long n = __Ad.size();

			for (long i = 0; i < n; ++i)
				if (__Ad[i] != value_type{})
					__b[i] = __b[i] / __Ad[i];
				else
					__b[i] = (value_type)nan("");

			return *this;
		}// end solveDiagonal


		 /** solveLowerTriangular - solves linear system when matrix A has non-zero values under its diagonal,
		 * and zero values above.
		 * @return	SLS& - reference to self (fluent interface)
		 * @throws	-
		 */
		self_reference solveLowerTriangular(void) {
			long n = __A.getRows();

			for (long i = 0; i < n; ++i) {
				value_type sum{};
				for (long j = 0; j < i; ++j)
					sum += __A[i][j] * __b[j];
				__b[i] = (__b[i] - sum) / __A[i][i];
			}
			return *this;
		}// end solveLowerTriangular


		 /** solveUpperTriangular - solves linear system when matrix A has non-zero values above its diagonal,
		 * and zero values under.
		 * @return	SLS& - reference to self (fluent interface)
		 * @throws	-
		 */
		self_reference solveUpperTriangular(void) {
			long n = __A.getRows();

			for (long i = n - 1; i >= 0; --i) {
				value_type sum{};
				for (long j = i + 1; j < n; ++j)
					sum += __A[i][j] * __b[j];
				__b[i] = (__b[i] - sum) / __A[i][i];
			}
			return *this;
		}// end solveUpperTriangular


		 /** solveGaussianElimination - solves linear system for practically all matrixes.
		 * @return	SLS& - reference to self (fluent interface)
		 * @throws	-
		 */
		self_reference solveGaussianElimination(std::vector<long>& p = std::vector<long>()) {
			// Size of matrix A:
			long n = __A.getRows();
			long m = __A.getColumns();
			p.resize(n);

			// LU decomposition of matrix A:  PA = LU
			LU2(__A, p, PARTIAL_SELECTION);

			// Solving system Ax = b  --> LUx = b
			// 1) Lz = Pb    (P = pI)
			// 2) Ux = z
			
			// Solving Lz = Pb:
			for (long k = 0; k < n - 1; ++k)
				for (long i = k + 1; i < n; ++i)
					__b[p[i]] = __b[p[i]] - __A[p[i]][k] * __b[p[k]];
				
			// Solving Ux = z:
			for (long i = n - 1; i >= 0; --i) {
				value_type sum{};
				for (long j = i + 1; j < n; ++j)
					sum += __A[p[i]][j] * __b[p[j]];

				__b[p[i]] = (__b[p[i]] - sum) / __A[p[i]][i];
			}

			return *this;
		}// end solveGaussianElimination


		 /** solveThomas - solves linear system when matrix A is band (tridiagonal) matrix.
		 * @return	SLS& - reference to self (fluent interface)
		 * @throws	-
		 */
		self_reference solveThomas() {
			__LUThomas();
			__bidiagonalThomas();

			return *this;
		}// end solveThomas
		

		self_reference setA(matrix_reference A) {
			if (!is_square(A))
				throw *(new SLSException("SLS::setA(): Matrix is not square or is empty"));
			__A = A;
			return *this;
		}// end setA

		self_reference setA(vector_reference A) {
			if (A.empty())
				throw *(new SLSException("SLS::setA(): Matrix can not be empty"));
			__Ad = A;
			return *this;
		}// end setA

		self_reference setb(vector_reference b) {
			//if ((A.size() != b.getRows()) || (A.getRows() != b.getRows()))
			//	throw *(new SLSException("SLS::setb(): Matrix b must have the same amount of rows as matrix A"));
			__b = b;
			return *this;
		}// end setb

		vector_reference getSolutions(void) {
			return __b;
		}// end getSolutions
};

#endif /* __SOLVINGLINEARSYSTEM */