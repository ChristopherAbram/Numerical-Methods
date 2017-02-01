#ifndef __SOLVINGLINEARSYSTEM
#define __SOLVINGLINEARSYSTEM

#include<vector>
#include "Matrix.h"

class SLSException : public std::runtime_error {
	public:
		explicit SLSException(const std::string& what_str): std::runtime_error(what_str) {}
};

template<typename T>
class SLS {
		using value_type = T;

		// Main A matrix:
		Matrix<value_type> __A;

		// b matrix:
		Matrix<value_type> __b;

		// Matrix A as vector - only for diagonal A:
		std::vector<value_type> __Ad;

		// 3 vectors for Thomas' algorithm:
		std::vector<value_type> __u;
		std::vector<value_type> __d;
		std::vector<value_type> __l;

		// Aliases:
		using self = SLS;
		using self_reference = SLS&;
		using const_self_reference = const SLS&;

	public:
		explicit SLS(const Matrix<value_type>& A, const Matrix<value_type>& b) : __A{ A }, __b{b} {
			if(!is_square(A))
				throw *(new SLSException("SLS::SLS(): Matrix is not square or is empty"));
			if(A.getRows() != b.getRows())
				throw *(new SLSException("SLS::SLS(): Matrix A must have the same amount of rows as matrix b"));
		}
		explicit SLS(const std::vector<value_type>& A, const Matrix<value_type>& b) : __Ad{ A }, __b{b} {
			if (A.size() != b.getRows())
				throw *(new SLSException("SLS::SLS(): Matrix A must have the same amount of rows as matrix b"));
		}
		explicit SLS(const std::vector<value_type>& l, const std::vector<value_type>& d, const std::vector<value_type>& u) : __l{ l }, __d{ d }, __u{ u } { }

		self_reference solveDiagonal(void) {
			// TODO: write implementation...

			return *this;
		}// end solveDiagonal

		self_reference solveLowerTriangular(void) {
			long n = __A.getRows();

			for (long i = 0; i < n; ++i) {
				value_type sum{};
				for (long j = 0; j < i; ++j)
					sum += __A[i][j] * __b[j][0];
				__b[i][0] = (__b[i][0] - sum) / __A[i][i];
			}
			return *this;
		}// end solveLowerTriangular

		self_reference solveUpperTriangular(void) {
			long n = __A.getRows();

			for (long i = n - 1; i >= 0; --i) {
				value_type sum{};
				for (long j = i + 1; j < n; ++j)
					sum += __A[i][j] * __b[j][0];
				__b[i][0] = (__b[i][0] - sum) / __A[i][i];
			}
			return *this;
		}// end solveUpperTriangular

		self_reference solveGaussianElimination(std::vector<long>& p = std::vector<long>()) {
			long n = __A.getRows();
			long m = __A.getColumns();

			p.resize(n);
			LU2(__A, p, PARTIAL_SELECTION);
			Matrix<value_type> x(n, 1);

			// Solving Lz = Pb:
			for (long k = 0; k < n - 1; ++k)
				for (long i = k + 1; i < n; ++i)
					__b[p[i]][0] = __b[p[i]][0] - __A[p[i]][k] * __b[p[k]][0];
				
			// Solving Ux = z:
			for (long i = n - 1; i >= 0; --i) {
				value_type sum{};
				for (long j = i + 1; j < n; ++j)
					sum += __A[p[i]][j] * x[j][0];

				x[i][0] = (__b[p[i]][0] - sum) / __A[p[i]][i];
			}
			// Saving results:
			__b = x;
			return *this;
		}// end solveGaussianElimination

		void LUThomas() {
			long n = __d.size();
			for (long i = 1; i < n; ++i)
				__d[i] = __d[i] - __l[i] * __u[i - 1] / __d[i - 1];
			return;
		}

		void bidiagonalThomas() {
			long n = __d.size();
			std::vector<value_type> x(n);

			for (long i = 1; i < n; ++i)
				__b[i] = __b[i] - __l[i] * __b[i - 1] / __d[i - 1];

			x[n - 1] = __b[n - 1] / __d[n - 1];
			for (long i = n - 2; i >= 0; --i)
				x[i] = (__b[i] - __u[i] * x[i + 1]) / __d[i];

			__b = x;
			return;
		}

		self_reference solveThomas() {
			LUThomas();
			bidiagonalThomas();
			return *this;
		}// end solveThomas

		self_reference setA(const Matrix<value_type>& A) {
			if (!is_square(A))
				throw *(new SLSException("SLS::setA(): Matrix is not square or is empty"));
			__A = A;
			return *this;
		}// end setA

		self_reference setA(const std::vector<value_type>& A) {
			if (A.empty())
				throw *(new SLSException("SLS::setA(): Matrix can not be empty"));
			__Ad = A;
			return *this;
		}// end setA

		self_reference setb(const Matrix<value_type>& b) {
			if ((A.size() != b.getRows()) || (A.getRows() != b.getRows()))
				throw *(new SLSException("SLS::setb(): Matrix b must have the same amount of rows as matrix A"));
			__b = b;
			return *this;
		}// end setb

		Matrix<value_type>& getSolutions(void) {
			return __b;
		}// end getSolutions
};

#endif /* __SOLVINGLINEARSYSTEM */