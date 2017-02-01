#ifndef __SOLVINGLINEARSYSTEM
#define __SOLVINGLINEARSYSTEM

#include<vector>
#include<cmath>  // nan()

#include "Matrix.h"
#include "LUDecomposition.h"
#include "Result.h"

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

		// Beginning point (only for itaration algorithm):
		vector_reference __x0;

		// Whether to save error estimation or not:
		bool __approximation = false;
		// Result with errors:
		std::vector<Result<value_type>> __results;

		// Iteration parameters:
		static value_type __delta;		// TOL
		static value_type __epsilon;	// TOLF
		static int __M;					// max amout of iteraions 

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

		using value_type = T;

		explicit SLS(matrix_reference A, vector_reference b, vector_reference x0 = __vtmp) : __A{ A }, __b{ b }, __l{ __vtmp }, __d{ __vtmp }, __u{ __vtmp }, __Ad{ __vtmp }, __x0{x0} {
			if(!is_square(A))
				throw *(new SLSException("SLS::SLS(): Matrix is not square or is empty"));
			if(A.getRows() != b.size())
				throw *(new SLSException("SLS::SLS(): Matrix A must have the same amount of rows as matrix b"));
		}


		explicit SLS(vector_reference A, vector_reference b) : __Ad{ A }, __b{b}, __A{ __tmp }, __l{ __vtmp }, __d{ __vtmp }, __u{ __vtmp }, __x0{ __vtmp } {
			if (A.size() != b.size())
				throw *(new SLSException("SLS::SLS(): Matrix A must have the same amount of rows as matrix b"));
		}


		explicit SLS(vector_reference l, vector_reference d, vector_reference u) : __l{ l }, __d{ d }, __u{ u }, __A{ __tmp }, __b{ __vtmp }, __Ad{__vtmp}, __x0{ __vtmp } { }


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
		self_reference solveThomas(void) {
			__LUThomas();
			__bidiagonalThomas();

			return *this;
		}// end solveThomas










		/* ****************************************************************************** */
		/* **************************   Iterative Algorithms   ************************** */
		/* ****************************************************************************** */
		
		/** solveJacobi - solves linear system using Jacbi method (iterative mathod).
		* @return	SLS&	- reference to self (fluent interface)
		* @throws	-
		*/
		self_reference solveJacobi(void) {
			long n = __A.getRows();
			std::vector<value_type> u(n);

			__results.clear();

			// Error Estimators:
			value_type EST{};
			value_type RES{};

			long k;
			for (k = 0; k < getMax(); ++k) {
				for (long i = 0; i < n; ++i) {
					value_type sum{};
					for (long j = 0; j < n; ++j)
						if (i != j)
							sum += __A[i][j] * __x0[j];
					u[i] = (__b[i] - sum) / __A[i][i];
				}

				// Calculating estimators:
				EST = norm_max(u - __x0);
				RES = norm_max(__A * __x0 - __b);

				for (long i = 0; i < n; ++i)
					__x0[i] = u[i];

				if (EST <= getDelta() && RES <= getEpsilon()) {
					if (__approximation)
						__results.push_back(Result<value_type>(k + 1, (EST <= getDelta() ? DELTA_REACHED : 0) | (RES <= getEpsilon() ? EPSILON_REACHED : 0), __x0, RES, EST));
					break;
				}
				else if(__approximation)
					__results.push_back(Result<value_type>(k + 1, 0, __x0, RES, EST));
			}

			if(k + 1 == getMax() && __approximation)
				__results.push_back(Result<value_type>(k + 1, ITER_REACHED, __x0, RES, EST));

			return *this;
		}// end solveJacobi


		 /** solveGaussSeidel - solves linear system using Gauss - Seidel method (iterative mathod).
		 * @return	SLS&	- reference to self (fluent interface)
		 * @throws	-
		 */
		self_reference solveGaussSeidel(void) {
			long n = __A.getRows();
			std::vector<value_type> u(n);

			__results.clear();

			// Error Estimators:
			value_type EST{};
			value_type RES{};

			long k;
			for (k = 0; k < getMax(); ++k) {
				
				for (long i = 0; i < n; ++i) {
					value_type sum{};
					for (long j = 0; j < n; ++j) {
						if (j != i)
							sum += __A[i][j] * __x0[j];
					}
					u[i] = __x0[i];
					__x0[i] = (__b[i] - sum) / __A[i][i];
				}

				// Calculating estimators:
				EST = norm_max(u - __x0);
				RES = norm_max(__A * __x0 - __b);

				if (EST <= getDelta() && RES <= getEpsilon()) {
					if (__approximation)
						__results.push_back(Result<value_type>(k + 1, (EST <= getDelta() ? DELTA_REACHED : 0) | (RES <= getEpsilon() ? EPSILON_REACHED : 0), __x0, RES, EST));
					break;
				}
				else if (__approximation)
					__results.push_back(Result<value_type>(k + 1, 0, __x0, RES, EST));
			}

			if (k + 1 == getMax() && __approximation)
				__results.push_back(Result<value_type>(k + 1, ITER_REACHED, __x0, RES, EST));

			return *this;
		}// end solveGaussSeidel


		 /** solveSOR - solves linear system using successive over-relaxation method (iterative mathod).
		 * @return	SLS&	- reference to self (fluent interface)
		 * @throws	-
		 */
		self_reference solveSOR(value_type w = value_type{}) {
			long n = __A.getRows();
			std::vector<value_type> u(n);

			__results.clear();

			// Error Estimators:
			value_type EST{};
			value_type RES{};

			long k;
			for (k = 0; k < getMax(); ++k) {
				
				for (long i = 0; i < n; ++i) {
					value_type sum{};
					for (long j = 0; j < n; ++j)
						if (j != i)
							sum += __A[i][j] * __x0[j];
					sum = (__b[i] - sum) / __A[i][i];
					u[i] = __x0[i] + w * (sum - __x0[i]);
				}

				// Calculating estimators:
				EST = norm_max(u - __x0);
				RES = norm_max(__A * __x0 - __b);

				for (long i = 0; i < n; ++i)
					__x0[i] = u[i];

				if (EST <= getDelta() && RES <= getEpsilon()) {
					if (__approximation)
						__results.push_back(Result<value_type>(k + 1, (EST <= getDelta() ? DELTA_REACHED : 0) | (RES <= getEpsilon() ? EPSILON_REACHED : 0), __x0, RES, EST));
					break;
				}
				else if (__approximation)
					__results.push_back(Result<value_type>(k + 1, 0, __x0, RES, EST));
			}

			if (k + 1 == getMax() && __approximation)
				__results.push_back(Result<value_type>(k + 1, ITER_REACHED, __x0, RES, EST));

			return *this;
		}// end solveSOR






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

		self_reference setBeginningPoint(vector_reference x0) {
			__x0 = x0;
			return *this;
		}// end setBeginnigPoint

		self_reference areErrorsAvailable(bool a) {
			__approximation = a;
			return *this;
		}// end areErrorsAvailable

		vector_reference getSolutions(void) {
			return __b;
		}// end getSolutions

		std::vector<Result<value_type>>& getResults(void){
			return __results;
		}// end getResults

		static void setMax(int max) {
			__M = max;
			return;
		}// end setMax

		static int getMax(void) {
			return __M;
		}// end getMax

		static void setDelta(value_type d) {
			__delta = d;
			return;
		}// end setDelta

		static void setEpsilon(value_type e) {
			__epsilon = e;
			return;
		}// end setEpsilon

		static value_type getDelta() {
			return __delta;
		}// end getDelata

		static value_type getEpsilon() {
			return __epsilon;
		}// end getEpsilon
};

template<typename T> typename SLS<T>::value_type SLS<T>::__delta = SLS<T>::value_type{};
template<typename T> typename SLS<T>::value_type SLS<T>::__epsilon = SLS<T>::value_type{};
template<typename T> int SLS<T>::__M = 100;

#endif /* __SOLVINGLINEARSYSTEM */