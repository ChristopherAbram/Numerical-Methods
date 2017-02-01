#ifndef __MATRIX
#define __MATRIX

#include<vector>
#include<initializer_list>

class MatrixException : public std::runtime_error {
	public:
		explicit MatrixException(const std::string& what_str): std::runtime_error(what_str) {}
};

template <typename T>
class Matrix {
	private:

		// Vector of vectors:
		std::vector<std::vector<T>> __matrix;
		// Dimensions:
		long __n, __m;
		// Not typical matrix:
		static Matrix<T> Zero;

		// Aliases:
		using Type = T;
		using self = Matrix<T>;

	public:

		Matrix(std::initializer_list<std::vector<Type>> list) {
			std::initializer_list<std::vector<Type>>::iterator b = list.begin();
			std::initializer_list<std::vector<Type>>::const_iterator e = list.end();

			__matrix.reserve(list.size()); // reserves rows
			long n = 0, max = b->size();
			while (b != e) {
				__matrix.push_back(*b);
				++b;
				++n;
				if (max < b->size())
					max = b->size();
			}
			__n = n;
			__m = max;
		}

		explicit Matrix(long n = 1, long m = 1, Type init = Type{}) {
			__matrix = std::vector<std::vector<Type>>(n, std::vector<Type>(m, init));
			__n = n;
			__m = m;
		}

		Matrix(const std::vector<Type>& ob) : Matrix(1, ob.size()) {
			for (long i = 0; i < getColumns(); i++)
				__matrix[0][i] = ob[i];
		}

		/* Position - struct represents position of one cell */
		struct Position {
			long row;
			long column;
		};

		std::vector<Type>& operator[](long n) {
			return __matrix[n];
		}// end operator[]

		std::vector<Type> operator[](long n) const {
			return __matrix[n];
		}// end operator[]

		self& operator=(const std::vector<Type>& ob) {
			resize(1, ob.size());
			for (long i = 0; i < getColumns(); i++)
				__matrix[0][i] = ob[i];
			return *this;
		}

		void resize(long n, long m, Type init = Type{}) {
			if (n == getRows() && m == getColumns())
				return;
			else {
				__matrix.resize(n);
				for (long i = 0; i < __matrix.size(); i++)
					__matrix[i].resize(m, init);
				__n = __matrix.size();
				__m = __matrix[0].size();
			}
			return;
		}

		void appendRows(const self& M, long rows = 0) {
			if (M.getColumns() > getColumns())
				throw *(new MatrixException("Matrix::appendRows: Unable to append Matrix."));

			rows = rows == 0 ? M.getRows() : abs(rows);
			long n = getRows();
			long m = getColumns();
			resize(n + rows, m);

			for (long i = n; i < getRows(); ++i)
				for (long j = 0; j < M.getColumns(); ++j)
					__matrix[i][j] = M[i - n][j];

			return;
		}// end appendRows

		self row(long n) const {
			self R = __matrix[n];
			return R;
		}// end row

		void appendColumns(const self& M, long columns = 0) {
			if (M.getRows() > getRows())
				throw *(new MatrixException("Matrix::appendColumns: Unable to append Matrix."));

			columns = columns == 0 ? M.getColumns() : abs(columns);
			long n = getRows();
			long m = getColumns();
			resize(n, m + columns);

			for (long i = 0; i < M.getRows(); ++i)
				for (long j = m; j < getColumns(); ++j)
					__matrix[i][j] = M[i][j - m];

			return;
		}// end appendColumns

		self column(long n) const {
			self C(getRows(), 1);
			for (long i = 0; i < getRows(); i++)
				C[i][0] = __matrix[i][n];
			return C;
		}// end row

		self& transpose() {
			if (getRows() == 1 && getColumns() == 1)
				return *this;
			else {
				self T(getColumns(), getRows());
				for (long i = 0; i < T.getRows(); i++) {
					for (long j = 0; j < T.getColumns(); j++) {
						T[i][j] = __matrix[j][i];
					}
				}
				*this = T;
			}
			return *this;
		}// end transpose

		self operator+(const self& ob) const {
			if (getRows() == ob.getRows() && getColumns() == ob.getColumns()) {
				self A(getRows(), getColumns());
				for (long i = 0; i < getRows(); i++)
					for (long j = 0; j < getColumns(); j++)
						A[i][j] = __matrix[i][j] + ob[i][j];
				return A;
			}
			return self(0, 0);
		}

		self operator-(const self& ob) const {
			if (getRows() == ob.getRows() && getColumns() == ob.getColumns()) {
				self A(getRows(), getColumns());
				for (long i = 0; i < getRows(); i++)
					for (long j = 0; j < getColumns(); j++)
						A[i][j] = __matrix[i][j] - ob[i][j];
				return A;
			}
			return self(0, 0);
		}

		self& operator+=(const self& ob) {
			if (getRows() == ob.getRows() && getColumns() == ob.getColumns()) {
				for (long i = 0; i < getRows(); i++)
					for (long j = 0; j < getColumns(); j++)
						__matrix[i][j] += ob[i][j];
				return *this;
			}
			return Zero;
		}

		self& operator-=(const self& ob) {
			if (getRows() == ob.getRows() && getColumns() == ob.getColumns()) {
				for (long i = 0; i < getRows(); i++)
					for (long j = 0; j < getColumns(); j++)
						__matrix[i][j] -= ob[i][j];
				return *this;
			}
			return Zero;
		}

		template<typename F>
		self operator*(const F arg) const {
			self A(getRows(), getColumns());
			for (long i = 0; i < getRows(); i++)
				for (long j = 0; j < getColumns(); j++)
					A[i][j] = __matrix[i][j] * arg;
			return A;
		}

		template<typename F>
		self& operator*=(const F arg) {
			for (long i = 0; i < getRows(); i++)
				for (long j = 0; j < getColumns(); j++)
					__matrix[i][j] *= arg;
			return *this;
		}

		self operator*(const self& ob) const {
			if (getColumns() == ob.getRows()) {
				self A(getRows(), ob.getColumns());
				for (long i = 0; i < getRows(); i++) {
					for (long j = 0; j < ob.getColumns(); j++) {
						Type sum = Type{};
						for (long k = 0; k < getColumns(); k++)
							sum += __matrix[i][k] * ob[k][j];
						A[i][j] = sum;
					}
				}
				return A;
			}
			return Zero;
		}

		bool operator<=(const self& ob) const {
			bool p = true;

			if (getRows() == ob.getRows() && getColumns() == ob.getColumns()) {
				for (long i = 0; i < ob.getRows(); i++) {
					for (long j = 0; j < ob.getColumns(); j++)
						if (__matrix[i][j] > ob[i][j]) {
							p = false;
							break;
						}
					if (!p) break;
				}
			}
			else
				p = false;
			return p;
		}

		bool empty() const {
			if (__n == 0 || __m == 0)
				return true;
			return false;
		}

		long getRows() const {
			return __n;
		}// end getRows

		long getColumns() const {
			return __m;
		}// end getColumns

};

template <typename T>
Matrix<T> Matrix<T>::Zero = Matrix<T>(0,0);

template<typename T, typename F>
Matrix<T> operator*(const F arg, const Matrix<T> ob) {
	Matrix<T> A(ob.getRows(), ob.getColumns());
	for (long i = 0; i < ob.getRows(); i++)
		for (long j = 0; j < ob.getColumns(); j++)
			A[i][j] = ob[i][j] * arg;
	return A;
}

template <typename T>
Matrix<T> abs(const Matrix<T>& M) {
	Matrix<T> A(M.getRows(), M.getColumns());
	for (long i = 0; i < M.getRows(); i++)
		for (long j = 0; j < M.getColumns(); j++)
				A[i][j] = abs(M[i][j]);
	return A;
}

template <typename T>
T norm_max(const Matrix<T>& ob) {
	T max = T{};
	Matrix<T> A(ob.getRows(), 1, T{});
	for (long i = 0; i < ob.getRows(); i++)
		for (long j = 0; j < ob.getColumns(); j++)
			A[i][0] += abs(ob[i][j]);
	
	max = A[0][0];
	for (long i = 0; i < A.getRows(); i++)
		if (max < A[i][0])
			max = A[i][0];
	return max;
}// end norm_max

template<typename T>
bool is_square(const Matrix<T>& m) {
	return (!m.empty() && (m.getColumns() == m.getRows()));
}// end is_square

template<typename T>
typename Matrix<T>::Position partialSelection(const Matrix<T>& M){
	if (M.empty())
		throw *(new MatrixException("Matrix is empty"));

	Matrix<T>::Position P = {0, 0};
	T max = abs(M[0][0]);
	for (long i = 0; i < M.getRows(); ++i)
		if (max < abs(M[i][0])) {
			max = abs(M[i][0]);
			P = {i, 0};
		}
	return P;
}// end partialSelection

template<typename T>
typename Matrix<T>::Position fullSelection(const Matrix<T>& M) {
	if (M.empty())
		throw *(new MatrixException("Matrix is empty"));

	Matrix<T>::Position P = {0, 0};
	T max = abs(M[0][0]);
	for (long i = 0; i < M.getRows(); ++i)
		for (long j = 0; j < M.getColumns(); ++j)
			if (max < abs(M[i][j])) {
				max = abs(M[i][j]);
				P = { i, j };
			}
	return P;
}// end fullSelection

#define BASIC_GAUSSIAN				1
#define PARTIAL_SELECTION			2
#define SCALABLE_CHOICE				3
#define FULL_SELECTION				4

/* LU - decomposes matrix L to LU representation.
* Matrix L has ones on its diagonal (L[i][i] = 1 for each i = 0, 1, ..., n).
* @params:
*	Matrix<T>& L - lower matrix,
*	Matrix<T>& U - upper matrix. Its also input matrix - which will be decomposed,
*   int method - special flag specifies algortihm for LU decomposition.
* @return: void
* @throws: -
*/
template<typename T>
void LU(Matrix<T>& L, Matrix<T>& U, std::vector<long>& p, int method = BASIC_GAUSSIAN) {
	if (U.empty())
		throw *(new MatrixException("Matrix is empty"));

	//if(L.getRows() != U.getRows() || L.getColumns() != U.getColumns())
	//	throw *(new MatrixException("LU decomposition can not be done. Matrixes L and U are different size"));
	
	// Amount of matrix rows:
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


 /* LU - decomposes matrix L to LU representation.
 * Matrix L has ones on its diagonal (L[i][i] = 1 for each i = 0, 1, ..., n).
 * @params:
 *	Matrix<T>& L - lower matrix,
 *	Matrix<T>& U - upper matrix. Its also input matrix - which will be decomposed,
 *   int method - special flag specifies algortihm for LU decomposition.
 * @return: void
 * @throws: -
 */
template<typename T>
void LU2(Matrix<T>& A, std::vector<long>& p, int method = BASIC_GAUSSIAN) {
	if (A.empty())
		throw *(new MatrixException("Matrix is empty"));

	//if(L.getRows() != U.getRows() || L.getColumns() != U.getColumns())
	//	throw *(new MatrixException("LU decomposition can not be done. Matrixes L and U are different size"));

	// Amount of matrix rows:
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
}// end LU

template <typename T>
Matrix<T> Kronecker(const std::vector<long>& i, const std::vector<long>& j) {
	Matrix<T> delta(i.size(), j.size(), T{});
	for (long k = 0; k < i.size(); ++k)
		for (long l = 0; l < j.size(); ++l)
			if (i[k] == j[l])
				delta[k][l] = T{ 1 };
	return delta;
}// end Kronecker

template<typename T>
void print(const Matrix<T>& M, int width = 8, int precison = 4) {
	using namespace std;
	for (long i = 0; i < M.getRows(); ++i) {
		for (long j = 0; j < M.getColumns(); ++j) {
			cout.width(width);
			cout.precision(precison);
			cout << M[i][j] << " ";
		}
		cout << endl;
	}
	return;
}// end print

#endif /* __MATRIX */