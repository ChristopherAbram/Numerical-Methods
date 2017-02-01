#ifndef __MATRIX
#define __MATRIX

#include<vector>
#include<initializer_list>

class MatrixException : public std::runtime_error {
	public:
		explicit MatrixException(const std::string& what_str): std::runtime_error(what_str) {}
};

/** Matrix - special class to represent math matrix objects.
* @author: Christopher Abram
*/
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

		using value_type = Type;

		Matrix(std::initializer_list<std::vector<Type>> list) {
			std::initializer_list<std::vector<Type>>::iterator b = list.begin();
			std::initializer_list<std::vector<Type>>::const_iterator e = list.end();

			__matrix.reserve(list.size()); // reserves rows
			long max = b->size();
			while (b != e) {
				__matrix.push_back(*b);
				if (!b && max < b->size())
					max = b->size();
				++b;
			}
			__n = __matrix.size();
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
		}// end operator*

		std::vector<Type> operator*(const std::vector<Type>& ob) const {
			long n = ob.size();
			std::vector<Type> m(n);
			if (getColumns() == n) {
				for (long i = 0; i < getRows(); i++) {
					//for (long j = 0; j < ob.getColumns(); j++) {
					Type sum{};
					for (long k = 0; k < getColumns(); k++)
						sum += __matrix[i][k] * ob[k];
					m[i] = sum;
					//}
				}
			}
			else
				throw *(new MatrixException("Warning: Unable to multiply matrix by vector, inconsistent dimensions."));
			return m;
		}// end operator*

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

template<typename T>
std::vector<T> operator -(const std::vector<T>& l, const std::vector<T>& r) {
	std::vector<T> m(l.size());
	for (long i = 0; i < l.size() && i < r.size(); ++i)
		m[i] = l[i] - r[i];
	return m;
}// end operator-

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
	std::vector<T> A(ob.getRows());
	for (long i = 0; i < ob.getRows(); i++)
		for (long j = 0; j < ob.getColumns(); j++)
			A[i] += abs(ob[i][j]);
	
	max = A[0];
	for (long i = 0; i < A.size(); i++)
		if (max < A[i])
			max = A[i];
	return max;
}// end norm_max

template <typename T>
T norm_max(const std::vector<T>& ob) {
	T max = abs(ob[0]);
	for (long i = 0; i < ob.size(); i++)
		if (max < abs(ob[i]))
			max = abs(ob[i]);
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