#ifndef __MATRIX
#define __MATRIX

#include<vector>
#include<initializer_list>

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

		self row(long n) const {
			self R = __matrix[n];
			return R;
		}// end row

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

		Type determinant() {}

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
}

template<typename T>
bool is_square(Matrix<T>& m) {
	bool p = false;

	return p;
}// end is_square

#endif /* __MATRIX */