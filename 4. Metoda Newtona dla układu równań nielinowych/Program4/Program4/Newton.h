#ifndef __NEWTON
#define __NEWTON

#include<cmath>
#include<vector>
#include<iostream>
#include "Matrix.h"
#include "Result.h"

#define DELTA		0.00001
#define EPSILON		0.00001
#define MAX_ITER	100

typedef Result::value_type(*math_function)(Result::value_type, Result::value_type, Result::value_type);

class Newton {
	private:

		int n;
		// Matrix of functions:
		Matrix<math_function> __F;
		// Reversed Jacobian matrix:
		Matrix<math_function> __J_R;
		// Begin point:
		Matrix<Result::value_type> __X;

		// Iteration parameters:
		static Result::value_type __delta;    // TOL
		static Result::value_type __epsilon;  // TOLF
		static int __M;			  // max amout of iteraions 

		std::vector<Result> __appromax;

	public:
		
		Newton(int n) : n{n}{}

		Result solve() {
			Matrix<Result::value_type> Zero(n, 1, Result::value_type{});
			Result R(n, 0, Zero, Zero, Zero, 0);

			__appromax.clear();

			Matrix<Result::value_type> X(n, 1), X1(n, 1), J(n, n), F(n, 1);
			Matrix<Result::value_type> E(n, 1, getEpsilon()), E_TMP(n, 1), D(n, 1, getDelta());

			X = __X; // beginig point
			F = calculateFunctions(X); // values of functions for begining point
			J = calculateJacobian(X);  // values of derivations for begining point

			for (int k = 0; k < getMax(); k++) {

				X1 = X - J * F;

				F = calculateFunctions(X1);
				J = calculateJacobian(X1);

				//std::cout << norm_max(F) << "    " << norm_max(E_TMP) << std::endl;

				E_TMP = X1 - X;
				__appromax.push_back(Result(n, 0, X1, F, E_TMP, k + 1));

				if (norm_max(E_TMP) <= getDelta() || norm_max(F) <= getEpsilon())
					return Result(n, ((norm_max(E_TMP) <= getDelta() ? EPSILON_REACHED : 0) | (norm_max(F) <= getEpsilon() ? DELTA_REACHED : 0)), X1, F, E_TMP, k + 1);

				X = X1;
			}

			return Result(n, MAX_ITER, X1, F, E_TMP, getMax());
		}// end solve

		Matrix<Result::value_type> calculateFunctions(Matrix<Result::value_type> X) {
			Matrix<Result::value_type> FX(n, 1, Result::value_type{});
			for (long i = 0; i < FX.getRows(); i++)
				FX[i][0] = __F[i][0](X[0][0], X[1][0], X[2][0]);
			
			return FX;
		}// end calculateFunctions

		Matrix<Result::value_type> calculateJacobian(Matrix<Result::value_type> X) {
			Matrix<Result::value_type> JX(n, n, Result::value_type{});
			for (long i = 0; i < JX.getRows(); i++)
				for (long j = 0; j < JX.getColumns(); j++)
					JX[i][j] = __J_R[i][j](X[0][0], X[1][0], X[2][0]);
			return JX;
		}// end calculateJacobian

		void setBeginPoint(const Matrix<Result::value_type>& X) {
			__X = X;
			return;
		}// end setBeginPoint

		Matrix<Result::value_type> getBeginPoint() {
			return __X;
		}// end getBeginPoint

		void setFunctions(const Matrix<math_function>& F) {
			__F = F;
			return;
		}// end setFunctions

		Matrix<math_function> getFunctions() {
			return __F;
		}// end getFunctions

		void setJacobian(const Matrix<math_function>& J) {
			__J_R = J;
			return;
		}// end setFunctions

		Matrix<math_function> getJacobian() {
			return __J_R;
		}// end getJacobian

		std::vector<Result> getApproximations();
		static void setMax(int max);
		static int getMax();
		static void setDelta(Result::value_type d);
		static void setEpsilon(Result::value_type e);
		static Result::value_type getDelta();
		static Result::value_type getEpsilon();
};

#endif /* __NEWTON */