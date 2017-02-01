#ifndef __NONLINEAREQUATIONSOLVING
#define __NONLINEAREQUATIONSOLVING

#include<cmath>
#include<vector>
#include"ResultsStructures.h"

#define DELTA		0.00001
#define EPSILON		0.00001
#define MAX_ITER	100

typedef Result::value_type(*math_function)(Result::value_type);

/* NonLinearEquationSolving - class contains methods that enable you 
* to calculate zero of function f(x) = 0 in given range [a,b]. 
* @author: Christopher Abram
*/
class NLES { 
	private:

		// Iteration parameters:
		static Result::value_type __delta;    // TOL
		static Result::value_type __epsilon;  // TOLF
		static int __M;			  // max amout of iteraions 

		// Function pointers:
		static math_function __f;
		static math_function __d;

		// Seeking range:
		static Result::value_type __a; // left side end
		static Result::value_type __b; // rigth side end

		// The series of successive approximations:
		static std::vector<Result> __appromax;

	public: 

		/* bisection - calculates zero of function f using bisection algorithm.
		* @param -
		* @return Result r - struct with results
		*/
		static Result bisection();

		/* picard - calculates zero of function f using Picard's method
		* @param Result::value_type x - start point.
		* @return Result r - struct with results
		*/
		static Result picard(Result::value_type x = __b);

		/* newton - caclulates zero of function f using Newton's method.
		* @param Result::value_type x - start point.
		* @return Result r - struct with results.
		*/
		static Result newton(Result::value_type x = __b);

		/* secant - caclulates zero of function f using secant method.
		* @param Result::value_type x - start point.
		* @return Result r - struct with results.
		*/
		static Result secant(Result::value_type x1, Result::value_type x2);

		 /* bestPoint - choose point x in range [a, b] for which f(x) is min.
		 * @param Result::value_type step - step value for scanning range [a, b]
		 * @return Result::value_type x - the best point
		 */
		static Result::value_type bestPoint(Result::value_type step);

		static std::vector<Result> getApproximations();
		static void setRange(Result::value_type a, Result::value_type b);
		static int signum(Result::value_type x);
		static void setFunction(math_function f);
		static math_function getFunction();
		static void setDerivative(math_function d);
		static math_function getDerivative();
		static void setMax(int max);
		static int getMax();
		static void setDelta(Result::value_type d);
		static void setEpsilon(Result::value_type e);
		static Result::value_type getDelta();
		static Result::value_type getEpsilon();
};

#endif /* __NONLINEAREQUATIONSOLVING */