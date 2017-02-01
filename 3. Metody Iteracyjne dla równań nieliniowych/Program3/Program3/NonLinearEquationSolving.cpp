#include<cmath>
#include<vector>
#include"ResultsStructures.h"
#include"NonLinearEquationSolving.h"

static Result::value_type f(Result::value_type x) { return x; }
static Result::value_type d(Result::value_type x) { return 1.; }

Result::value_type NLES::__delta = DELTA;
Result::value_type NLES::__epsilon = EPSILON;
int NLES::__M = MAX_ITER;

math_function NLES::__f = f;
math_function NLES::__d = d;

Result::value_type NLES::__a = 0.;
Result::value_type NLES::__b = 1.;

std::vector<Result> NLES::__appromax = std::vector<Result>(MAX_ITER);

Result NLES::bisection() {
	Result::value_type u, v, e, c, w, a, b;

	// Start values:
	b = __b;
	a = __a;
	u = getFunction()(a);
	v = getFunction()(b);
	e = b - a;

	if (signum(u) == signum(v))
		return Result(BAD_INPUT);

	__appromax.clear();
	__appromax.reserve(getMax());

	for (int k = 0; k < getMax(); k++) {
		e /= 2;
		c = a + e;
		w = getFunction()(c);

		// Adding next point to an array:
		__appromax.push_back(Result(0, c, w, e, k + 1));

		if ((abs(e) <= getDelta()) || (abs(w) <= getEpsilon()))
			return Result((((abs(e) <= getDelta()) ? DELTA_REACHED : 0) | ((abs(w) <= getEpsilon()) ? EPSILON_REACHED : 0)), c, w, e, k + 1);

		if (signum(w) != signum(u)) {
			b = c;
			v = w;
		}
		else {
			a = c;
			u = w;
		}
	}
	return Result(ITER_REACHED, c, w, e, getMax());
}// end bisection

Result NLES::picard(Result::value_type x) {
	Result::value_type x1, w;

	__appromax.clear();
	__appromax.reserve(getMax());

	w = getFunction()(x);
	if (abs(w) <= getEpsilon())
		return Result(EPSILON_REACHED, x, w);

	for (int k = 0; k < getMax(); k++) {
		w = getFunction()(x) + x;
		x1 = w;

		// Adding next point to an array:
		__appromax.push_back(Result(0, x, w - x, x1 - x, k + 1));

		if (abs(x1 - x) <= getDelta() || abs(w - x) <= getEpsilon())
			return Result((((abs(x1 - x) <= getDelta()) ? DELTA_REACHED : 0) | ((abs(w - x) <= getEpsilon()) ? EPSILON_REACHED : 0)), x, w - x, x1 - x, k + 1);
		x = w;
	}
	return Result(ITER_REACHED, x, w - x, x1 - x, getMax());
}// end picard

Result NLES::newton(Result::value_type x) {
	Result::value_type v, x1, d;

	if (signum(getFunction()(__a)) == signum(getFunction()(__b)))
		return Result(BAD_INPUT);

	__appromax.clear();
	__appromax.reserve(getMax());

	v = getFunction()(x);
	if (abs(v) <= getEpsilon())
		return Result(EPSILON_REACHED, x, v);

	for (int k = 0; k < getMax(); k++) {
		d = getDerivative()(x);
		if (d == 0.) {
			return Result(BAD_INPUT);
			//x += getDelta();
			//v = getFunction()(x);
			//continue;
		}
		else
			x1 = x - v / d;
		v = getFunction()(x1);

		// Adding next point to an array:
		__appromax.push_back(Result(0, x1, v, x1 - x, k + 1));

		if (abs(x1 - x) <= getDelta() || abs(v) <= getEpsilon())
			return Result((((abs(x1 - x) <= getDelta()) ? DELTA_REACHED : 0) | ((abs(v) <= getEpsilon()) ? EPSILON_REACHED : 0)), x1, v, x1 - x, k + 1);
		x = x1;
	}
	return Result(ITER_REACHED, x1, v, x1 - x, getMax());
}// end newton

Result NLES::secant(Result::value_type x1, Result::value_type x2) {
	Result::value_type u, v, c, w;

	if (signum(getFunction()(__a)) == signum(getFunction()(__b)))
		return Result(BAD_INPUT);

	__appromax.clear();
	__appromax.reserve(getMax());

	u = getFunction()(x1);
	if (abs(u) <= getEpsilon())
		return Result(EPSILON_REACHED, x1, u);
	v = getFunction()(x2);
	if (abs(v) <= getEpsilon())
		return Result(EPSILON_REACHED, x2, v);

	for (int k = 0; k < getMax(); k++) {
		if (abs(u - v) > getDelta())
			c = x1 - u * ((x1 - x2) / (u - v));
		else
			return Result(BAD_INPUT);
		w = getFunction()(c);

		// Adding next point to an array:
		__appromax.push_back(Result(0, c, w, c - x1, k + 1));

		if (abs(c - x1) <= getDelta() || abs(w) <= getEpsilon())
			return Result((((abs(x1 - c) <= getDelta()) ? DELTA_REACHED : 0) | ((abs(w) <= getEpsilon()) ? EPSILON_REACHED : 0)), c, w, c - x1, k + 1);
		x2 = x1;
		v = u;
		x1 = c;
		u = w;
	}
	return Result(ITER_REACHED, c, w, c - x1, getMax());
}// end secant

Result::value_type NLES::bestPoint(Result::value_type step) {
	Result::value_type a, b, min;
	a = __a;
	b = __b;
	min = a;
	for (; a <= b; a += step) {
		if (abs(getFunction()(min)) > abs(getFunction()(a)))
			min = a;
	}
	return min;
}// end bestPoint

std::vector<Result> NLES::getApproximations() {
	return __appromax;
}// end getApproximations

void NLES::setRange(Result::value_type a, Result::value_type b) {
	if (a > b) std::swap(a, b);
	__a = a;
	__b = b;
	return;
}// end setRange

int NLES::signum(Result::value_type x) {
	return x >= 0 ? 1 : -1;
}// end signum

void NLES::setFunction(math_function f) {
	__f = f;
	return;
}// end setFunction

math_function NLES::getFunction() {
	return __f;
}// getFunction

void NLES::setDerivative(math_function d) {
	__d = d;
	return;
}// end setDerivative

math_function NLES::getDerivative() {
	return __d;
}// getDerivative

void NLES::setMax(int max) {
	__M = max;
	return;
}// end setMax

int NLES::getMax() {
	return __M;
}// end getMax

void NLES::setDelta(Result::value_type d) {
	__delta = d;
	return;
}// end setDelta

void NLES::setEpsilon(Result::value_type e) {
	__epsilon = e;
	return;
}// end setEpsilon

Result::value_type NLES::getDelta() {
	return __delta;
}// end getDelta

Result::value_type NLES::getEpsilon() {
	return __epsilon;
}// end getEpsilon