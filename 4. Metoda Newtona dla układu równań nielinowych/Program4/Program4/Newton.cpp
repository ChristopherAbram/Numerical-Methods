#include "Newton.h"

Result::value_type Newton::__delta = DELTA;
Result::value_type Newton::__epsilon = EPSILON;
int Newton::__M = MAX_ITER;

std::vector<Result> Newton::getApproximations() {
	return __appromax;
}// end getApproximations

void Newton::setMax(int max) {
	__M = max;
	return;
}// end setMax

int Newton::getMax() {
	return __M;
}// end getMax

void Newton::setDelta(Result::value_type d) {
	__delta = d;
	return;
}// end setDelta

void Newton::setEpsilon(Result::value_type e) {
	__epsilon = e;
	return;
}// end setEpsilon

Result::value_type Newton::getDelta() {
	return __delta;
}// end getDelta

Result::value_type Newton::getEpsilon() {
	return __epsilon;
}// end getEpsilon