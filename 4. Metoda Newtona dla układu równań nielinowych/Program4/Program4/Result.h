#ifndef __RESULT
#define __RESULT

#include "Matrix.h"

#define DELTA_REACHED	0x01
#define EPSILON_REACHED 0x02
#define ITER_REACHED	0x04
#define BAD_INPUT		0x08

struct Result {
	typedef double value_type;

	long n;
	// X <=> F(x) = 0:
	Matrix<value_type> X;
	// Values of functions fn for argument xn:
	Matrix<value_type> F;
	// Error estimations:
	Matrix<value_type> E;
	// Bit mask
	int end_condition;
	// Amount of iterations:
	int k;

	Result(long n, int end_condition, Matrix<value_type> X, Matrix<value_type> F, Matrix<value_type> E, int k) : n{ n }, end_condition{ end_condition }, X{ X }, F{ F }, E{ E }, k{k} {}
};

#endif /* __RESULT */