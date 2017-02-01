#ifndef __RESULTSSTRUCTURES
#define __RESULTSSTRUCTURES

#include<iostream>
#include<iomanip>

#define DELTA_REACHED	0x01
#define EPSILON_REACHED 0x02
#define ITER_REACHED	0x04
#define BAD_INPUT		0x08

struct Result {
	typedef double value_type;

	// x <=> f(x) = 0:
	value_type x;
	// Value of function f for argument x:
	value_type f;
	// Error estimation:
	value_type est;
	// Amount of iterations:
	int k;
	// Bit mask
	int end_condition;

	Result(int end = 0, value_type x = 0., value_type f = 0., value_type e = 0., int k = 0) : x{ x }, f{ f }, est{ e }, k{k}, end_condition { end } {}
	friend std::ostream& operator <<(std::ostream& os, const Result& r);
};

std::ostream& operator <<(std::ostream& os, const Result& r);

#endif /* __RESULTSSTRUCTURES */