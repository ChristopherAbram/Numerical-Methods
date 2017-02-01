#include<iostream>
#include<iomanip>
#include"ResultsStructures.h"

std::ostream& operator <<(std::ostream& os, const Result& r) {
	os.width(5);
	os << r.k;

	os.width(25);
	os << std::setprecision(15) << r.x;

	os.width(25);
	os << std::setprecision(15) << r.est;

	os.width(25);
	os << std::setprecision(15) << r.f;
	return os;
}