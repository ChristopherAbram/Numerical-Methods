#ifndef __RESULT
#define __RESULT

#include<vector>

#define DELTA_REACHED	0x01
#define EPSILON_REACHED 0x02
#define ITER_REACHED	0x04
#define BAD_INPUT		0x08

template<typename T>
struct Result {
	// X <=> F(x) = 0:
	std::vector<T> X;
	// | Ax - b |
	T RES;
	// Error estimation:
	T EST;
	// Bit mask
	int end_condition;
	// Amount of iterations:
	int k;

	Result(int k = 0, int end_condition = 0, std::vector<T>& X = std::vector<T>(), T RES = T{}, T EST = T{}) : end_condition{ end_condition }, X{ X }, RES{ RES }, EST{ EST }, k{ k } {}
};

template<typename T>
std::ostream& operator <<(std::ostream& os, const Result<T>& r) {
	os.width(5);
	os << r.k;

	for (int i = 0; i < r.X.size(); ++i) {
		os.width(15);
		os << std::setprecision(12) << r.X[i];
	}

	os.width(15);
	os << std::setprecision(8) << r.EST;

	os.width(15);
	os << std::setprecision(8) << r.RES;
	return os;
}// end operator <<

#endif /* __RESULT */