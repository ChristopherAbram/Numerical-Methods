#include<iostream>
#include<conio.h>

#include "SolvingLinearSystem.h"

/*
Dodatkowe æwiczenie, dla chêtnych; prawid³owe rozwi¹zanie bêdzie premiowane dodatkow¹ ocen¹ 5.0.
Dany jest uk³ad równañ z macierz¹ A:

	1 + e     1       1        1
	  1     1 + e     1        1
	  1       1     1 + e      1
	  1       1       1      1 + e

oraz wektorem b:

	6 + e
	6 + 2e
	6 + 2e
	6 + e

. ZnajdŸ najpierw rozwi¹zanie analityczne, a nastêpnie numeryczne, przyjmuj¹c coraz mniejsze wartoœci e = 10e-5, 10e-6, 10e-7, itd.
Porównaj rozwi¹zania numeryczne z analitycznym, i wyjaœnij ewentualne obserwowane zmiany b³êdu rozwi¹zañ numerycznych.*/

int main(int argc, char** argv) {

	using namespace std;
	using type = double;

	Matrix<type> A(4, 4);
	vector<type> b(4);

	type e = 10.e-5;

	for (int k = 0; k < 15; ++k) {

		/* Matrix A */
		A = {
			{ 1. + e, 1., 1., 1. },
			{ 1., 1. + e, 1., 1. },
			{ 1., 1., 1. + e, 1. },
			{ 1., 1., 1., 1. + e }
		};

		b = {
			6. + e,
			6. + 2. * e,
			6. + 2. * e,
			6. + e
		};

		try {
			/* Solving linear system */
			SLS<type> solver(A, b);
			b = solver.solveGaussianElimination().getSolutions();
		}
		catch (SLSException& slse) {
			cout << slse.what() << endl;
		}
		catch (MatrixException& me) {
			cout << me.what() << endl;
		}

		/* Presentation of results */
		cout << endl << "Solutions of the system Ax = b when e = " << scientific << e << defaultfloat << endl;
		for (long i = 0; i < b.size(); ++i) {
			cout.precision(19);
			cout << " x" << i + 1 << " = " << b[i] << endl;
		}
		cout << endl;

		e /= 10.;
	}

	std::cout << endl << endl << "Press any button to close...";
	_getch();
	return EXIT_SUCCESS;
}// end main