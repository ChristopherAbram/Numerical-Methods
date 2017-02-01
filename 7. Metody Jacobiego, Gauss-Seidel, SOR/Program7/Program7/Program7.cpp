#include<iostream>
#include<conio.h>
#include<windows.h>
#include<iomanip>

#include "SolvingLinearSystem.h"

/**
Napisz program w jêzyku „C/C++”, rozwi¹zuj¹cy uk³ad czterech równañ liniowych metodami
iteracyjnymi: 
	(a) Jacobiego, 
	(b) Gaussa-Seidela, 
	(c) SOR z parametrem w = 1/2, 
a nastêpnie zastosuj ten program do rozwi¹zania uk³adu równañ liniowych Ax = b, gdzie:

macierz A:

	100 -1  2   -3
	1   200 -4  5
	-2  4   300 -6
	3   -5  6   400

macierz b:

	116
	-226
	912
	-1174

Przyjmij przybli¿enie pocz¹tkowe x0:

	2
	2
	2
	2

Zastosuj trzy niezale¿ne kryteria zakoñczenia iteracji. Zadbaj o to, aby wyprowadzaæ na konsolê
wyniki poœrednie obliczeñ dla ka¿dej iteracji, tak aby mo¿liwe by³o obserwowanie zbie¿noœci
kolejnych przybli¿eñ pierwiastków i porównanie liczby iteracji niezbêdnych do uzyskania
rozwi¹zania o zadanej dok³adnoœci bezwzglêdnej. Oblicz jak zmienia siê residuum uk³adu w trakcie
kolejnych iteracji.
*/

HANDLE  hConsole;

template<typename T>
void print(std::vector<Result<T>>& r) {
	using namespace std;
	
	cout << "------------------------------------------------------------------------------------------------" << endl;
	cout.width(5);
	cout << "n";
	if(r.size() != 0)
		for (long i = 0; i < r[0].X.size(); ++i) {
			cout.width(14);
			cout << setprecision(10) << "x" << i;
		}
	cout.width(15);
	cout << setprecision(10) << "||est||";
	cout.width(15);
	cout << setprecision(10) << "||Rn||";
	cout << endl;
	cout << "------------------------------------------------------------------------------------------------" << endl;

	/*if (r[0].end_condition & BAD_INPUT) {
		SetConsoleTextAttribute(hConsole, FOREGROUND_RED);
		cout << endl << "No itaration done, bad input data: function, derivative or range." << endl << endl;
		SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_RED);
	}
	else {*/
		vector<Result<T>>::iterator b = r.begin();
		vector<Result<T>>::const_iterator e = r.end();
		while (b != e) {
			cout << *b << endl;
			++b;
		}
	//}

	cout << "------------------------------------------------------------------------------------------------" << endl;

	if (r.empty() || r.back().end_condition & BAD_INPUT) {
		SetConsoleTextAttribute(hConsole, FOREGROUND_RED);
		cout << "No results..." << endl;
		SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_RED);
	}
	else {
		SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN);
		cout << r.back() << endl;
		SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_RED);
	}

	cout << "------------------------------------------------------------------------------------------------" << endl;
	cout << "End condition: ";
	if (r.empty() || r.back().end_condition & BAD_INPUT) {
		SetConsoleTextAttribute(hConsole, FOREGROUND_RED);
		cout << "BAD INPUT DATA" << endl;
		SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_RED);
	}
	else if (!r.empty() && r.back().end_condition == 0)
		cout << "NONE" << endl;
	else if(!r.empty()) {
		if (r.back().end_condition & ITER_REACHED)
			cout << "MAX AMOUT OF ITERATION, ";
		if (r.back().end_condition & EPSILON_REACHED)
			cout << "EPSILON, ";
		if (r.back().end_condition & DELTA_REACHED)
			cout << "DELTA";
		cout << endl;
	}

	return;
}// end print

int main(int argc, char** argv) {

	using namespace std;
	using type = double;

	hConsole = GetStdHandle(STD_OUTPUT_HANDLE);

	// Initialization {

		/*Matrix<type> A = {
			{100., -1., 2., -3.},
			{1., 200., -4., 5.},
			{-2., 4., 300., -6.},
			{3., -5., 6., 400.}
		};
		Matrix<type> tmpA = A;

		vector<type> b = {
			116.,
			-226.,
			912.,
			-1174.
		};
		vector<type> tmpb = b;

		vector<type> x0 = {
			2.,
			2.,
			2.,
			2.
		};
		vector<type> tmpx = x0;*/

		Matrix<type> A = {
			{ 1., 2., 3. },
			{ 2., 2., 3. },
			{ 3., 3., 3. }
		};
		Matrix<type> tmpA = A;

		vector<type> b = {
			0.,
			0.,
			0.
		};
		vector<type> tmpb = b;

		vector<type> x0 = {
			1., 1., 1.
		};
		vector<type> tmpx = x0;

	// } Solving linear system {

		SLS<type>::setDelta(1e-15);
		SLS<type>::setEpsilon(1e-15);
		SLS<type>::setMax(100);

		try {
			SLS<type> solver(tmpA, tmpb, tmpx);
			solver.areErrorsAvailable(true); // enables error storing

			vector<Result<type>> results;

			// solving using Jacob's method {
				
				solver.solveJacobi();
				results = solver.getResults();

				cout << endl << "Results for Jacobi's method: " << endl;
				print(results);

			// } solving using Gauss-Seidel's method {

				tmpA = A;
				tmpb = b;
				tmpx = x0;
				solver.setA(tmpA).setb(tmpb).setBeginningPoint(tmpx);

				solver.solveGaussSeidel();
				results = solver.getResults();

				cout << endl << endl << "Results for Gauss-Seidel's method: " << endl;
				print(results);

			// } solving using SOR method {

				tmpA = A;
				tmpb = b;
				tmpx = x0;
				solver.setA(tmpA).setb(tmpb).setBeginningPoint(tmpx);

				solver.solveSOR(0.5);
				results = solver.getResults();

				cout << endl << endl << "Results for SOR method: " << endl;
				print(results);

			// }

		}
		catch (SLSException& slse) {
			cout << slse.what() << endl;
		}
		catch (std::exception& e) {
			cout << e.what() << endl;
		}

	// }


	cout << endl << endl << "Press any button to finish...";
	_getch();
	return EXIT_SUCCESS;
}