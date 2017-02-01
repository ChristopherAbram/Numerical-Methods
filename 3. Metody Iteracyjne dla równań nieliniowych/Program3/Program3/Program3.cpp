#include<iostream>
#include<cmath>
#include<windows.h>
#include<conio.h>

#include "NonLinearEquationSolving.h"

/* Napisz program w jêzyku „C/C++”, realizuj¹cy metody:
(a) Picarda
(b) bisekcji
(c) Newtona
(d) siecznych
rozwi¹zywania pojedynczych algebraicznych równañ nieliniowych. Zastosuj program do
przyk³adów z zadania 1. Zastosuj trzy niezale¿ne kryteria zakoñczenia iteracji. Zadbaj o to, aby
wyprowadzaæ na konsolê wyniki poœrednie obliczeñ dla ka¿dej iteracji, tak aby mo¿liwe by³o
obserwowanie zbie¿noœci kolejnych przybli¿eñ pierwiastków i porównanie liczby iteracji
niezbêdnych do uzyskania rozwi¹zania o zadanej dok³adnoœci przez ka¿d¹ z metod. 

1) cos^2(x/4) - x = 0				-> f1(x) = cos^2(x/4) - x
2) exp(x) - exp(-x) + x - 1 = 0		-> f2(x) = exp(x) - exp(-x) + x - 1
*/

double f1(double);
double f2(double);
double df1(double);
double df2(double);

void displayTable(Result&, std::vector<Result>&);
HANDLE  hConsole;

int main(int argc, char** argv) {

	using namespace std;

	vector<Result> app;
	Result r;

	hConsole = GetStdHandle(STD_OUTPUT_HANDLE);

	// Accuracy of caculations:
	NLES::setEpsilon(0.000000000000001);
	NLES::setDelta(0.000000000000001);
	NLES::setMax(50);

	cout << endl << "Function f1(x) = 0" << endl;
	cout << "_________________________________________________________________________________" << endl << endl;
	
	// Function f1:
	NLES::setFunction(f1);
	NLES::setDerivative(df1);

	// Caculating with bisection:
	NLES::setRange(0., 2.); // [0, 2]
	r = NLES::bisection();
	app = NLES::getApproximations();

	cout << endl;
	cout << "Bisection method" << endl;
	displayTable(r, app);

	// Calculating with Picard:
	r = NLES::picard();
	app = NLES::getApproximations();

	cout << endl;
	cout << "Picard's method" << endl;
	displayTable(r, app);

	// Calculating with Newton:
	r = NLES::newton();
	app = NLES::getApproximations();

	cout << endl;
	cout << "Newton's method" << endl;
	displayTable(r, app);

	// Calculating with secant:
	r = NLES::secant(1.9, 2.);
	app = NLES::getApproximations();

	cout << endl;
	cout << "Secant method" << endl;
	displayTable(r, app);



	cout << endl << endl << "Function f2(x) = 0" << endl;
	cout << "_________________________________________________________________________________" << endl << endl;

	// Function f2:
	NLES::setFunction(f2);
	NLES::setDerivative(df2);

	// Caculating with bisection:
	NLES::setRange(-1., 3.); // [-1, 3]
	r = NLES::bisection();
	app = NLES::getApproximations();

	cout << endl;
	cout << "Bisection method" << endl;
	displayTable(r, app);

	// Calculating with Picard:
	r = NLES::picard();
	app = NLES::getApproximations();

	cout << endl;
	cout << "Picard's method" << endl;
	displayTable(r, app);

	// Calculating with Newton:
	r = NLES::newton();
	app = NLES::getApproximations();

	cout << endl;
	cout << "Newton's method" << endl;
	displayTable(r, app);

	// Calculating with secant:
	r = NLES::secant(2.9, 3.);
	app = NLES::getApproximations();

	cout << endl;
	cout << "Secant method" << endl;
	displayTable(r, app);

	_getch();
	return EXIT_SUCCESS;
}// end main

double f1(double x) {
	double tmp = cos(x / 4.);
	return (tmp * tmp - x);
}// end f1

double f2(double x) {
	return (exp(x) - exp(-x) + x - 1.);
}// end f1

double df1(double x) {
	return (-0.5 * sin(x / 4.) * cos(x / 4.) - 1.);
}// end df1

double df2(double x) {
	return (exp(x) + exp(-x) + 1.);
}// end f1

void displayTable(Result& r, std::vector<Result>& app) {
	using namespace std;

	cout << "---------------------------------------------------------------------------------" << endl;
	cout.width(5);
	cout << "n";
	cout.width(25);
	cout << setprecision(15) << "x";
	cout.width(25);
	cout << setprecision(15) << "est";
	cout.width(25);
	cout << setprecision(15) << "f(x)";
	cout << endl;
	cout << "---------------------------------------------------------------------------------" << endl;

	if (r.end_condition & BAD_INPUT) {
		SetConsoleTextAttribute(hConsole, FOREGROUND_RED);
		cout << endl << "No itaration done, bad input data: function, derivative or range." << endl << endl;
		SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_RED);
	}
	else {
		vector<Result>::iterator b = app.begin();
		vector<Result>::const_iterator e = app.end();
		while (b != e) {
			cout << *b << endl;
			++b;
		}
	}

	cout << "---------------------------------------------------------------------------------" << endl;

	if (r.end_condition & BAD_INPUT) {
		SetConsoleTextAttribute(hConsole, FOREGROUND_RED);
		cout << "No results..." << endl;
		SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_RED);
	}
	else {
		SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN);
		cout << r << endl;
		SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_RED);
	}

	cout << "---------------------------------------------------------------------------------" << endl;
	cout << "End condition: ";
	if (r.end_condition & BAD_INPUT) {
		SetConsoleTextAttribute(hConsole, FOREGROUND_RED);
		cout << "BAD INPUT DATA" << endl;
		SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_RED);
	}
	else if (r.end_condition == 0)
		cout << "NONE" << endl;
	else {
		if (r.end_condition & ITER_REACHED)
			cout << "MAX AMOUT OF ITERATION, ";
		if (r.end_condition & EPSILON_REACHED)
			cout << "EPSILON, ";
		if (r.end_condition & DELTA_REACHED)
			cout << "DELTA";
		cout << endl;
	}
}