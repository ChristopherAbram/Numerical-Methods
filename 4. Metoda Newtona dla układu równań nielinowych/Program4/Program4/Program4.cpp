#include<iostream>
#include<windows.h>
#include<conio.h>
#include<iomanip>

#include "Matrix.h"
#include "Newton.h"

#define PI 3.14159265358979323846

/* Napisz program w jêzyku „C/C++”, realizuj¹cy metodê Newtona rozwi¹zywania uk³adu trzech 
algebraicznych równañ nieliniowych, i zastosuj ten program do przyk³adu z zadania 1. 
Przyjmij takie przybli¿enie pocz¹tkowe, aby uzyskaæ zbie¿noœæ metody. Zastosuj trzy niezale¿ne kryteria 
zakoñczenia iteracji. Zadbaj o to, aby wyprowadzaæ na konsolê wyniki poœrednie obliczeñ dla ka¿dej iteracji, 
tak aby mo¿liwe by³o obserwowanie zbie¿noœci kolejnych przybli¿eñ pierwiastków i porównanie liczby iteracji 
niezbêdnych do uzyskania rozwi¹zania o zadanej dok³adnoœci. 
Oblicz jak zmienia siê residuum uk³adu w trakcie iteracji.

(1)   xy - 2 = 0
(2)   y/2 - sin(pi/4 - z) = 0
(3)   x^2 + y^2 + z^2 - 4 = 0
*/

double f1(double, double, double);
double f2(double, double, double);
double f3(double, double, double);

double f1dx(double, double, double);
double f1dy(double, double, double);
double f1dz(double, double, double);

double f2dx(double, double, double);
double f2dy(double, double, double);
double f2dz(double, double, double);

double f3dx(double, double, double);
double f3dy(double, double, double);
double f3dz(double, double, double);

HANDLE  hConsole;

int main(int argc, char** argv) {

	using namespace std;
	using value_type = double;

	hConsole = GetStdHandle(STD_OUTPUT_HANDLE);

	Matrix<int> A = {
		{ 1, 3 },
		{ 2, 6 },
		{ 3, -7 },
		{ -4, 3 },
		{ -1, 1 },
	};

	cout << norm_max(A) << endl;

	Newton NEWTON(3);

	Newton::setDelta(0.000000000000001);
	Newton::setEpsilon(0.000000000000001);
	Newton::setMax(50);
	
	// Function Matrix:
	Matrix<math_function> F = {
		{ f1 },
		{ f2 },
		{ f3 },
	};

	// Reverse Jacobian Matrix:
	Matrix<math_function> J = {
		{ f1dx, f1dy, f1dz },
		{ f2dx, f2dy, f2dz },
		{ f3dx, f3dy, f3dz },
	};

	// Start point for Newton's method:
	Matrix<Result::value_type> X = {
		{ 1. },
		{ 1. },
		{ 1. },
	};

	// Setting NEWTON obj:
	NEWTON.setFunctions(F);
	NEWTON.setJacobian(J);
	NEWTON.setBeginPoint(X);

	// Let's solve at last this system:
	Result R = NEWTON.solve();

	vector<Result> app = NEWTON.getApproximations();
	int width = 17;

	cout << "----------------------------------------------------------------------------------------------------------------------------" << endl;
	cout.width(4);
	cout << "n";
	cout.width(width);
	cout << setprecision(10) << "x";
	cout.width(width);
	cout << setprecision(10) << "y";
	cout.width(width);
	cout << setprecision(10) << "z";
	cout.width(width);
	cout << setprecision(10) << "f1(X)";
	cout.width(width);
	cout << setprecision(10) << "f2(X)";
	cout.width(width);
	cout << setprecision(10) << "f3(X)";
	cout.width(width);
	cout << setprecision(10) << "|f(X)|";
	cout << endl;
	cout << "----------------------------------------------------------------------------------------------------------------------------" << endl;

	for (int i = 0; i < app.size(); i++) {
		cout.width(4);
		cout << app[i].k;
		cout.width(width);
		cout << setprecision(10) << app[i].X[0][0];
		cout.width(width);
		cout << setprecision(10) << app[i].X[1][0];
		cout.width(width);
		cout << setprecision(10) << app[i].X[2][0];
		cout.width(width);
		cout << setprecision(10) << app[i].F[0][0];
		cout.width(width);
		cout << setprecision(10) << app[i].F[1][0];
		cout.width(width);
		cout << setprecision(10) << app[i].F[2][0];
		cout.width(width);
		cout << setprecision(10) << norm_max(app[i].F);
		cout << endl;
	}

	cout << "----------------------------------------------------------------------------------------------------------------------------" << endl;
	SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN);
	cout.width(4);
	cout << R.k;
	cout.width(width);
	cout << setprecision(10) << R.X[0][0];
	cout.width(width);
	cout << setprecision(10) << R.X[1][0];
	cout.width(width);
	cout << setprecision(10) << R.X[2][0];
	cout.width(width);
	cout << setprecision(10) << R.F[0][0];
	cout.width(width);
	cout << setprecision(10) << R.F[1][0];
	cout.width(width);
	cout << setprecision(10) << R.F[2][0];
	cout.width(width);
	cout << setprecision(10) << norm_max(R.F);
	cout << endl;
	SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_RED);
	cout << "----------------------------------------------------------------------------------------------------------------------------" << endl;

	cout << "End condition: ";
	if (R.end_condition & MAX_ITER)
		cout << "Max amout of iteration ";
	if (R.end_condition & EPSILON_REACHED)
		cout << "EPSILON";
	if (R.end_condition & DELTA_REACHED)
		cout << "DELTA";

	cout << endl;

	_getch();
	return EXIT_SUCCESS;
}

double f1(double x, double y, double z) {return (x*y - 2);}
double f2(double x, double y, double z) { return (0.5*y - sin(PI / 4. - z)); }
double f3(double x, double y, double z) { return (x*x + y*y + z*z - 4.); }

double detJ(double x, double y, double z) { return (y*z + 2.*cos(PI/4. - z)*(x*x - y*y)); }

double f1dx(double x, double y, double z) { return ((z - 2*y*cos(PI/4. - z)) / detJ(x,y,z)); }
double f1dy(double x, double y, double z) { return ((-2.*x*z) / detJ(x,y,z)); }
double f1dz(double x, double y, double z) { return ((x*cos(PI / 4. - z)) / detJ(x,y,z)); }

double f2dx(double x, double y, double z) { return ((2.*x*cos(PI / 4. - z)) / detJ(x,y,z)); }
double f2dy(double x, double y, double z) { return ((2.*y*z) / detJ(x, y, z)); }
double f2dz(double x, double y, double z) { return ((-y * cos(PI / 4. - z)) / detJ(x,y,z)); }

double f3dx(double x, double y, double z) { return (-x / detJ(x,y,z)); }
double f3dy(double x, double y, double z) { return ((2.*x*x - 2.*y*y) / detJ(x,y,z)); }
double f3dz(double x, double y, double z) { return (0.5 * y / detJ(x,y,z)); }