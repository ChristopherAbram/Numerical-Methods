#include<iostream>
#include<cfloat>
#include<conio.h>

/* Zaimplementuj w języku „C/C++” algorytm obliczający przybliżone wartości funkcji f(x) = exp(x) dla x  [-30, 30], 
poprzez sumowanie N wyrazów rozwinięcia tej funkcji w szereg Taylora wokół x = 0. 
Zbadaj jak zmieniają się błędy względne przybliżenia funkcji w tym algorytmie, przy wzrastającej liczbie N  [1, 1000]. 
Wyjaśnij przyczyny obserwowanych błędów i ich zmian ze wzrostem N. 
Następnie dokonaj takiej modyfikacji algorytmu, i wybierz takie N, aby uzyskać dokładność maszynową dla dowolnego x  [-30, 30]. 
W obliczeniach zastosuj zmienne podwójnej precyzji.
*/

const int N_MAX = 1000;
const int N_MIN = 1;

double exp_first(double, int);
double exp_second(double, int);
double blad_wzgl(double, double);

int main(int argc, char** argv){
	using namespace std;

	double x, fx1, fx2;

	double arr[6] = { -30, -20, -10, 10, 20, 30 };

	// Pierwsza metoda, badanie zależności błędu względnego od ilości iteracji N
	cout << "Pierwsza metoda, badanie zaleznosci bledu wzglednego od ilosci iteracji N:" << endl << endl;
	cout << "  x      N       f(x)           Blad wzgledny  " << endl;
	for (int j = 0; j < 6; j++) {
		for (int k = 10; k <= N_MAX/10; k += 15) {
			fx1 = exp_first(arr[j], k);
			cout.width(5);
			cout << arr[j];
			cout.width(7);
			cout << k;
			cout.width(15);
			cout.precision(7);
			cout << fx1;
			cout.width(15);
			cout.precision(4);
			cout << blad_wzgl(arr[j], fx1) << endl;
		}
		cout << endl;
	}

	cout << endl << endl;

	// Druga metoda, badanie zależności błędu względnego od ilości iteracji N
	cout << "Druga metoda, badanie zaleznosci bledu wzglednego od ilosci iteracji N:" << endl << endl;
	cout << "  x      N       f(x)           Blad wzgledny  " << endl;
	for (int j = 0; j < 6; j++) {
		for (int k = 10; k <= N_MAX/10; k += 15) {
			fx2 = exp_second(arr[j], k);
			cout.width(5);
			cout << arr[j];
			cout.width(7);
			cout << k;
			cout.width(15);
			cout.precision(7);
			cout << fx2;
			cout.width(15);
			cout.precision(4);
			cout << blad_wzgl(arr[j], fx2) << endl;
		}
		cout << endl;
	}

	cout << endl << endl;

	cout << "Wartosci liczone druga metoda:" << endl << endl;
	cout << "  x       f(x)                     Blad wzgledny  " << endl;
	for (x = -30.L; x <= 30.L; x += 1.L) {
		fx2 = exp_second(x, N_MAX);

		cout.width(5);
		cout << x;
		cout.width(25);
		cout.precision(15);
		cout << fx2;
		cout.width(15);
		cout.precision(4);
		cout << blad_wzgl(x, fx2) << endl;
	}

	cout << endl;

	cout << "Nacisnij dowolny klawisz...";
	_getch();
	return EXIT_SUCCESS;
}// end main

double exp_first(double x, int c) {
	double r = 1.L;
	int i = 0;
	double k = 1.L;
	do {
		k *= x / ((double)i + 1.L);
		r += k;
		i++;
	} while (i <= c);

	return r;
}// end exp_first

double exp_second(double x, int c) {
	double r = 1.L;
	int i = 0;
	double k = 1.L;
	bool p = x >= 0;
	x = abs(x);
	do {
		k *= x / ((double)i + 1.L);
		r += k;
		i++;
	} while (i <= c);

	return p ? r : 1.L / r;
}// end exp_second

double blad_wzgl(double x, double fx) {
	double f = exp(x);
	return abs((fx - f) / f);
}// end blad_wzgl