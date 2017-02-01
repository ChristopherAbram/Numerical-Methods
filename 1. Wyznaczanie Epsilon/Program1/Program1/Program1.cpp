#include<iostream>
#include <conio.h>
#include<cfloat>

/*
Napisz program w jêzyku „C/C++”, umo¿liwiaj¹cy „doœwiadczalne” wyznaczenie liczby bitów
mantysy oraz tzw. epsylona maszynowego, dla zmiennych typu float i double, tj. najmniejszej liczby
e takiej, ¿e fl(e + 1) > 1. Jaki jest zwi¹zek e z precyzj¹ arytmetyki?
*/
int main(int argc, char** argv) {
	using namespace std;
	
	int m_d = 0, m_f = 0;	// bit count for double and float
	float e_f;				// float epsilon
	double e_d;				// double epsilon

	// FLOAT
	e_f = 1.f;
	while ((.5f*e_f + 1.f) > 1) {
		e_f *= .5f;
		m_f++;
	}
	cout << "Liczba bitow mantysy dla typu float: " << m_f << endl 
		<< "epsilon dla typu float: " << e_f << endl 
		<< "epsilon biblioteczny: " << FLT_EPSILON << endl 
		<< endl;

	// DOUBLE:
	e_d = 1.L;
	while ((.5L*e_d + 1.L) > 1) {
		e_d *= .5L;
		m_d++;
	}
	cout << "Liczba bitow mantysy dla typu double: " << m_d << endl 
		<< "epsilon dla typu double: " << e_d << endl 
		<< "epsilon biblioteczny: " << DBL_EPSILON << endl 
		<< endl;

	cout << "Nacisnij dowolny klawisz...";
	_getch();
	return EXIT_SUCCESS;
}// end main