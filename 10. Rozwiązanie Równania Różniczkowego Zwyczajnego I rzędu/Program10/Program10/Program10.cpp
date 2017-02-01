#include<iostream>
#include <conio.h>
#include<fstream>
#define _USE_MATH_DEFINES
#include<math.h>
#include<string>

#include "Difference.h"
#include "NonLinearEquationSolving.h"

/* Napisz program w języku „C/C++”, rozwiązujący równanie różniczkowe zwyczajne pierwszego rzędu:
	dy(t)/dt + (10t^2 + 20)/(t^2 + 1) * [y(t) - 1] = 0, określone dla zmiennej t >= 0,
z warunkiem początkowym y(0) = 0, za pomocą metod:

	(a) bezpośredniej Eulera,
	(b) pośredniej Eulera,
	(c) metody trapezów.

Dla metod (b) i (c) wykonaj oddzielne rysunki przedstawiające po dwa wykresy: 
wykres przykładowego rozwiązania numerycznego oraz (dla porównania) wykres rozwiązania
analitycznego: y(t) = 1 - exp{-10[t + arctg(t)]}. Oba wykresy winny przedstawiać zależność y od
zmiennej niezależnej t. Rozwiązania analityczne zaznacz linią ciągłą, a numeryczne punktami. W
przypadku metody (a) wykonaj dwa takie rysunki: jeden uzyskany w warunkach numerycznej
stabilności metody, a drugi w warunkach numerycznej niestabilności. Wyjaśnij różnice pomiędzy
uzyskanymi wykresami.

Pokaż, że rząd dokładności uzyskanych stabilnych rozwiązań numerycznych jest zgodny z
przewidywaniami teoretycznymi. W tym celu wykonaj (na jednym rysunku) wykresy
przedstawiające zależności maksymalnych błędów bezwzględnych rozwiązań uzyskanych trzema
metodami, od kroku sieci czasowej (delta)t, posługując się skalą logarytmiczną (tzn. wykresy zależności
log10|błędu| od log10(delta)t ). Na podstawie wykresów wyznacz doświadczalnie rzędy dokładności
rozwiązań uzyskanych za pomocą różnych metod i porównaj je z rzędami teoretycznymi. O ile to
możliwe, zidentyfikuj też wartości kroku sieci poniżej których pojawia się wpływ błędów
maszynowych.
*/

using type = double;
using iterate = __int64;

type t_k_1 = {}, t_k = {}; // poziom czasowy: k + 1, k
type y_k = {};
type __dt = {};

type y(type t) { return (1. - exp(-10.*(t + atan(t)))); }
type f(type tk, type yk) { return (-((10.*tk*tk + 20.)/(tk*tk + 1)) * (yk - 1.)); }
type _fmpe(type yk_1) { return ((yk_1 - y_k)/__dt - f(t_k_1, yk_1)); }
type _fmt(type yk_1) { return ((yk_1 - y_k) / __dt - (f(t_k, y_k) + f(t_k_1, yk_1))/2.); }
type error(type yi, type y) { return abs(yi - y); }

bool S(type x, type s) { return ((x >= (s - 10.e-16)) && (x <= (s + 10.e-16))); }

int main(int argc, char** argv) {

	using namespace std;

	// Przedział [0, inf]:
	register type a_ = 0, _b = 0.5;

	// Warunek początkowy y(0) = 0:
	register type y0 = 0.;

	// Siatka:
	register type dt_step = 0.1;
	register type dt_ = 10.e-7, _dt = 10.e-1, _dt_ = pow(10., dt_step);

	register type tk, yk_1, yk;
	register type max_err, err;

	char file_name[32];

	// Pliki z błędami:
	ofstream error_mbe, error_mpe, error_mt;
	//error_mbe.open("errors_mbe.dat", ofstream::out);
	//error_mpe.open("errors_mpe.dat", ofstream::out);
	//error_mt.open("errors_mt.dat", ofstream::out);
	ofstream file;

	// Ustawienia dla rozwiązywania równania nieliniowego:
	NLES::setDelta(0.001);
	NLES::setEpsilon(0.001);
	NLES::setMax(100);
	NLES::saveEachResult(false);

	if (error_mbe && error_mpe && error_mt) {
		int pl = 0;
		//for (type dt = dt_; dt <= _dt; dt *= _dt_) {
			type dt = 0.2;
			iterate N;
			if (S(dt, 0.2))
				N = (iterate)((5. - a_) / dt);
			else
				N = (iterate)((_b - a_) / dt);
			cout << "Obliczenia dla dt: " << dt << " Liczba punktow: " << N << endl;

			// Warunek do zapisywwania wyników (tylko wyniki dla poszczególnych dt):
			//bool p = S(dt, 10.e-4) || S(dt, 10.e-3) || S(dt, 10.e-1);






			/* *********** Metoda bezpośrednia Euler'a: *********** */
			//if (p) {
				memset((void*)file_name, 0, 32);
				sprintf_s(file_name, 32, "results_mbe%d.dat", 3);
				file.open(file_name, ofstream::out);
			//}

			tk = a_; yk = y0;
			max_err = 0.;
			for (iterate i = 1; i <= (N+1)/5; ++i) {

				yk_1 = yk + dt*f(tk, yk);

				// Zapisywanie wyników w formacie: (t0 y0 t1 y1 ... ti yi ... tn yn):
				//if(p && file) 
					file << tk << " " << yk << " ";

				err = error(yk, y(tk));
				if (err > max_err)
					max_err = err;

				yk = yk_1;
				tk += dt;




				yk_1 = yk + dt*f(tk, yk);

				// Zapisywanie wyników w formacie: (t0 y0 t1 y1 ... ti yi ... tn yn):
				//if (p && file)
					file << tk << " " << yk << " ";

				err = error(yk, y(tk));
				if (err > max_err)
					max_err = err;

				yk = yk_1;
				tk += dt;




				yk_1 = yk + dt*f(tk, yk);

				// Zapisywanie wyników w formacie: (t0 y0 t1 y1 ... ti yi ... tn yn):
				//if (p && file)
					file << tk << " " << yk << " ";

				err = error(yk, y(tk));
				if (err > max_err)
					max_err = err;

				yk = yk_1;
				tk += dt;



				yk_1 = yk + dt*f(tk, yk);

				// Zapisywanie wyników w formacie: (t0 y0 t1 y1 ... ti yi ... tn yn):
				//if (p && file)
					file << tk << " " << yk << " ";

				err = error(yk, y(tk));
				if (err > max_err)
					max_err = err;

				yk = yk_1;
				tk += dt;




				yk_1 = yk + dt*f(tk, yk);

				// Zapisywanie wyników w formacie: (t0 y0 t1 y1 ... ti yi ... tn yn):
				//if (p && file)
					file << tk << " " << yk << " ";

				err = error(yk, y(tk));
				if (err > max_err)
					max_err = err;

				yk = yk_1;
				tk += dt;
			}

			//if (p && file) file.close();
			file.close();

			// Zapisywanie błędów:
			//error_mbe << log10(dt) << " " << log10(max_err) << " ";*/









			/* *********** Metoda pośrednia Euler'a: *********** 
			
			if (p) {
				memset((void*)file_name, 0, 32);
				sprintf_s(file_name, 32, "results_mpe%d.dat", pl);
				file.open(file_name, ofstream::out);
			}

			tk = a_; yk = y0;
			max_err = 0.;
			__dt = dt;
			for (iterate i = 1; i <= (N + 1)/5; ++i) {
				
				y_k = yk;
				t_k_1 = tk + dt;
				NLES::setFunction(_fmpe);
				Result r = NLES::secant(0., 1.5);

				yk_1 = r.x;

				// Zapisywanie wyników w formacie: (t0 y0 t1 y1 ... ti yi ... tn yn):
				if (p && file)
					file << tk << " " << yk << " ";

				err = error(yk, y(tk));
				if (err > max_err)
					max_err = err;

				yk = yk_1;
				tk += dt;




				y_k = yk;
				t_k_1 = tk + dt;
				NLES::setFunction(_fmpe);
				r = NLES::secant(0., 1.5);

				yk_1 = r.x;

				// Zapisywanie wyników w formacie: (t0 y0 t1 y1 ... ti yi ... tn yn):
				if (p && file)
					file << tk << " " << yk << " ";

				err = error(yk, y(tk));
				if (err > max_err)
					max_err = err;

				yk = yk_1;
				tk += dt;





				y_k = yk;
				t_k_1 = tk + dt;
				NLES::setFunction(_fmpe);
				r = NLES::secant(0., 1.5);

				yk_1 = r.x;

				// Zapisywanie wyników w formacie: (t0 y0 t1 y1 ... ti yi ... tn yn):
				if (p && file)
					file << tk << " " << yk << " ";

				err = error(yk, y(tk));
				if (err > max_err)
					max_err = err;

				yk = yk_1;
				tk += dt;



				y_k = yk;
				t_k_1 = tk + dt;
				NLES::setFunction(_fmpe);
				r = NLES::secant(0., 1.5);

				yk_1 = r.x;

				// Zapisywanie wyników w formacie: (t0 y0 t1 y1 ... ti yi ... tn yn):
				if (p && file)
					file << tk << " " << yk << " ";

				err = error(yk, y(tk));
				if (err > max_err)
					max_err = err;

				yk = yk_1;
				tk += dt;



				y_k = yk;
				t_k_1 = tk + dt;
				NLES::setFunction(_fmpe);
				r = NLES::secant(0., 1.5);

				yk_1 = r.x;

				// Zapisywanie wyników w formacie: (t0 y0 t1 y1 ... ti yi ... tn yn):
				if (p && file)
					file << tk << " " << yk << " ";

				err = error(yk, y(tk));
				if (err > max_err)
					max_err = err;

				yk = yk_1;
				tk += dt;
			}

			if (p && file) file.close();

			// Zapisywanie błędów:
			error_mpe << log10(dt) << " " << log10(max_err) << " ";






			/* *********** Metoda pośrednia trapezów: *********** */

			/*if (p) {
				memset((void*)file_name, 0, 32);
				sprintf_s(file_name, 32, "results_mt%d.dat", pl);
				file.open(file_name, ofstream::out);
			}

			tk = a_; yk = y0;
			max_err = 0.;
			__dt = dt;
			for (iterate i = 1; i <= N + 1; ++i) {

				y_k = yk;
				t_k = tk;
				t_k_1 = tk + dt;
				NLES::setFunction(_fmt);
				Result r = NLES::secant(0., 1.5);

				yk_1 = r.x;

				// Zapisywanie wyników w formacie: (t0 y0 t1 y1 ... ti yi ... tn yn):
				if (p && file)
					file << tk << " " << yk << " ";

				err = error(yk, y(tk));
				if (err > max_err)
					max_err = err;

				yk = yk_1;
				tk += dt;
			}

			if (p && file) file.close();

			// Zapisywanie błędów:
			error_mt << log10(dt) << " " << log10(max_err) << " ";

			*/

			//if (p) ++pl;
		//}
		//error_mbe.close();
		//error_mpe.close();
		//error_mt.close();
	}

	cout << "Press any button to finish...";
	_getch();
	return EXIT_SUCCESS;
}// end main