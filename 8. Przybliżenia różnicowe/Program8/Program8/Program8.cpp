#include<iostream>
#include<conio.h>
#include<cmath>
#include<windows.h>
#include<fstream>

#define PI 3.14159265358979323846

#include "Difference.h"

/**
Napisz program w jêzyku „C/C++”, obliczaj¹cy przybli¿one wartoœci pierwszych pochodnych
funkcji f(x) = cos(x) w punktach koñcowych i œrodkowym przedzia³u [0, pi/2] zmiennej x. Zastosuj
wszystkie omawiane na wyk³adzie i na æwiczeniach przybli¿enia ró¿nicowe dwupunktowe i
trzypunktowe (jednostronne b¹dŸ centralne, w zale¿noœci od po³o¿enia punktu w przedziale) na sieci
jednorodnej o kroku h. 

Wykonaj (na jednym rysunku) wykresy przedstawiaj¹ce zale¿noœci b³êdów
bezwzglêdnych przybli¿eñ ró¿nicowych od kroku sieci, pos³uguj¹c siê skal¹ logarytmiczn¹ (tzn.
wykresy zale¿noœci log10|b³êdu| od log10 h ). Na podstawie wykresów wyznacz doœwiadczalnie rzêdy
dok³adnoœci przybli¿eñ ró¿nicowych. SprawdŸ, czy tak wyznaczone rzêdy dok³adnoœci pokrywaj¹
siê z rzêdami teoretycznymi i wyjaœnij ewentualne rozbie¿noœci. Ponadto zidentyfikuj wartoœci
kroku sieci poni¿ej których pojawia siê wp³yw b³êdów maszynowych. Obliczenia powtórz dla
dwóch typów zmiennych rzeczywistych (float, i double) i porównaj wyniki.

Uwaga: najwygodniej jest zastosowaæ wzorzec funkcji (function template) z typem zmiennych jako
parametrem wzorca.
*/

HANDLE  hConsole;

float f_float(float x) { return cos(x); }
float df_float(float x) { return -sin(x); }

double f_double(double x) { return cos(x); }
double df_double(double x) { return -sin(x); }

int main(int argc, char** argv) {

	using namespace std;
	hConsole = GetStdHandle(STD_OUTPUT_HANDLE);

	try {
		/* FLoating point calculations */
		using type = float;

		type error;
		Difference<type> diff;
		diff.setFunction(f_float);

		ofstream out_file;
		out_file.open("results_float.dat", ofstream::out);

		if (out_file) {
			type begin = 10e-17f, end = 10e-2f;
			type x;

			/* For x = 0 */
			x = 0.f;
			// Forward difference two-points:
			for (type h = begin; h <= end; h *= 10.f) {
				diff.setIncrease(h);
				error = log10(abs(diff.forward(x) - df_float(x)));
				if(isfinite(error))
					out_file << log10(abs(h)) << " " << error << " ";
			}
			out_file << endl;
			// Forward difference three-points:
			for (type h = begin; h <= end; h *= 10.f) {
				diff.setIncrease(h);
				error = log10(abs(diff.forward(x, Difference<type>::Points::THREE) - df_float(x)));
				if (isfinite(error))
					out_file << log10(abs(h)) << " " << error << " ";
			}
			out_file << endl;



			/* For x = pi/4 */
			x = PI / 4.f;
			// Forward difference two-points:
			for (type h = begin; h <= end; h *= 10.f) {
				diff.setIncrease(h);
				error = log10(abs(diff.forward(x) - df_float(x)));
				if (isfinite(error))
					out_file << log10(abs(h)) << " " << error << " ";
			}
			out_file << endl;
			// Central difference two-points:
			for (type h = begin; h <= end; h *= 10.f) {
				diff.setIncrease(h);
				error = log10(abs(diff.central(x) - df_float(x)));
				if (isfinite(error))
					out_file << log10(abs(h)) << " " << error << " ";
			}
			out_file << endl;
			// Backward difference two-points:
			for (type h = begin; h <= end; h *= 10.f) {
				diff.setIncrease(h);
				error = log10(abs(diff.backward(x) - df_float(x)));
				if (isfinite(error))
					out_file << log10(abs(h)) << " " << error << " ";
			}
			out_file << endl;




			/* For x = pi/2 */
			x = PI / 2.f;
			// Backward difference two-points:
			for (type h = begin; h <= end; h *= 10.f) {
				diff.setIncrease(h);
				error = log10(abs(diff.backward(x) - df_float(x)));
				if (isfinite(error))
					out_file << log10(abs(h)) << " " << error << " ";
			}
			out_file << endl;
			// Backward difference three-points:
			for (type h = begin; h <= end; h *= 10.f) {
				diff.setIncrease(h);
				error = log10(abs(diff.backward(x, Difference<type>::Points::THREE) - df_float(x)));
				if (isfinite(error))
					out_file << log10(abs(h)) << " " << error << " ";
			}
			out_file << endl;

			out_file.close();
		}
	}
	catch (DifferenceException& de) {
		cout << "FLOAT: " << de.what() << endl;
	}



	try {
		/* Double point calculations */
		using type = double;

		type error;
		Difference<type> diff;
		diff.setFunction(f_double);

		ofstream out_file;
		out_file.open("results_double.dat", ofstream::out);

		if (out_file) {
			type begin = 10.e-17, end = 10.e-2;
			type x;

			/* For x = 0 */
			x = 0.;
			// Forward difference two-points:
			for (type h = begin; h <= end; h *= 10.) {
				diff.setIncrease(h);
				error = log10(abs(diff.forward(x) - df_float(x)));
				if (isfinite(error))
					out_file << log10(abs(h)) << " " << error << " ";
			}
			out_file << endl;
			// Forward difference three-points:
			for (type h = begin; h <= end; h *= 10.) {
				diff.setIncrease(h);
				error = log10(abs(diff.forward(x, Difference<type>::Points::THREE) - df_float(x)));
				if (isfinite(error))
					out_file << log10(abs(h)) << " " << error << " ";
			}
			out_file << endl;



			/* For x = pi/4 */
			x = PI / 4.;
			// Forward difference two-points:
			for (type h = begin; h <= end; h *= 10.) {
				diff.setIncrease(h);
				error = log10(abs(diff.forward(x) - df_float(x)));
				if (isfinite(error))
					out_file << log10(abs(h)) << " " << error << " ";
			}
			out_file << endl;
			// Central difference two-points:
			for (type h = begin; h <= end; h *= 10.) {
				diff.setIncrease(h);
				error = log10(abs(diff.central(x) - df_float(x)));
				if (isfinite(error))
					out_file << log10(abs(h)) << " " << error << " ";
			}
			out_file << endl;
			// Backward difference two-points:
			for (type h = begin; h <= end; h *= 10.) {
				diff.setIncrease(h);
				error = log10(abs(diff.backward(x) - df_float(x)));
				if (isfinite(error))
					out_file << log10(abs(h)) << " " << error << " ";
			}
			out_file << endl;



			/* For x = pi/2 */
			x = PI / 2.;
			// Backward difference two-points:
			for (type h = begin; h <= end; h *= 10.) {
				diff.setIncrease(h);
				error = log10(abs(diff.backward(x) - df_float(x)));
				if (isfinite(error))
					out_file << log10(abs(h)) << " " << error << " ";
			}
			out_file << endl;
			// Backward difference three-points:
			for (type h = begin; h <= end; h *= 10.) {
				diff.setIncrease(h);
				error = log10(abs(diff.backward(x, Difference<type>::Points::THREE) - df_float(x)));
				if (isfinite(error))
					out_file << log10(abs(h)) << " " << error << " ";
			}
			out_file << endl;

			out_file.close();
		}
	}
	catch (DifferenceException& de) {
		cout << "Double: " << de.what() << endl;
	}

	cout << endl << endl << "Press any button to finish...";
	_getch();
	return EXIT_SUCCESS;
}// end main