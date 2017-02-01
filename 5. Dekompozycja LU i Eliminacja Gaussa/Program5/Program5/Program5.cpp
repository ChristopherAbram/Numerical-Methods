#include<iostream>
#include<conio.h>

#include "SolvingLinearSystem.h"

/*
Dana jest macierz A:

	1   -20   30  -4
	2   -40  -6    50
	9   -180  11  -12
   -16   15  -140  13

oraz wektor b:

	 35
	 104
	-366
	-354

Napisz program w jêzyku „C/C++”, realizuj¹cy dekompozycjê LU macierzy A, przy zastosowaniu
eliminacji Gaussa z czêœciowym wyborem elementu podstawowego, a nastêpnie rozwi¹zuj¹cy uk³ad
równañ Ax = b. Uwaga: nale¿y zrealizowaæ wariant dekompozycji omawiany na wyk³adzie. */

int main(int argc, char** argv) {

	using namespace std;
	using type = double;

	Matrix<type> A(4, 4),
				 L(4, 4),   // Lower Matrix
				 U(4, 4),   // Upper Matrix
				 P(4, 4);   // Permutated I matrix according to p. It is useful to PA = LU

	vector<type> b(4);

	// Set of n indexes which could be permutated:
	vector<long> p(A.getRows());
	vector<long> j = { 0,1,2,3 }; // natural order

	/* Initialization of matrixes */
	A = {
		{1., -20., 30., -4.},
		{2., -40., -6., 50.},
		{9., -180., 11., -12.},
		{-16., 15., -140., 13.}
	};

	b = {
		35.,
		104.,
		-366.,
		-354.
	};

	/* LU decomposition */
	U = A; // Our matrix
	try {
		LU(L, U, p, PARTIAL_SELECTION);
	}
	catch (MatrixException e) {
		cout << e.what() << endl;
	}

	/* Getting P matrix */
	P = Kronecker<type>(p, j);

	/* Presentation of results */
	cout << endl << "Lower matrix:" << endl;
	print(P*L);
	cout << endl << "Upper matrix:" << endl;
	print(P*U);
	cout << endl << "PA = LU:" << endl;
	print((P*L)*(P*U));
	cout << endl;

	/* Upper and Lower matrix in one matrix */
	U = A;
	LU2(U, p, PARTIAL_SELECTION);

	cout << "Lower and Upper matrix in one matrix:" << endl;
	print(Kronecker<type>(p, j)*U);
	cout << endl;

	/* Showing step by step how to solve system */
	U = A;
	U.appendColumns(Matrix<type>(b).transpose());
	cout << "Given linear system:" << endl;
	print(U);
	cout << endl;

	/* Solving system Ax = b */
	SLS<type> solver(A, b);
	solver.solveGaussianElimination(p);
	b = solver.getSolutions();

	cout << endl << "Solutions of the system Ax = b:" << endl;
	for (long i = 0; i < b.size(); ++i) {
		cout.precision(6);
		cout << " x" << i + 1 << " = " << b[p[i]] << endl;
	}

	cout << endl << endl << "Press any button to close...";
	_getch();
	return EXIT_SUCCESS;
}// end main