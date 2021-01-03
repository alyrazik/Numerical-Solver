#include "univar_regressor.h"
#include <iostream>
#include <array>
#include <vector>
#include <cmath>

using namespace std;
#include "Linear_system.h"
#include "Matrix.h"


univar_regressor::univar_regressor(const double x[], const double y[], const int m)
	:n(sizeof(x)+1)
	, m(m)
	, x(new double [n])
	, y(new double [n])
{
}


univar_regressor::~univar_regressor()
{
}


Matrix univar_regressor::fit(const double x[], const double y[], const int n, const int m) {

	int i, j;
	double *as = new(double[2 * m + 1]);              //Array that will store the values of sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
	for (i = 0; i < 2 * m + 1; i++)
	{
		as[i] = 0;
		for (j = 0; j < n; j++)
			as[i] = as[i] + pow(x[j], i);        //consecutive positions of the array will store N,sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
	}
	double *C = new double[(m+1)*(m+1)];        //the final coefficients array
	int cntr = 0;
	for (i = 0; i <= m; i++)
		for (j = 0; j <= m; j++)
		{
			C[cntr] = as[i + j];
			cntr++;
		}
			

	double* bs = new double[m + 1];                    //Array to store the values of sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
	for (i = 0; i < m + 1; i++)
	{
		bs[i] = 0;
		for (j = 0; j < n; j++)
			bs[i] = bs[i] + pow(x[j], i)*y[j];        //consecutive positions will store sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
	}

	Matrix A(m + 1, m + 1, C);
	Matrix b(m + 1, 1, bs);
	Matrix out = A.augment(b); //don't forget to multiply with -1 to the b elements.
	cout << "\nAugmented Matrix:\n" << endl;
	cout << out << endl;
	
	Linear_system S(A, b);
	Matrix sol = S.solve();

	return sol;
}

double univar_regressor::predict(const double xi, const Matrix coeff, const int m)
{
	double yi=0;
	for (int idx = 0; idx <= m; idx++)
		yi += coeff.at(idx, 0) * pow(xi, idx);
	return yi;
}
