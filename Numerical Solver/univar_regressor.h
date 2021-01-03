#pragma once
#include <array>
#include "Matrix.h"

class univar_regressor
{
private: 
	int n, m;
	double *x, *y;


public:
	univar_regressor(const double x[], const double y[], const int m);
	~univar_regressor();

	Matrix fit(const double x[], const double y[], const int n, const int m);
	double predict(const double x, const Matrix coeff, const int m);
};

