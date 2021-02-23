#pragma once
#include <array>
#include "Matrix.h"

class univar_regressor
{
private: 
	int n, m;
	double *x, *y;
	bool valid_solution;


public:
	int seidel_iterations;
	univar_regressor(const double x[], const double y[], const int m, const int s);
	~univar_regressor();

	Matrix fit(const double x[], const double y[], const int n, const int m, const int s);
	double predict(const double x, const Matrix coeff, const int m);
	bool is_valid_solution() { return valid_solution; }
};

