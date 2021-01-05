#pragma once
#include "Matrix.h"
class Newton_interpolator
{
	double* x;
	double* y; 
	int n_points;
	double* coefficients;
	double finite_difference(const int& first, const int& last) const;


public:
	Newton_interpolator(double* x, double* y, int size);
	~Newton_interpolator() { delete x, y, coefficients; }
	double interpolate(const double& x) ;
	double* fit() const;
};

