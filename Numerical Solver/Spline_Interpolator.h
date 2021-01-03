#pragma once
class Spline_interpolator
{
	double* x;
	double* y;
	int n_points;
	Matrix solutions;

public:
	Spline_interpolator(double* x, double* y, int size);


	double interpolate(const double& x_new);

	void fit() ;
	
};

