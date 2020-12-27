#include "Matrix.h"
#include "Linear_system.h"
#include "Spline_Interpolator.h"


Spline_interpolator::Spline_interpolator(double* x, double* y, int size)
	: x(x)
	, y(y)
	, n_points(size)

{
}


double Spline_interpolator::interpolate(const double& x_new)
{
	//identify interval
	int interval=0;
	for (int i = 0; i < n_points-2; i = i ++)
	{
		if (x_new < x[i+1] && x_new > x[i])
		{
			interval = i;
			break;
		}
	}
	//interval = 1;
	//std::cout << "\nInterval " << interval << "\n";

	//calculate the polynomial and substitute with x at this interval
	double common_term = x[interval + 1] - x[interval];
	double term1 = (solutions.at(interval, 0)) * (x[interval + 1] - x_new)* (x[interval + 1] - x_new)* (x[interval + 1] - x_new) / (6 * common_term);
	double term2 = (solutions.at(interval + 1, 0)) * (x_new - x[interval])* (x_new - x[interval])* (x_new - x[interval]) / (6 * common_term);
	double term3 = ((y[interval]/(common_term)) - (solutions.at(interval, 0) *(common_term))/6.0) * (x[interval + 1] - x_new);
	double term4 =( (y[interval + 1] / common_term) - solutions.at(interval + 1, 0) * (common_term / 6)) * (x_new - x[interval]);

	return (term1 + term2 + term3 + term4); 

	
}

void Spline_interpolator::fit() 
{
	double* eqn = new double[n_points*n_points];
	double* b = new double[n_points];
	for (int k = 0; k < n_points-2; k++)
	{
		int i = 0;
		for (i = 0; i < k; i++)
			eqn[n_points*k+i] = 0;
		eqn[n_points*k+i] = x[i + 1] - x[i];
		eqn[n_points*k+i + 1] = 2 * (x[i + 2] - x[i]);
		eqn[n_points*k+i + 2] = x[i + 2] - x[i + 1];
		for (int j = i + 3; j < n_points; j++)
			eqn[n_points*k+j] = 0;
		b[k] = -(6 * (y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]) + 6 * (y[i] - y[i + 1]) / (x[i + 1] - x[i]));
	}
	for (int s = 0; s < n_points; s++)
	{
		eqn[n_points * (n_points - 2) + s] = 0;
		eqn[n_points * (n_points - 1) + s] = 0;
	}
	eqn[n_points * (n_points - 2) + 0] = 1; //first element (for first knot)
	eqn[n_points * (n_points - 1) + n_points-1] = 1; //last element (for last knot)
	b[n_points - 2] = 0;
	b[n_points - 1] = 0; 

	Matrix A(n_points, n_points, eqn);
	Matrix B(n_points, 1, b);
	delete [] eqn, b; 
	Linear_system sys(A, B);
	solutions = sys.solve();

}
