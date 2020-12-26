#include "Newton_interpolator.h"

double Newton_interpolator::finite_difference(const int& first, const int& last) const
{
	int diff = last - first;
	if (diff == 1)
	{
		return (y[last] - y[first]) / (x[last] - x[first]);
	}
	return ((finite_difference(first+1, last)) - (finite_difference(first, last-1))) / (x[last] - x[first]);
}

Newton_interpolator::Newton_interpolator(double* x, double* y, int size)
	:x(x)
	,y(y)
	,n_points(size)
	,coefficients(new double[n_points])
{
}

double Newton_interpolator::interpolate(const double& x_new)
{
	double answer=0.0;
	for (int n = 1; n < n_points; n++)
	{
		double miqdar = coefficients[n];
		for (int i = 0; i < n; i++)
		{
			miqdar *= (x_new - x[i]);
		}
		answer += miqdar;
	}
	return answer;
}

void Newton_interpolator::fit() const
{
	for (int i = 1; i < n_points; i++)
	{
		coefficients[i] = finite_difference(0, i);
	}
	coefficients[0] = y[0];
	for (int i = 0; i < n_points; i++)
		std::cout << coefficients[i] << "\n";
}
