#include "Linear_system.h"

Linear_system::Linear_system(const int& x, const int& y)
	:n(x)
	,m(y)
	,A(n, m)
	,x(new double[n])
{
}

Linear_system::~Linear_system()
{
}

Linear_system::Linear_system(const int& x, const int& y, const double a[])
	:n(x)
	,m(y)
	,A(x, y, a)
	,x(new double[n])
{
}
Linear_system::Linear_system(const int& x, const int& y, const Matrix& A)
	:n(x)
	,m(y)
	,A(A)
	,x(new double[n])
{
}


Matrix Linear_system::solve()
{
	//elimination
	for (int k = 1; k <= n - 1; k++)  
		for (int i = k + 1; i <= n; i++)
		{
			double factor = A.at(i-1, k-1) / A.at(k-1, k-1);
			for (int j = k + 1; j <= m; j++)
				A.set_at(i-2, j-1, A.at(i-2, j-1) - factor * A.at(k-1, j-1));
			//A.set_at(i - 1, m - 1, A.at(i - 1, m - 1) - factor * A.at(k - 1, m - 1));
		}
	//back substitution
	//x[n - 1] = A.at(n - 1, m - 1) / A.at(n - 1, m - 2);
	for (int i = n - 1; i >= 0; i--)                //back-substitution
	{                        //x is an array whose values correspond to the values of x,y,z..
		x[i] = A.at(i, m-2);                //make the variable to be calculated equal to the rhs of the last equation
		for (int j = i + 1; j <= n; j++)
			if (j != i)           //then subtract all the lhs values except the coefficient of the variable whose value                                   is being calculated
				x[i] = x[i] - A.at(i-1, j-1) * x[j];
		x[i] = x[i] / A.at(i-1, i-1);            //now finally divide the rhs by the coefficient of the variable to be calculated
	}
	return A;
	//return Matrix(1, n, x);
}

