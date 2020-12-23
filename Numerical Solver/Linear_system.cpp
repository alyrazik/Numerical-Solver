#include "Linear_system.h"
double ILL_CONDITIONING = 0.001;

Linear_system::Linear_system(const int& x, const int& y)
	:valid_solution(false)
	,n(x)
	,m(y)
	,A(n, m)
	,x(new double[n])
{
}

Linear_system::~Linear_system()
{
}

Linear_system::Linear_system(const int& x, const int& y, const double a[])
	:valid_solution(false)
	,n(x)
	,m(y)
	,A(x, y, a)
	,x(new double[n])
{
}
Linear_system::Linear_system(const int& x, const int& y, const Matrix& A)
	:valid_solution(false)
	,n(x)
	,m(y)
	,A(A)
	,x(new double[n])
{
}

Linear_system::Linear_system(const Matrix& A, const Matrix& b)
	:valid_solution(false)
	, n(A.n_rows())
	, m(1+A.n_rows())
	, A(A.augment(b))
	, x(new double[n])
{
}


Matrix Linear_system::solve()
{
	//elimination
	for (int k = 0; k < n; k++)  
		for (int i = k + 1; i < n; i++)
		{
			double factor = A.at(i, k) / A.at(k, k);
			for (int j = k ; j <m+1; j++) //m+1 since our A matrix is augmented, we start at k+1 since at k, answer is known as 0.
				A.set_at(i, j, A.at(i, j) - factor * A.at(k, j));
		}
	//back substitution
	for (int i = n - 1; i >= 0; i--)                
	{                       
		x[i] = -(A.at(i, m-1));                
		for (int j = i ; j <n; j++)
			if (j != i)          
				x[i] = x[i] - A.at(i, j) * x[j];
		x[i] = x[i] / A.at(i, i);
	}
	double diagonal=1;
	for (int i = 0; i < n; i++)
		diagonal = diagonal * A.at(i, i);
	valid_solution = (diagonal >= ILL_CONDITIONING);
	return Matrix(n,1,x);
}

