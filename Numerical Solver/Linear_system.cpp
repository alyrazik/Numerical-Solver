#include "Linear_system.h"
#include "stdlib.h"
double ILL_CONDITIONING = 0.0001;

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
	//pivoting
	//initialize the S and L vectors
	double* maxs = new double[n]; //scale vector
	for (int i = 0; i < n; i++)
		maxs[i] = 0;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			if (abs(A.at(i, j)) > maxs[i])
				maxs[i] = abs(A.at(i, j));
	int* L = new int[n]; //index vector
	for (int i = 0; i < n; i++)
		L[i] = i;
	
	double max;
	int max_index;

	//elimination
	for (int k = 0; k < n; k++) // for each pivot equation
	{
		//update the L vector
		max = 0;
		max_index = 0;
		for (int i = k; i < n; i++)
		{
			double scaled = abs(A.at(L[i], k)) / maxs[i];
			if (scaled > max)
			{
				max = scaled;
				max_index = i;
			}
		}
		int temp = L[k];
		L[k] = L[max_index];
		L[max_index] = temp;
		//checked that L and maxs are calculated correctly for first iteration, k=0;
		for (int i = k + 1; i < n; i++) //for all subsequent equations
		{
			double factor = A.at(L[i], k) / A.at(L[k], k);
			for (int j = k; j < m ; j++) //for all elements in the equation
				//m+1 since our A matrix is augmented, we start at k+1 since at k, answer is known as 0.
				A.set_at(L[i], j, A.at(L[i], j) - factor * A.at(L[k], j));
		}
	}

	//back substitution
	for (int i = n - 1; i >= 0; i--)                
	{                       
		x[i] = -(A.at(L[i], m-1));                
		for (int j = i ; j <n; j++)
			if (j != i)          
				x[i] = x[i] - A.at(L[i], j) * x[j];
		x[i] = x[i] / A.at(L[i], i);
	}

	//check for ill-conditioned systems
	double diagonal = 1;
	for (int i = 0; i < n; i++)
		diagonal = diagonal * A.at(L[i], i);
	valid_solution = (diagonal > ILL_CONDITIONING);
	
	return Matrix(n,1,x);
	//return A;
}

