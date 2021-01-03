#include <iostream>
using namespace std;
#include "Linear_system.h"
#include "Matrix.h"
#include "univar_regressor.h"
int main()
{
	//double temp[20] = { 1, 1, -3, 1, -2, -3, 3, -4, 1, 0, 1, 0, 2, -1, -1, 1, 2, 0, 0, -12 };
	//double t[3] = { 1, 2, 3 };
	//double as[16] = { 1, 1, -3, 1, -3, 3, -4, 1, 1, 0, 2, -1, 1, 2, 0, 0 };
	//double bs[4] = { -2, 0, -1, -12 };

	int i, m, N;
	cout.precision(4);                        //set precision
	cout.setf(ios::fixed);
	cout << "\nEnter the no. of data pairs to be entered:\n";
	cin >> N;
	double *x = new(double[N]);
	double *y = new(double[N]);
	cout << "\nEnter the x-axis values:\n";                
	for (i = 0; i<N; i++)
		cin >> x[i];
	cout << "\nEnter the y-axis values:\n";                
	for (i = 0; i<N; i++)
		cin >> y[i];
	cout << "\nEnter the degree of Polynomial you want to fit:\n";
	cin >> m;                                
	if (N < m + 1)
	{
		cout << "\nYou need more samples to fit a polynomial of this degree";
		return -1;
	}
	

	univar_regressor poly_regressor = univar_regressor(*x, *y, m);
	Matrix sol = poly_regressor.fit(x,y,m);

	//Matrix out = A.augment(b); //don't forget to multiply with -1 to the b elements.
	//cout << out << endl;
	//cout << A << endl;
	//cout << b << endl;
	//cout << x << endl;
	//Linear_system x(A, b);
	//Matrix sol = x.solve();
	cout << "solution" << endl << sol;
	//cout << "system is " << x.is_valid_solution();
	//Matrix tt(3, 1, t);
	//Matrix x(3, 4, temp);
	//cout << x;
	//cout << " and the b matrix is" << endl;
	//cout << tt; 
	//Linear_system A(3, 3, x, tt);
	//Linear_system B(4, 5, temp);
	//Matrix x = B.solve();
	//cout << x;
	//cout << tt;

	int ss;
	cin >> ss;
	return 0;
}