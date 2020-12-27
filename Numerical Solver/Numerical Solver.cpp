// Numerical Solver.cpp : This file contains the 'main' function. Program execution begins and ends there.
#include <iostream>
using namespace std;
#include "Linear_system.h"
#include "Matrix.h"
#include "Newton_interpolator.h"
#include "Spline_Interpolator.h"

int main()
{
	double* x = new double[4];
	x[0] = 3; x[1] = 4.5; x[2] = 7; x[3] = 9;
	double* y = new double[4];
	y[0] = 2.5; y[1] = 1.0; y[2] = 2.5; y[3] = 0.5;
	Newton_interpolator n(x, y, 4);
	n.fit();
	double answer = n.interpolate(5);
	cout << "answer: " << answer;
	

	// spline interpolation:
	Spline_interpolator s(x, y, 4);
	s.fit();
	double answer2 = s.interpolate(5);
	cout << "answer using splines is: " << answer2; 

	int temp;
	cin >> temp;
}
