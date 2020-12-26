// Numerical Solver.cpp : This file contains the 'main' function. Program execution begins and ends there.
#include <iostream>
using namespace std;
#include "Linear_system.h"
#include "Matrix.h"
#include "Newton_interpolator.h"

int main()
{
	double* x = new double[4];
	x[0] = 0; x[1] = 1; x[2] = 2; x[3] = 3;
	double* y = new double[4];
	y[0] = 0; y[1] = 1; y[2] = 8; y[3] = 27;
	Newton_interpolator n(x, y, 4);
	n.fit();
	double answer = n.interpolate(5);
	cout << "answer: " << answer;
	int i; 
	cin >> i;
}
