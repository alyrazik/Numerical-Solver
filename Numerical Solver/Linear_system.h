#pragma once
#include <array>
#include "Matrix.h"
class Linear_system
{
private:
	bool valid_solution;
	const int n;
	const int m;
	Matrix A;
	double* x;

public:
	//default constructor
	Linear_system(const int& x, const int& y);
	Linear_system(const int& x, const int& y, const double a[]);
	Linear_system(const int& x, const int& y, const Matrix& A);
	~Linear_system();
	
	Matrix solve();

};

