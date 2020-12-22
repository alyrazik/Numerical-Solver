#pragma once
#include <iostream>
class Matrix
{
private:
	int n;
	int m;
	double* array;
public:
	Matrix(const int& n, const int& m);
	Matrix(const int& n, const int& m, const double a[]);
	~Matrix();
	double at(const int& i, const int& j) const;
	void set_at(const int& i, const int& j, const double& number);
	friend std::ostream& operator<<(std::ostream& out, const Matrix& m);

	

};

