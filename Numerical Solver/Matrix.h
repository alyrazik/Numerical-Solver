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
	Matrix(const Matrix& m2);
	Matrix(const int& n, const int& m, const double a[]);
	~Matrix();
	int n_rows() const;
	int n_cols() const;
	Matrix augment(const Matrix& y) const;
	double at(const int& i, const int& j) const;
	void set_at(const int& i, const int& j, const double& number);
	friend std::ostream& operator<<(std::ostream& out, const Matrix& m);

	

};

