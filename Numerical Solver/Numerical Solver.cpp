// Numerical Solver.cpp : This file contains the 'main' function. Program execution begins and ends there.
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
using namespace std;
#include "Linear_system.h"
#include "Matrix.h"
#include "univar_regressor.h"
#include "multivar_regressor.h"
#include "Newton_interpolator.h"
#include "Spline_Interpolator.h"

#include "gnuplot_i.hpp"



//Normalization functions >>> Used for regression

Matrix norm_x(const Matrix x, const int n)
{
	int i, j;
	Matrix x_tmp = Matrix(3, n);
	for (i = 0; i < n; i++)
		x_tmp.set_at(0, i, 1);
	//for (i = 1; i < 3; i++)
	//{
	//	double x_min = x.at(i, 0);
	//	double x_max = x.at(i, 0);
	//	for (int j = 1; j < n; j++) {
	//		if (x.at(i, j) < x_min) {
	//			x_min = x.at(i, j);
	//		}
	//		if (x.at(i, j) > x_max) {
	//			x_max = x.at(i, j);
	//		}
	//	}
	//	for (int j = 0; j < n; j++)
	//		x_norm.set_at(i, j, (x.at(i, j) - x_min) / (x_max - x_min));
	//}

	return x_tmp;
}
double* norm_y(const double y[], const int n) {
	int i;
	double* y_norm = new double[n];
	double y_min = y[0];
	double y_max = y[0];
	for (i = 0; i < n; i++) {
		if (y[i] < y_min)
			y_min = y[i];
		if (y[i] > y_max)
			y_max = y[i];
	}
	for (i = 0; i < n; i++)
		y_norm[i] = (y[i] - y_min) / (y_max - y_min);
	return y_norm;
}

int regress_line() {
	int i, n, m;
	ifstream  trainData;
	try {
		trainData.open("2_a_dataset_1.csv");
		cout << "Opened file: 2_a_dataset_1.csv";
	}
	catch (const char* cstr) {
		cerr << cstr << '\n';
		exit(1);
	}
	string line;
	vector<vector<string> > parsedCsv;
	while (getline(trainData, line))
	{
		stringstream lineStream(line);
		string cell;
		vector<string> parsedRow;
		while (getline(lineStream, cell, ','))
		{
			parsedRow.push_back(cell);
		}

		parsedCsv.push_back(parsedRow);
	}
	n = parsedCsv.size() - 1;
	double *x = new double[n];
	double *y = new double[n];

	for (size_t i = 1; i <= n; ++i) // start at i = 1 to skip header line
	{
		x[i - 1] = stod(parsedCsv[i][0]);
		y[i - 1] = stod(parsedCsv[i][1]);
	}

	cout << "\nWhat degree of Polynomial do you want to use for the fit?\n";
	cin >> m;

	univar_regressor polyRegress = univar_regressor(x, y, m);
	Matrix coeff = polyRegress.fit(x, y, n, m);


	//cout <<"system is "<< S.is_valid_solution();

	cout << "\nThe values of the solution coefficients are:\n";
	for (int idx = 0; idx <= m; idx++)
		cout << "x^" << idx << "=" << coeff.at(idx, 0) << endl;            // Print the values of x^0,x^1,x^2,x^3,....    
	cout << "\nThe fitted Polynomial is given by:\ny=";
	for (int idx = 0; idx <= m; idx++)
		cout << " + (" << coeff.at(idx, 0) << ")" << "x^" << idx;
	cout << "\n";


	double *y_pred = new double[n];
	for (i = 0; i < n; i++)
		y_pred[i] = polyRegress.predict(x[i], coeff, m);

	cout << "\n Predictions: \n" << endl;
	cout << "\nx   y_pred    y" << endl;
	cout << "=================" << endl;

	for (i = 0; i < n; i++) {
		cout << x[i] << " " << y_pred[i] << " " << y[i] << endl;
	}



	int ss;
	cin >> ss;
	return 0;
}

int regress_2d() {
	int i, j, n, m;
	m = 2;    // two dimensional input x: x_1, x_2
	ifstream  trainData;
	try {
		trainData.open("Housing.csv");
		cout << "Opened file: Housing.csv";
	}
	catch (const char* cstr) {
		cerr << cstr << '\n';
		exit(1);
	}
	string line;
	vector<vector<string> > parsedCsv;
	while (getline(trainData, line))
	{
		stringstream lineStream(line);
		string cell;
		vector<string> parsedRow;
		while (getline(lineStream, cell, ','))
		{
			parsedRow.push_back(cell);
		}

		parsedCsv.push_back(parsedRow);
	}
	n = parsedCsv.size() - 1;


	Matrix x = Matrix(m + 1, n);
	double *y = new double[n];

	for (i = 1; i <= n; ++i) // start at i = 1 to skip header line
	{
		x.set_at(0, i - 1, 1);        // put 1's in x0 for the algorithm to work, data points will be stored in x1, x2
		for (j = 1; j < m + 1; j++)
			x.set_at(j, i - 1, stod(parsedCsv[i][j - 1]));
		y[i - 1] = stod(parsedCsv[i][4]);
	}

	//Matrix x_norm = Matrix(m + 1, n);
	double *y_norm = new double[n];

	//x_norm = norm_x(x, n);
	y_norm = norm_y(y, n);

	//for (i = 0; i < 20; i++)
	//	cout << "\n" << y_norm[i];


	multivar_regressor polyRegress2D = multivar_regressor(x, y_norm, n, m);
	Matrix coeff = polyRegress2D.fit(x, y_norm, n, m);


	////cout <<"system is "<< S.is_valid_solution();

	cout << "\nThe values of the solution coefficients are:\n";
	for (int idx = 0; idx <= m; idx++)
		cout << "x_" << idx << "=" << coeff.at(idx, 0) << endl;            // Print the values of x^0,x^1,x^2,x^3,....    
	cout << "\nThe fitted Polynomial is given by:\ny=";
	for (int idx = 0; idx <= m; idx++)
		cout << " + (" << coeff.at(idx, 0) << ")" << "x_" << idx;
	cout << "\n";

	// Prediction
	double* xi = new double[3];
	double *y_pred = new double[n];
	for (i = 0; i < n; i++)
		for (j = 0; j < 3; j++)
			xi[j] = x.at(j, i);

	y_pred[i] = polyRegress2D.predict(xi, coeff, m);

	cout << "\n Predictions: \n" << endl;
	cout << "\nx1     x2      y_pred    y" << endl;
	cout << "============================" << endl;

	for (i = 0; i < 11; i++) {
		cout << x.at(1, i) << "       " << x.at(2, i) << "    " << y_pred[i] << " " << y[i] << endl;
	}



	int ss;
	cin >> ss;
	return 0;
}


void interpolate() {
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
}
int main()
{
	
	int number = 0;
	cout << "Chose what do you want to do: \n \t1: Single variable Polynomial Regression. \n\t2: 2D Ploynomial Regression. \n\t3: Interpolation. \n";
	cin >> number;
	
	if (number == 1)
	{
		cout << "Linear Regression" << endl;
		int out = regress_line();
		
	}
	else if (number == 2) {
		cout << "2D Ploynomial Regression" << endl;
		int out = regress_2d();
	}
	else if (number == 3) {
		cout << "Interpolation" << endl;
		interpolate();
	}
	else {
		cout << "Wrong Choice !" << endl;
	}
	int ss;
	cin >> ss;
	return 0;


}


