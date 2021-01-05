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

#include <cmath>

#include "gnuplot_i.hpp" //Gnuplot class handles POSIX-Pipe-communikation with Gnuplot


//Normalization functions >>> Used for normalizing housing data 

Matrix norm_x(const Matrix x, const int n)
{
	int i, j;
	Matrix x_tmp = Matrix(3, n);
	for (i = 0; i < n; i++)
		x_tmp.set_at(0, i, 1);
	for (i = 1; i < 3; i++)
	{
		double x_min = x.at(i, 0);
		double x_max = x.at(i, 0);
		for (int j = 1; j < n; j++) {
			if (x.at(i, j) < x_min) {
				x_min = x.at(i, j);
			}
			if (x.at(i, j) > x_max) {
				x_max = x.at(i, j);
			}
		}
		for (int j = 0; j < n; j++)
			x_tmp.set_at(i, j, (x.at(i, j) - x_min) / (x_max - x_min));
	}

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

double denorm_y(const double y_norm, const double y_min, const double y_max) {
	double y_denorm = y_norm * (y_max - y_min) + y_min;
	return y_denorm;
}

double calc_mse(const double y[], const double y_pred[], const int n) {
	double mse = 0;
	for (int i = 0; i < n; i++)
		mse = mse + pow((y[i] - y_pred[i]),2) /n;
	return mse;
}
void doSomePlotting(vector<double> x, vector<double> y)
{
	// if path-variable for gnuplot is not set, do it with:
// Gnuplot::set_GNUPlotPath("C:/program files/gnuplot/bin/");

// set a special standard terminal for showonscreen (normally not needed),
//   e.g. Mac users who want to use x11 instead of aqua terminal:
// Gnuplot::set_terminal_std("x11");

	cout << "*** example of gnuplot control through C++ ***" << endl << endl;

	//
	// Using the GnuplotException class
	//
	try
	{
		Gnuplot g1("lines");

		//
		// Slopes
		//
		//cout << "*** plotting slopes" << endl;
		//g1.set_title("Slopes\\nNew Line");

		//cout << "y = x" << endl;
		//g1.plot_slope(1.0, 0.0, "y=x");

		//cout << "y = 2*x" << endl;
		//g1.plot_slope(2.0, 0.0, "y=2x");

		//cout << "y = -x" << endl;
		//g1.plot_slope(-1.0, 0.0, "y=-x");
		//g1.unset_title();

		////
		//// Equations
		////
		//g1.reset_plot();
		//cout << endl << endl << "*** various equations" << endl;

		//cout << "y = sin(x)" << endl;
		//g1.plot_equation("sin(x)", "sine");

		//cout << "y = log(x)" << endl;
		//g1.plot_equation("log(x)", "logarithm");

		//cout << "y = sin(x) * cos(2*x)" << endl;
		//g1.plot_equation("sin(x)*cos(2*x)", "sine product");

		////
		//// Styles
		////
		//g1.reset_plot();
		//cout << endl << endl << "*** showing styles" << endl;

		//cout << "sine in points" << endl;
		//g1.set_pointsize(0.8).set_style("points");
		//g1.plot_equation("sin(x)", "points");

		//cout << "sine in impulses" << endl;
		//g1.set_style("impulses");
		//g1.plot_equation("sin(x)", "impulses");

		//cout << "sine in steps" << endl;
		//g1.set_style("steps");
		//g1.plot_equation("sin(x)", "steps");

		////
		//// Save to ps
		////
		//g1.reset_all();
		//cout << endl << endl << "*** save to ps " << endl;

		//cout << "y = sin(x) saved to test_output.ps in working directory" << endl;
		////      g1.savetops("test_output");
		////g1.savetofigure("test_output.ps", "postscript color");
		//g1.set_style("lines").set_samples(300).set_xrange(0, 5);
		//g1.plot_equation("sin(12*x)*exp(-x)").plot_equation("exp(-x)");

		//g1.showonscreen(); // window output


		//
		// User defined 1d, 2d and 3d point sets
		//
		//std::vector<double> x, y, y2, dy, z;
		//int NPOINTS = 50;
		//int SLEEP_LGTH = 2;
		//for (unsigned int i = 0; i < NPOINTS; i++)  // fill double arrays x, y, z
		//{
		//	x.push_back((double)i);             // x[i] = i
		//	y.push_back((double)i * (double)i); // y[i] = i^2
		//	z.push_back(x[i] * y[i]);           // z[i] = x[i]*y[i] = i^3
		//	dy.push_back((double)i * (double)i / (double)10); // dy[i] = i^2 / 10
		//}
		//y2.push_back(0.00);
		//y2.push_back(0.78);
		//y2.push_back(0.97);
		//y2.push_back(0.43);
		//y2.push_back(-0.44);
		//y2.push_back(-0.98);
		//y2.push_back(-0.77);
		//y2.push_back(0.02);


		//g1.reset_all();
		//cout << endl << endl << "*** user-defined lists of doubles" << endl;
		//g1.set_style("impulses").plot_x(y, "user-defined doubles");

		//g1.reset_plot();
		cout << endl << endl << "*** user-defined lists of points (x,y)" << endl;
		g1.set_grid();
		g1.set_style("points").plot_xy(x, y, "user-defined points 2d");

		//g1.reset_plot();
		//cout << endl << endl << "*** user-defined lists of points (x,y,z)" << endl;
		//g1.unset_grid();
		//g1.plot_xyz(x, y, z, "user-defined points 3d");

		//g1.reset_plot();
		//cout << endl << endl << "*** user-defined lists of points (x,y,dy)" << endl;
		//g1.plot_xy_err(x, y, dy, "user-defined points 2d with errorbars");


		////
		//// Multiple output screens
		////
		//cout << endl << endl;
		//cout << "*** multiple output windows" << endl;

		//g1.reset_plot();
		//g1.set_style("lines");
		//cout << "window 1: sin(x)" << endl;
		//g1.set_grid().set_samples(600).set_xrange(0, 300);
		//g1.plot_equation("sin(x)+sin(x*1.1)");

		//g1.set_xautoscale().replot();

		//Gnuplot g2;
		//cout << "window 2: user defined points" << endl;
		//g2.plot_x(y2, "points");
		//g2.set_smooth().plot_x(y2, "cspline");
		//g2.set_smooth("bezier").plot_x(y2, "bezier");
		//g2.unset_smooth();

		//Gnuplot g3("lines");
		//cout << "window 3: log(x)/x" << endl;
		//g3.set_grid();
		//g3.plot_equation("log(x)/x", "log(x)/x");

		//Gnuplot g4("lines");
		//cout << "window 4: splot x*x+y*y" << endl;
		//g4.set_zrange(0, 100);
		//g4.set_xlabel("x-axis").set_ylabel("y-axis").set_zlabel("z-axis");
		//g4.plot_equation3d("x*x+y*y");

		//Gnuplot g5("lines");
		//cout << "window 5: splot with hidden3d" << endl;
		//g5.set_isosamples(25).set_hidden3d();
		//g5.plot_equation3d("x*y*y");

		//Gnuplot g6("lines");
		//cout << "window 6: splot with contour" << endl;
		//g6.set_isosamples(60).set_contour();
		//g6.unset_surface().plot_equation3d("sin(x)*sin(y)+4");

		//g6.set_surface().replot();

		//Gnuplot g7("lines");
		//cout << "window 7: set_samples" << endl;
		//g7.set_xrange(-30, 20).set_samples(40);
		//g7.plot_equation("besj0(x)*0.12e1").plot_equation("(x**besj0(x))-2.5");

		//g7.set_samples(400).replot();

		//Gnuplot g8("filledcurves");
		//cout << "window 8: filledcurves" << endl;
		//g8.set_legend("outside right top").set_xrange(-5, 5);
		//g8.plot_equation("x*x").plot_equation("-x*x+4");

		////
		//// Plot an image
		////
		//Gnuplot g9;
		//cout << "window 9: plot_image" << endl;
		//const int unsigned uiWidth = 255U;
		//const int unsigned uiHeight = 255U;
		//g9.set_xrange(0, uiWidth).set_yrange(0, uiHeight).set_cbrange(0, 255);
		//g9.cmd("set palette gray");
		//unsigned char ucPicBuf[uiWidth * uiHeight];
		//// generate a greyscale image
		//for (unsigned int uiIndex = 0; uiIndex < uiHeight * uiWidth; uiIndex++)
		//{
		//	ucPicBuf[uiIndex] = static_cast<unsigned char>(uiIndex % 255U);
		//}
		//g9.plot_image(ucPicBuf, uiWidth, uiHeight, "greyscale");

		//g9.set_pointsize(0.6).unset_legend().plot_slope(0.8, 20);

		//
		// manual control
		//
		//cout << "echo \"plot x; pause =1\" | gnuplot ";
		//Gnuplot g10;
		//cout << "window 10: manual control" << endl;

		////g10.cmd("set samples 400").cmd("plot abs(x)/2 ; pause -1"); // either with cmd()
		//g10 << "plot sqrt(x)";
		//g10.showonscreen();

		//g10 << "replot sqrt(x)" << "replot sqrt(-x)";    // or with <<

		//wait_for_key();
	}
	catch (GnuplotException ge)
	{
		cout << ge.what() << endl;
	}


	cout << endl << "*** end of gnuplot example" << endl;

}


void regress_line(string filename, int m, int s) {
	int i, n;
	// Reading training data file

	ifstream  trainData;
	trainData.open(filename);
	cout << "Opened file: " << filename;


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



	univar_regressor polyRegress = univar_regressor(x, y, m, s);
	Matrix coeff = polyRegress.fit(x, y, n, m, s);

	bool valid_sol = polyRegress.is_valid_solution();
	if (valid_sol)
		cout << "System has valid Solution" << endl;
	else
		cout << "System is ill conditioned" << endl;

	cout << "\nThe values of the solution coefficients are:\n";
	for (int idx = 0; idx <= m; idx++)
		cout << "x^" << idx << "=" << coeff.at(idx, 0) << endl;            // Print the values of x^0,x^1,x^2,x^3,....    
	cout << "\nThe fitted Polynomial is given by:\ny=";
	for (int idx = 0; idx <= m; idx++)
		cout << " + (" << coeff.at(idx, 0) << ")" << "x^" << idx;
	cout << "\n";

	// Running Prediction on training data
	double *y_pred = new double[n];
	for (i = 0; i < n; i++)
		y_pred[i] = polyRegress.predict(x[i], coeff, m);

	// Printing Actual vs. prediction on console
	cout << "\n Dataset 2_a_dataset_1.csv \n" << endl;
	cout << "\n Predictions: \n" << endl;
	cout << "\nx   y_pred    y" << endl;
	cout << "=================" << endl;


	for (i = 0; i < n; i++) {
		cout << x[i] << " " << y_pred[i] << " " << y[i] << endl;
	}


	//Plotting original points vs. fitted line
	std::vector<double> x_vec(x, x + n);
	std::vector<double> y_vec(y, y + n);
	

	//Formulate the equation with the coefficients
	string eqn = "";
	for (i = 0; i < m + 1; i++)
		eqn = eqn + "+" + to_string(coeff.at(i, 0)) + "*(x**" + to_string(i) + ")";

	//Use GNUPLOT
	Gnuplot g1("lines");
	g1.set_grid();
	g1.plot_equation(eqn, "Fitted Line");  //plot the fitted line
	g1.set_style("points").plot_xy(x_vec, y_vec, "Original Points"); //plot the original points
	
}

int regress_2d(const int x1, const int x2) {  //Accepts the indexes of the needed two columns for regression
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

	Matrix x = Matrix(3, n);
	double *y = new double[n];

	for (i = 1; i <= n; ++i) // start at i = 1 to skip header line
	{
		x.set_at(0, i - 1, 1);        // put 1's in x0 for the algorithm to work, data points will be stored in x1, x2
		x.set_at(1, i - 1, stod(parsedCsv[i][x1]));
		x.set_at(2, i - 1, stod(parsedCsv[i][x2]));
		y[i - 1] = stod(parsedCsv[i][4]);
	}


	// Traget column "Price" is in a very large scale
	//Normalizing "Price"
	double *y_norm = new double[n];
	y_norm = norm_y(y, n);

	Matrix x_norm(m + 1, m + 1);
	x_norm = norm_x(x, n);


	multivar_regressor polyRegress2D = multivar_regressor(x_norm, y_norm, n, m);
	Matrix coeff = polyRegress2D.fit(x_norm, y_norm, n, m);

	cout << "\nThe values of the solution coefficients are:\n";
	for (int idx = 0; idx <= m; idx++)
		cout << "x_" << idx << "=" << coeff.at(idx, 0) << endl;            // Print the values of x^0,x^1,x^2,x^3,....    
	cout << "\nThe fitted Polynomial is given by:\ny=";
	for (int idx = 0; idx <= m; idx++)
		cout << " + (" << coeff.at(idx, 0) << ")" << "x_" << idx;
	cout << "\n";

	//Getting Min, Max price for de-normalizing prediction
	double y_min = y[0];
	double y_max = y[0];
	for (i = 0; i < n; i++) {
		if (y[i] < y_min)
			y_min = y[i];
		if (y[i] > y_max)
			y_max = y[i];
	}

	// Prediction
	double* xi = new double[3];
	double *y_pred = new double[n];
	double* y_norm_pred = new double[n];
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < 3; j++)
			xi[j] = x_norm.at(j, i);
		y_norm_pred[i] = polyRegress2D.predict(xi, coeff, m);
	}

	for (i = 0; i < n; i++)
	{
		y_pred[i] = denorm_y(y_norm_pred[i], y_min, y_max);
	}

	double mse = calc_mse(y_norm, y_norm_pred,n);
	cout << "\n\n=========================================" << endl;
	cout << "Mean Square Error= " << mse << endl;
	cout << "\n\n=========================================" << endl;

	cout << "\n Sample Predictions: \n" << endl;
	cout << "\nx1 x2  Pred_Price Price" << endl;
	cout << "=========================================" << endl;

	for (i = 0; i < 11; i++) {
		cout << x.at(1, i) << "  " << x.at(2, i) << "       " << y_pred[i] << "     " << y[i] << endl;
	}
	cout << "=========================================" << endl;
	cout << "=========================================" << endl;
	return 0;
}


void interpolate(string filename) {

	//Read input file:

	int i, n;
	// Reading training data file

	ifstream  trainData;
	trainData.open(filename);
	cout << "Opened file: " << filename;


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

	cout << "\nNewton Interpolation" << endl;
	Newton_interpolator ni(x, y, n);
	ni.fit();
	// Interpolating 2 points and 4 points
	double* ip2 = new double[2];
	double* ip4 = new double[4];

	for (int p=1; p<3; p++)
		ip2[p] = ni.interpolate(p);
	for (int p = 2; p < 6; p++)
		ip4[p] = ni.interpolate(p);
	//cout << "\nanswer: " << answer;


	// spline interpolation:
	cout << "\nSpline Interpolation: " << endl;
	Spline_interpolator s(x, y, n);
	s.fit();
	double answer2 = s.interpolate(5);
	cout << "\nanswer: " << answer2;
}



int main()
{


	int number = 0;
	cout << "Chose what do you want to do: \n \t1: Single variable Polynomial Regression. \n\t2: 2D Ploynomial Regression. \n\t3: Interpolation. \n";
	cin >> number;

	if (number == 1)
	{
		int m;
		cout << "Linear Regression" << endl;
		cout << "\nWhat degree of Polynomial do you want to use for the fit?\n";
		cin >> m;

		cout << "\n Chose to Sove using: \n\t1.Gauss Scaled Elimination. \n\t2.Iterative Gauss-Siedel\n";
		int s;
		cin >> s;

		cout << "===========================" << endl;
			cout << "\n\nUsing First Dataset:" << endl;
			cout << "===========================" << endl;

			regress_line("2_a_dataset_1.csv", m, s);

			cout << "\n===========================" << endl;
			cout << "\n\nUsing Second Dataset:" << endl;
			cout << "===========================" << endl;

			regress_line("2_a_dataset_2.csv", m, s);
	

	}
	else if (number == 2) {
		cout << "2D Ploynomial Regression" << endl;
		cout << "///////////////////////////////////////" << endl;
		cout << "\nUsing \"bedrooms\" and \"bathrooms\"" << endl;
		int out = regress_2d(0, 1);
		cout << "\nUsing \"bedrooms\" and \"stories\"" << endl;
		out = regress_2d(0, 2);
		cout << "\nUsing \"bedrooms\" and \"lotsize\"" << endl;
		out = regress_2d(0, 3);
		cout << "\nUsing \"bathrooms\" and \"stories\"" << endl;
		out = regress_2d(1, 2);
		cout << "\nUsing \"bathrooms\" and \"lotsize\"" << endl;
		out = regress_2d(1, 3);
		cout << "\nUsing \"stories\" and \"lotsize\"" << endl;
		out = regress_2d(2, 3);

	}
	else if (number == 3) {
		cout << "Interpolation" << endl;

		cout << "\n===========================" << endl;
		cout << "\n\nUsing First Dataset:" << endl;
		cout << "===========================" << endl;

		interpolate("3_dataset_1.csv");

		cout << "===========================" << endl;
		cout << "\n\nUsing Second Dataset:" << endl;
		cout << "===========================" << endl;
		
		interpolate("3_dataset_2.csv");

	}
	else {
		cout << "Wrong Choice !" << endl;
	}
	int ss;
	cin >> ss;
	return 0;


}

