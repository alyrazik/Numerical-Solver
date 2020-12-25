// Numerical Solver.cpp : This file contains the 'main' function. Program execution begins and ends there.
#include <iostream>
using namespace std;
#include "Linear_system.h"
#include "Matrix.h"

int main()
{
    //double temp[20] = { 1, 1, -3, 1, -2, -3, 3, -4, 1, 0, 1, 0, 2, -1, -1, 1, 2, 0, 0, -12};
    //double t[3] = { 1, 2, 3 };
    /*
    int size, i;
    cin >> size;
    double* as = new double[size];
    for (i = 0; i < size; i++)
    {
        cout << "Enter value of matrix element " <<i<< endl;
        cin >> as[i];
    }
    */
    //double as [16] = { 1, 1, -3, 1, -3, 3, -4, 1, 1, 0, 2, -1, 1, 2, 0, 0 };
    //double bs [4] = { -2, 0, -1, -12 };
    //double as[9] = { 1, 2, 1, 1, 1, -1, 1, 1, 1 }; // sol is 1, 2, 3
    //double bs[3] = { -8, 0, -6 };
    //double as[16] = { 1, 1,1,1, 1, -2, 0.5, 1, 2, -1, 3, 0.5, 1, 1, -1, -1 }; //sol is 1, 2, 3, 4
    //double bs[4] = { -10, -2.5, -11, 4 };
    //double as[9] = { 2, 1, -1, 1, 4, 3, -1, 2, 7 };
    //double bs[3] = { 0, -14, -30 };
    double as[9] = { 3, -0.1, -0.2, 0.1, 7, -0.3, 0.3, -0.2, 10 };
    double bs[3] = { -7.85, 19.3, -71.4 };
    Matrix A(3, 3, as); 
    //double bs[3] = { 0, -14, -30 };
    Matrix b(3, 1, bs);
    Matrix out = A.augment(b); //don't forget to multiply with -1 to the b elements.
    Linear_system x(A, b);
    Matrix sol = x.solve();
    cout << "solution" << endl<< sol;
    //double arr[3] = { 0 };
    //Matrix sol2 = x.solve_iteratively(arr, 30);
    //cout << "solution iteratively" << endl << sol2;
    cout <<"system is "<< x.is_valid_solution();
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

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
