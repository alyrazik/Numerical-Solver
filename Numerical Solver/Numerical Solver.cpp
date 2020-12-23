// Numerical Solver.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
using namespace std;
#include "Linear_system.h"
#include "Matrix.h"

int main()
{
    double temp[20] = { 1, 1, -3, 1, -2, -3, 3, -4, 1, 0, 1, 0, 2, -1, -1, 1, 2, 0, 0, -12};
    //double temp[12] = { 1,1,1,-6,1,2,-1,-2,2,-0.5,1,-4 };
    double t[3] = { 1, 2, 3 };
    Matrix tt(3, 1, t); //always pass list of rows, not columns.
    //Matrix x(3, 4, temp);
    //cout << x;
    //cout << " and the b matrix is" << endl;
    //cout << tt; 
    //Linear_system A(3, 3, x, tt);
    Linear_system B(4, 5, temp);
    Matrix x = B.solve();
    cout << x;
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
