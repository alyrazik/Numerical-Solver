#include "Matrix.h"
#include <iostream>

Matrix::Matrix(const int& x, const int& y)
    :n(x)
    ,m(y)
    ,array(new double[(x) * (y)])
{
}

Matrix::~Matrix()
{
    delete array;
}

double Matrix::at(const int& i, const int& j) const
{
    return array[(i) * m + (j)];
}

void Matrix::set_at(const int& i, const int& j, const double& number)
{
    array[(i)*m + (j)] = number;

}

Matrix::Matrix(const int& x, const int& y, const double a[])
    :n(x)
    ,m(y)
    ,array(new double[(x) * (y)])
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            array[i*m+j] = a[i*m+j];
}

std::ostream& operator<<(std::ostream& out, const Matrix& m)
{
    for (int i = 0; i < m.n; i++)
    {
        for (int j = 0; j < m.m; j++)
            out << m.at(i, j) << " ";
        out << "\n";
    }
    return out;
}
