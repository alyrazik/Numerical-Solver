#include "Matrix.h"
#include <iostream>

Matrix::Matrix(const int& x, const int& y)
    :n(x)
    ,m(y)
    ,array(new double[(x) * (y)])
{
}

Matrix::Matrix(const Matrix& m2)
{
    n = m2.n;
    m = m2.m;
    array = new double[n*m];
    for (int i = 0; i < n * m; i++)
        array[i] = m2.array[i];
}

Matrix::~Matrix()
{
    delete array;
}

int Matrix::n_rows() const
{
    return n;
}

int Matrix::n_cols() const
{
    return m;
}

Matrix Matrix::augment(const Matrix& y) const
{
    Matrix temp(n, m + y.n_cols());
    for (int k = 0; k < n; k++)
    {
        for (int i = 0; i < m; i++)
            temp.set_at(k, i, this->at(k, i));
        for (int j = 0; j < y.n_cols(); j++)
            temp.set_at(k, j+m, y.at(k, j));
    }
    //std::cout << temp;
    return temp;
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
