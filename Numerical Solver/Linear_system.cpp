#include "Linear_system.h"

Linear_system(const int* n , const int* m)
{
	A = new double(n * m);
}

~Linear_system()
{
	delete A;
}

Linear_system::Linear_system()
{
	A = new double(1);
}


double* Linear_system::solve()
{

	return vector<double>();
}

