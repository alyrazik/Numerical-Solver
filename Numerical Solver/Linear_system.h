#pragma once
#include <vector>
#include <array>
class Linear_system
{
private:
	bool valid_solution;
	int n;
	int m;
	double* A;
	double* b;

public:
	//default constructor
	Linear_system(); 
	Linear_system(const int* n, const int* m);
	~Linear_system();
	double* solve();

};

