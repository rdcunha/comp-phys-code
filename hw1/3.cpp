#include <cmath>
#include <vector>
#include <iostream>

using std::cin;
using std::cout;
using std::vector;

//defining needed variables and vectors
vector <double> func;
double integ=0;
double maxm = 1;
double minm = 0;
int nslice=8;
double h=(maxm-minm)/nslice;
double x;

//obtaining the value of the function at a point
double function(double y)
{
	double val=exp(-y);
	return val;
}

//using Simpson's rule to numerically integrate
double integrate()
{
	for(int i=0;i<(nslice/2);i++)
	{
		integ+=((func[2*i])+(4*func[(2*i)+1])+(func[(2*i)+2]));
	}
	integ=(h/3)*integ;
	return integ;
}

int main()
{
	x=minm;
	func.push_back(function(x));
//storing the function values in a vector called func
	for(int i=0;i<nslice+1;i++)
	{
		x+=h;
		func.push_back(function(x));
	}
//outputting the value of the integral
	cout<<"The value of the integral over {"<<minm<<", "<<maxm<<"} is "<<integrate()<<'\n';
	return 0;
}
