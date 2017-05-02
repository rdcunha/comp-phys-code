#include "function.h"
#include <vector>
#include <cmath>
#include <iostream>

using namespace std;

double Function::function(double y)
{
	double val=y*y*y*cos(y);
	return val;
}

//function to use the recursive formulas to get C and D arrays
void Function::makeCD(double y)
{
	for(int m=0;m<(nval-1);++m)
	{
		for(int i=0;i<(nval-m);++i)
		{
			C[m+1][i]=((y-func[i][0])*(C[m][i+1]-D[m][i] ))/(func[i+m+1][0]-func[i][0]);
			D[m+1][i]=((y-func[i+m+1][0])*(D[m][i]-C[m][i+1] ))/(func[i][0]-func[i+m+1][0]);
		}
	}
}

//function to return the value needed
double Function::value(double z)
{
	makeCD(z);
	double val=func[0][1]+C[1][0]+C[2][0]+C[3][0]+C[4][0];
	return val;
}

double Function::derive1(int j)
{
	double der=(func[j+1]-func[j-1])/(2*h);
	//cout<<func[j+1]-func[j-1]<<'\n';
	return der;
}

double Function::derive2(int k)
{
	double dder=(func[k+1]-func[k-1]-(2*func[k]))/(h*h);
	return dder;
}

Function::Function(){}
Function::~Function(){}
