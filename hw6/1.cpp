#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Core>

using std::cin;
using std::cout;
using std::vector;
using std::ofstream;
using std::string;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

//ni=number of x grid points
//nj=number of t grid points
double ni=50;
double dx=2.0/ni;
double dt=0.001;
double nj=100;

//defining necessary variables
Matrix phi;
vector <double> x;
vector <double> t;
string file1 ="phi_0.005.dat";

//declaring necessary functions
double phi_t0(double x);
void initialize();

int main()
{
	phi=phi.Zero(ni+2,nj+2);
	//cout<<phi<<'\n';
	initialize();
	//calculating phi
	for(int j=0;j<nj;j++)
	{
		for(int i=1;i<ni+1;i++)
		{
			phi(i,j+1)=phi(i,j)+(dt/(dx*dx))*(phi(i+1,j)+phi(i-1,j)-2*phi(i,j));
		}
	}
	ofstream data1(file1);
	data1.precision(6);
	for(int i=0;i<ni+2;i++)
	{
		data1<<i<<'\t'<<phi(i,5)<<'\n';
	}
	//data1<<phi<<'\n';
        data1.close();
	return 0;
}

//initializing x, t and phi
void initialize()
{
	x.resize(ni+2);
	x[0]=0;
	for(int i=0;i<ni+2;i++)
	{
		x[i+1]=x[i]+dx;
	}
	t.resize(nj+2);
	t[0]=0;
	for(int j=0;j<nj+2;j++)
	{
		t[j+1]=t[j]+dt;
	}
	for(int i=0;i<ni+2;i++)
	{
		phi(i,0)=phi_t0(x[i]);
	}
}

double phi_t0(double x)
{
	double val;
	if(x>=0 && x<=1)
	{
		val=x;
	}
	else if(x>1 && x<=2)
	{
		val=-x+2.0;
	}
	return val;
}
