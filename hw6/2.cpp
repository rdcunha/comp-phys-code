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
Matrix A;
vector <double> x;
vector <double> t;
Eigen::VectorXd phi_n;
Eigen::VectorXd b;
Eigen::VectorXd phi_a;
string file1 ="phi_0.005.dat";
string file2 ="phia_0.005.dat";

//the time we are interested in evaluating at
int tim=5;

//declaring necessary functions
double phi_t0(double x);
void initialize();
void analytical(int j);

int main()
{
	phi=phi.Zero(ni+2,nj+2);
	A=A.Zero(ni+2,nj+2);
	//cout<<phi<<'\n';
	initialize();
	//calculating phi
	phi_n.resize(ni+2);
	b=phi.col(0);
	//cout<<b<<'\n';
	for(int j=1;j<nj+2;j++)
	{
		phi_n=A.fullPivLu().solve(b);
		for(int i=1;i<ni+1;i++)
		{
			phi(i,j)=phi_n(i);
			b(i)=phi_n(i);
		}
	}
	analytical(tim);
	ofstream data1(file1);
	ofstream data2(file2);
	data1.precision(6);
	data2.precision(6);
	for(int i=0;i<ni+2;i++)
	{
		data1<<i<<'\t'<<phi(i,tim)<<'\n';
		data2<<i<<'\t'<<phi_a(i)-phi(i,tim)<<'\n';
	}
	//cout<<phi<<'\n';
        data1.close();
        data2.close();
	return 0;
}

//initializing x, t, phi and A
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
	double l=dt/(dx*dx);
	A(0,0) =1.0-(2.0*l);
	for(int i=1;i<ni+1;i++)
	{
		A(i+1,i)=-l;
		A(i-1,i)=-l;
		A(i,i)=1.0+(2.0*l);
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

void analytical(int j)
{
	double fun=0;
	phi_a.resize(ni+2);
	for(int i=0;i<ni+2;i++)
	{
		fun=0;
		for(int n=0;n<2000;n++)
		{
			fun += std::pow(-1.0,n)*std::pow(((2.0*n)+1.0)*std::acos(-1.0),-2.0)*std::sin((n+std::acos(-1.0)/2.0)*x[i])*std::exp(-(std::pow(((2.0*n)+1.0)*std::acos(-1.0),2.0)*t[j]/4.0));
		}
		fun*=8.0;
		phi_a(i)=fun;
	}
}
