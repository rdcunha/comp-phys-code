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
//nj=number of y grid points
int ni=20;
double h=1.0/ni;
int nj=20000;

//defining necessary variables
Matrix phi;
Matrix phi_old;
vector <double> x;
vector <double> y;
vector <double> E;
double w=0.5;
string file1 ="phi_w_0.5_h_0.075.dat";
string file2 ="E_w_0.5_h_0.075.dat";


//declaring necessary functions
double phi_y0(double x);
void initialize();
bool converge();

int main()
{
	phi=phi.Zero(ni+1,nj+1);
	//cout<<phi<<'\n';
	initialize();
	int iter=0;
	phi_old=phi_old.Zero(ni+1,nj+1);
	while(converge())
	{
	//calculating phi
		phi_old=phi;
		for(int i=1;i<ni;i++)
		{
			for(int j=1;j<nj;j++)
			{
				phi(i,j)=(1.0-w)*phi_old(i,j)+(w/4.0)*(phi(i+1,j)+phi(i-1,j)+phi(i,j+1)+phi(i,j-1));
			}
		}
	//calculating the Energy functional
		double Eval=0;
		for(int i=0;i<ni;i++)
		{
			for(int j=0;j<nj;j++)
			{
				Eval+=std::pow((phi(i+1,j)-phi(i,j))/h,2.0)+std::pow((phi(i,j+1)-phi(i,j))/h,2.0);
			}
		}
		Eval=Eval*h*h/2.0;
		E.push_back(Eval);
		iter++;
		cout<<"Iteration: "<<iter<<'\n';
	}
	ofstream data1(file1);
	ofstream data2(file2);
	data1.precision(6);
	data2.precision(6);
	data1<<phi<<'\n';
	//Eigen::Map<Matrix,0,Eigen::Stride<Eigen::Dynamic,10> > M4(phi.data(), phi.rows(), (phi.cols()+9)/10,
      	//Eigen::Stride<Eigen::Dynamic,10>(phi.outerStride(),10));
	//data1<< M4 << "\n";
	for(int i=0;i<iter;i++)
	{
		data2<<E[i]<<'\n';
	}
        data1.close();
	return 0;
}

//initializing x, y and phi
void initialize()
{
	x.resize(ni);
	x[0]=0;
	for(int i=0;i<ni;i++)
	{
		x[i+1]=x[i]+h;
	}
	y.resize(nj);
	y[0]=0;
	for(int j=0;j<nj;j++)
	{
		y[j+1]=y[j]+h;
	}
	for(int i=0;i<ni+1;i++)
	{
		phi(i,0)=phi_y0(x[i]);
	}
}

double phi_y0(double x)
{
	double val=x*(1.0-x);
	return val;
}

//checking for convergence
bool converge()
{
	bool check=true;
	double del=0;
	for(int i=0;i<ni+1;i++)
	{
		for(int j=0;j<nj+1;j++)
		{
			del+=std::abs(phi(i,j)-phi_old(i,j));
		}
	}
	del*=(1.0/(ni+1.0))*(1.0/(nj+1.0));
	if(del<pow(10,-7))
	{
		check=false;
	}
	cout<<"Difference: "<<del<<'\n';
	return check;
}
