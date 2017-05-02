#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

using std::cin;
using std::cout;
using std::vector;
using std::ofstream;

//defining the number of timesteps
const int ntot=1000;

//defining other necessary variables
vector <double> theta(ntot);
vector <double> der1(ntot);
vector <double> tim(ntot);
double r[4]={0,0,0,0};
double s[4]={0,0,0,0};
double tau=(3.0*acos(-1))/100.0;
double l=0.3;
double q=0.5;
double factor=2.0/3.0;
double omega0=factor;

//defining the driving force value
double fd=0.89;
double b=fd;

//int count;

//initializing theta, first derivative of theta and the time array
void initialize()
{
	theta[0]=0;
	der1[0]=11.431*sqrt(0.3/9.8);
	tim[0]=0;
	for(int i=0; i<ntot;i++)
	{
		tim[i+1]=tim[i]+tau;
	}
	//cout<<theta[0]<<" "<<der1[0]<<tim[9]<<'\n';
}

//given theta, the first derivative and the time, returns the value of the second derivative
double eqn(double th,double der,double ti)
{
	double val=0;
	val=-sin(th)-(q*der)+(b*cos(omega0*ti));
	//cout<<val<<'\n';
	//count++;
	return val;
}

//implementing the fourth order runge-kutta to find theta and the first derivative of theta
//finds coefficients r1-4 for theta and s1-4 for first derivative of theta
void rk()
{
	for(int k=0;k<ntot;k++)
	{
		r[0]=(tau*der1[k]);
		s[0]=(tau*eqn(theta[k],der1[k],tim[k]));
		//cout<<r[0]<<s[0]<<'\n';
		r[1]=tau*(der1[k]+(s[0]/2.0));
		s[1]=tau*eqn(theta[k]+(r[0]/2.0),der1[k]+(s[0]/2.0),tim[k]+tau/2.0);
		//cout<<r[1]<<s[1]<<'\n';
		r[2]=tau*(der1[k]+(s[1]/2.0));
		s[2]=tau*eqn(theta[k]+(r[1]/2.0),der1[k]+(s[1]/2.0),tim[k]+tau/2.0);
		//cout<<r[2]<<s[2]<<'\n';
		r[3]=tau*(der1[k]+s[2]);
		s[3]=tau*eqn(theta[k]+(r[2]),der1[k]+(s[2]),tim[k]+tau);
		der1[k+1]=der1[k]+(s[0]+2*s[1]+2*s[2]+s[3])*(1.0/6.0);
		theta[k+1]=theta[k]+(r[0]+2*r[1]+2*r[2]+r[3])*(1.0/6.0);
		cout<<theta[k+1]<<" "<<der1[k+1]<<'\n';
	}
}

int main()
{
	initialize();
	rk();
	/*for(int j=0;j<ntot;j++)
	{
	cout<<tim[j]<<'\n';
	}*/

//writing theta and the first derivative to file
	ofstream data1("fd_89.dat");
        for(int i=0;i<ntot;i++)
        {
                data1.precision(6);
                data1<<theta[i]<<" "<<der1[i]<<'\n';
        }
        data1.close();
	return 0;
}
