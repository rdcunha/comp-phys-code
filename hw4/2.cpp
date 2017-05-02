#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

//defining the number of timesteps
const int ntot=1000;

//defining other necessary variables
vector <double> y_1(ntot);
vector <double> der1(ntot);
vector <double> g;
vector <double> E;
double r[4];
double s[4];
double tau=1.0/ntot;

//initializing y_1, first derivative of y_1 and the time array
void initialize(double val)
{
	y_1[0]=1.0;
	der1[0]=val;
	g.push_back(val);
	cout<<y_1[0]<<" "<<der1[0]<<'\n';
}

//given y_1, the first derivative and the time, returns the value of the second derivative
double eqn(double y_1)
{
	double val=0;
	val=-4.0*(pow(acos(-1),2.0)*y_1);
	//cout<<val<<'\n';
	return val;
}

//implementing the fourth order runge-kutta to find y_1 and the first derivative of y_1
//finds coefficients r1-4 for y_1 and s1-4 for first derivative of y_1
double rk()
{
	double valu=0;
	for(int k=0;k<ntot;k++)
	{
		r[0]=(tau*der1[k]);
		s[0]=(tau*eqn(y_1[k]));
		//cout<<r[0]<<s[0]<<'\n';
		r[1]=tau*(der1[k]+(s[0]/2.0));
		s[1]=tau*eqn(y_1[k]+(r[0]/2.0));
		//cout<<r[1]<<s[1]<<'\n';
		r[2]=tau*(der1[k]+(s[1]/2.0));
		s[2]=tau*eqn(y_1[k]+(r[1]/2.0));
		//cout<<r[2]<<s[2]<<'\n';
		r[3]=tau*(der1[k]+s[2]);
		s[3]=tau*eqn(y_1[k]+(r[2]));
		der1[k+1]=der1[k]+(s[0]+2*s[1]+2*s[2]+s[3])*(1.0/6.0);
		y_1[k+1]=y_1[k]+(r[0]+2*r[1]+2*r[2]+r[3])*(1.0/6.0);
		//cout<<y_1[k+1]<<" "<<der1[k+1]<<'\n';
	}
	valu=der1[ntot];
	return valu;
}

//using the secant method to obtain the next guess
double check(int j)
{
	double val=0;	
	val=g[j]-E[j]*((g[j-1]-g[j])/(E[j-1]-E[j]));
	return val;
}
	
int main()
{
	initialize(6.0);
	E.push_back(rk()-2*acos(-1));
	cout<<"RK 1 done.\n";
	cout<<"Error: "<<E[0]<<'\n';
	initialize(7.0);
	E.push_back(rk()-2*acos(-1));
	cout<<"RK 2 done.\n";
	cout<<"Error: "<<E[1]<<'\n';
	int i=1;
	//iterating over RK methods for y and y', to obtain convergence for the errors
	while(abs(E[i-1]-E[i])>pow(10,-20))
	{
		initialize(check(i));
		E.push_back(rk()-2*acos(-1));
		cout<<"RK "<<i+2<<" done.\n";
		cout<<"Error: "<<E[i+1]<<'\n';
		i++;
	}	
	cout<<"Final derivative value at 1 : "<<der1[ntot]<<'\n';
	cout<<"Function and derivative values stored in files Guess*.dat\n";
//writing y_1 and the first derivative to file
	ofstream data1("Guess67.dat");
	ofstream data2("Guess67_1.dat");
        for(int i=0;i<=ntot;i++)
        {
                data1.precision(6);
                data1<<y_1[i]<<" "<<der1[i]<<'\n';
        }
        for(int i=0;i<=ntot;i++)
        {
                data2.precision(6);
                data2<<double(i)/ntot<<" "<<y_1[i]<<'\n';
        }
        data1.close();
        data2.close();
	return 0;
}
