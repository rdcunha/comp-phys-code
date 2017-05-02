#include <cmath>
#include <vector>
#include <iostream>
#include <random>
#include <fstream>

using namespace std;

//using a large number of x values
int m=10000000;
//setting the number of values skipped initially
double n0=1000;
//setting the value of delta
double del=1.0;

//declaring necessary vectors
vector <double> func(m);
vector <double> x;
vector <double> w;
vector <double> c(2);

//declaring necessary variables
double x0=0,x1=0,p=0,xi=0,eta=0;
double expn=0,expnn=0,expn2=0,expnl=0;
double std_dev=0;

//declaring function to create function values
double function(double ex);
//declaring function to create values of W
double dubs(double ex);
//declaring function to find convergence
bool converge();

int main()
{
	random_device rd;				//random number generation
	mt19937 mt(rd());
	uniform_real_distribution<double> dist(0.0,1.0);
	x0=dist(mt);
	cout<<"X0 = "<<x0<<'\n';
	int count =0;
	for(int i=0;count<m;i++)			//using count to keep track of x values stored
	{
		eta=dist(mt);
		//cout<<eta<<'\n';
		x1=x0+del*(2.0*eta-1.0);
		if(x1>=1 || x1<=0)
		{}					//discarding unreasonable values
		else
		{
			//cout<<x1<<'\n';
			p=dubs(x1)/dubs(x0);
			xi=dist(mt);
			//cout<<xi<<'\n';
			if(xi<=p)
			{
				x.push_back(x1);
				x0=x1;
			}
			else
			{
				x.push_back(x0);
			}
			count++;
		}
	}
	/*for(auto i:x)
	{
		cout<<i<<'\n';
	}*/
	c[0]=0;
	c[1]=10;				//setting initial values of C(l)
	
	bool con=converge();
	ofstream cs("C_vals.dat");
	cs.precision(8);
	double cnt=0.0;
	for(int l=1;;l++)
	{
		cnt=0.0;
		expn=0;
		expnn=0;
		expnl=0;
		expn2=0;
		for(int n=n0;n<m;n++)
		{
			//cout<<function(x[n])<<'\t'<<dubs(x[n])<<'\n';
			if(n%l==0)
			{
				expn+=function(x[n])/dubs(x[n]);
				expnl+=function(x[n+l])*function(x[n])/(dubs(x[n+l])*dubs(x[n]));
				expnn+=function(x[n])*function(x[n])/(dubs(x[n])*dubs(x[n]));;
				cnt+=1.0;
			}
		}
		expn2=expn*expn;
		expn2/=cnt*cnt;
		//cout<<expn2<<'\n';
		expnl/=cnt;
		//cout<<expnl<<'\n';
		expnn/=cnt;
		//cout<<expnn<<'\n';
		c[0]=c[1];
		c[1]=((expnl-expn2)/(expnn-expn2));
		std_dev=pow((expnn-expn2)/cnt,0.5);
		con=converge();
		//cout<<con<<'\n';
		if(con)
		{
			cout<<"\nFinal l = "<<l<<'\n';
			cout<<"Final integral value = "<<expn/cnt<<'\n';
			cout<<"Final std Deviation = "<<std_dev<<'\n';
			break;
		}
		cout<<"l = "<<l<<'\n';
		cout<<"C = "<<c[1]<<'\n';
		cout<<"Integral value = "<<expn/cnt<<'\n';
		cout<<"Std Deviation = "<<std_dev<<'\n';
		cs<<l<<'\t'<<c[1]<<'\n';
	}	
	cs.close();
	return 0;
}

double function(double ex)
{
	double fun=0;
	fun=ex*ex;
	return fun;
}

double dubs(double ex)
{
	double dub=(3.0/31.0)*(ex*ex +10.0);
	return dub;
}

bool converge()
{
	bool works=false;
	if(abs(c[1]-c[0])<pow(10,-4.0))
	{
		works=true;
	}
	return works;
}	
