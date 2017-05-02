#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>

using namespace std;

//defining necessary variables
vector <vector <double> > s;
vector <int> ip;
vector <int> im;
vector <double> w;
vector <int> ses;
vector <double> c(2);

//temperatures at which the magnetization, susceptibility are calculated
vector <double> temp={1.6,2.0,2.16,2.34,2.43,2.76,3.35};

//lattice size
int L=8;

double T=0;

//number of MCSS steps to skip initially
int n0=500;

//index for W lookup
int ien=0;

//max number of MCSS steps
int m=100000;

double cnt=0.0,expn=0,expnn=0,expnl=0,expn2=0,expn4=0,cu=0,chi=0,stddev=0;

//declaring functions

//function to print matrix
void print_matrix(vector < vector<double> > vec);
//function to create look up table for indices
void look_up_mod();
//function to create look up table for W
void look_up_table();
void initialize();
//function to do the monte carlo steps
void mc();
//function to sample and find the values of magnetization, susceptibility
void sample();
//function to test convergence of the autocorrelation function
bool converge();

int main()
{
	ofstream data("lattice_8.dat");
	data.precision(8);
	look_up_mod();
	cout<<"Look up tables for indices:\n";
	for(int i=0;i<L;i++)
	{
		cout<<ip[i]<<'\t'<<im[i]<<'\t';
		cout<<'\n';
	}
//looping over all temperatures in the vector
	for(int t=0;t<7;t++)
	{
		T=temp[t];
		cout<<"\nTemperature = "<<T<<endl;
		look_up_table();
		cout<<"Look up table for W:\n";
		for(int i=0;i<5;i++)
		{
			cout<<w[i]<<'\t';
		}
		initialize();
		mc();
		sample();
//printing data to file
		if(t==0)
			data<<"Temp\t\tMagnetization\t\tSusceptibility\t\tCumulant\t\tStd Deviation"<<endl;
		data<<T<<"\t\t"<<expn/cnt<<"\t\t"<<chi<<"\t\t"<<cu<<"\t\t"<<stddev<<endl;
		expn=0;
	}
	data.close();
	return 0;
}

//function to print a matrix of size LxL
void print_matrix(vector < vector<double> > vec)
{
	for(int i=0;i<L;i++)
	{
		for(int j=0;j<L;j++)
		{
			cout<<'\t'<<s[i][j]<<'\t';
		}
		cout<<'\n';
	}
	cout<<'\n';
}

void look_up_mod()
{
	ip.resize(L);
	im.resize(L);
	for(int i=0;i<L;i++)
	{
		ip[i]=i+1;
		im[i]=i-1;
	}
	ip[L-1]=0;
	im[0]=L-1;
}

void look_up_table()
{
	w.resize(5);
	for(int i=0;i<4;i++)
	{
		w[i]=1;
		if(i>2)
		{
			w[i]=exp(-4.0/T);
			w[i+1]=exp(-8.0/T);
		}
	}
}

void initialize()
{
	s.resize(L,vector <double> (L));
	for(int i=0;i<L;i++)
	{
		for(int j=0;j<L;j++)
		{
			s[i][j]=1;
		}
	}
}

void mc()
{
	random_device rd;				//random number generation
	mt19937 mt(rd());
	uniform_real_distribution<double> dist(0.0,1.0);
	for(int mcss=0;mcss<m;mcss++)
	{
		//print_matrix(s);
		int count=0;
		for(int i=0;i<L;i++)
		{
			for(int j=0;j<L;j++)
			{
				ien=2.0*s[i][j]*(s[ip[i]][j]+s[i][ip[j]]+s[im[i]][j]+s[i][im[j]]);
				//cout<<"W index is:"<<ien<<'\n';
				ien+=8;
				ien=ien/4;
				//cout<<"W index is:"<<ien<<'\n';
				//cout<<"W[ien] = "<<w[ien]<<endl;
				if(dist(mt) < w[ien])
				{
					s[i][j]*=-1;
				}
				ien=0;
			}
		}
		//print_matrix(s);
		if(mcss>n0)
		{
			count++;
			int sum=0;
			//cout<<"\nWorks here.\n";
			for(int i=0;i<L;i++)
			{
				for(int j=0;j<L;j++)
				{
					sum+=s[i][j];
				}
			}
			//cout<<"configuration sum is: "<<sum<<'\n';
			ses.push_back(sum);
		}
	}
	//print_matrix(s);
}
			
void sample()
{
	for(int l=1;;l++)
	{
		cnt=0.0;
		expn=0;
		expnn=0;
		expnl=0;
		expn2=0;
		expn4=0;
		//cout<<ses.size()<<endl;
		for(int n=0;n<m;n++)
		{
			if(n%l==0)
			{
				expn+=abs(ses[n]);
				expnl+=abs(ses[n+l])*abs(ses[n]);
				expnn+=abs(ses[n])*abs(ses[n]);
				expn4+=pow(ses[n],4.0);
				cnt+=1.0;
			}
		}
		ses.clear();
		expn/=(double)L*(double)L;
		expn2=expn*expn;
		expn2/=cnt*cnt;
		//cout<<expn2<<'\n';
		expnl/=((double)L*(double)L*(double)L*(double)L);
		expnl/=cnt;
		//cout<<expnl<<'\n';
		expnn/=(cnt*(double)L*(double)L*(double)L*(double)L);
		//cout<<"<m^2>^2 : "<<expnn*expnn<<endl;
		//cout<<"count : "<<cnt<<'\n';
		expn4/=(cnt*(double)L*(double)L*(double)L*(double)L*(double)L*(double)L*(double)L*(double)L);
		//cout<<"<m^4> : "<<expn4<<endl;
		c[0]=c[1];
		c[1]=((expnl-expn2)/(expnn-expn2));
		bool con=converge();
		if(T>2.269)
			chi=expnn*(double)L*(double)L/T;
		else
			chi=(expnn-expn2)*(double)L*(double)L/T;
		cu=1.0 - expn4/(3.0*expnn*expnn);
		stddev=pow((expnn-expn2)/cnt,0.5);
		if(con)
		{
			cout<<"\nFinal l = "<<l<<'\n';
			cout<<"Magnetization per site = "<<expn/cnt<<endl;
			cout<<"Susceptibility = "<<chi<<endl;
			cout<<"Standard deviation = "<<stddev<<endl;
			break;
		}
	}
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
	
