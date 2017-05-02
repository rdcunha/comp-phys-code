#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <Eigen/Dense>
#include <Eigen/Core>

using namespace std;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::VectorXf Vector;

//########################defining necessary variables##########################//

//l = size of the box
int l=9;
//dim = dimension of the box
int dim=2;

//ni = number of x grid points
//nj = number of t grid points
int ni=l-1;
int nj=l-1;
int ntot=ni*nj;

//step_size  = step size;
//nsteps = total # of steps;
int nsteps=60000;
double step_size=0.01;
double time_tot=nsteps*step_size;

double dis=0;
double dx=0;
double dy=0;

double temp=0;
double P0=0;
double P1=0;
double P10=0;
double P11=0;

//important matrices
Matrix R(ntot,dim);
Matrix V(ntot,dim);
Matrix G_old(ntot,dim);
Matrix G_new(ntot,dim);
Vector E(nsteps);
Vector KE(nsteps);
Vector PDF(450);
Vector Prob(32);

//#######################declaring necessary functions##########################//
void initialize();
void compute_r(Matrix &mat,double t);
void compute_E(Vector &vec,int n);
void compute_g(Matrix &mat);
void compute_v(Matrix &mat,double t);
void conserve_p(Matrix &mat);
double compute_temp(Vector &vec,int m,int n);
void compute_P(double &sump0, double &sump1);
void compute_gr(Vector &vec);
void compute_prob(Vector &vec);

int main()
{
	initialize();
	//cout<<endl<<R<<endl;
	//cout<<V<<endl;
	compute_g(G_old);
	//cout<<G_old<<endl;
	compute_E(E,0);
	//cout<<"E : "<<E(0)<<endl;
	int i=0;
	ofstream energy("energy.dat");
	ofstream gr("pdf.dat");
	ofstream prob("prob.dat");
	energy.precision(6);
	gr.precision(6);
	prob.precision(6);
	PDF.setZero(450);
	cout<<"Beginning MD simulation."<<endl;
	for(int i=0;i<nsteps-1;i++)
	{
		compute_r(R,step_size);
		compute_g(G_new);
		//cout<<"G : "<<G_new<<endl;
		compute_v(V,step_size);
		G_old=G_new;
		//cout<<"E : "<<E(i)<<endl;
		compute_E(E,i+1);
		if(i%200==0)
		{
			conserve_p(V);
		}
		compute_P(P0,P1);
		//cout<<"P at "<<i<<" :"<<P0<<endl;
		if(i>nsteps-10000)
		{
			compute_P(P10,P11);
		}
		energy<<step_size*i<<'\t'<<E(i)<<endl;
		int ir=0;
		compute_gr(PDF);
		if(i>nsteps-1000)
		{
			compute_prob(Prob);
		}
		if(i%(nsteps/10)==0)
			cout<<100*i/(nsteps)<<"% done"<<endl;
	}
	cout<<"MD simulation finished."<<endl;
	Prob/=(1000*ntot);
	cout<<"Probability density computed. Storing in file 'prob.dat'."<<endl;
	prob<<Prob<<endl;
	PDF/=nsteps;
	PDF*=(double)(l*l);
	PDF/=ntot*ntot;
	cout<<"Pair distribution function computed. Storing in file 'pdf.dat'."<<endl;
	gr<<PDF<<endl;
	prob.close();
	gr.close();
	cout<<"Energy computed. Storing in file 'energy.dat'."<<endl;
	energy.close();
	//cout<<"R final : "<<R<<endl;
	temp=compute_temp(KE,0,nsteps);
	cout<<"Avg Temp: "<<temp<<endl;
	cout<<"Avg Temp(10k): "<<compute_temp(KE,nsteps-10000,nsteps)<<endl;
	P0+=P1;
	//cout<<"P"<<P0;
	P0=(1.0-(P0/((double)nsteps*(double)ntot*2.0*temp)))*((double)ntot*temp)/((double)l*(double)l);
	cout<<"Avg Pressure: "<<P0<<endl;
	P10+=P11;
	P10=(1.0-(P10/(nsteps*ntot*2.0*temp)))*(ntot*temp)/((double)l*(double)l);
	cout<<"Avg Pressure(10k): "<<P10<<endl;
	return 0;
}

//#######################defining necessary functions##########################//
void initialize()
{
	random_device rd;				//random number generation
	mt19937 mt(rd());
	uniform_real_distribution<double> dist(-0.5,0.5);
	double sum0=0;
	double sum1=0;
	for(int i=0;i<ntot;i++)
	{
		R(i,0)=i/8;
		R(i,1)=(i%8);
		V(i,0)=dist(mt);
		V(i,1)=dist(mt);
		if(i==ntot-1)
		{
			V(i,0)=-sum0;
			V(i,1)=-sum1;
		}
		else
		{
			sum0+=V(i,0);
			sum1+=V(i,1);
		}
	}
}
			
void compute_g(Matrix &mat)
{
	double sum0=0;
	double sum1=0;
	dx=0;dy=0;mat.setZero();
	for(int i=0;i<ntot;i++)
	{
		sum0=0;
		sum1=0;
		for(int j=0;j<ntot;j++)
		{
			if(j!=i)
			{
				dx=(double)(R(i,0)-R(j,0));
				dy=(double)(R(i,1)-R(j,1));
				if(dx>(double)l/2.0)
				{
					dx-=(double)l;
				}
				if(dx<(double)l*(-1.0)/2.0)
				{
					dx+=(double)l;
				}
				if(dy>(double)l/2.0)
				{
					dy-=(double)l;
				}
				if(dy<(double)l*(-1.0)/2.0)
				{
					dy+=(double)l;
				}
				dis=sqrt(pow(dx,2.0)+pow(dy,2.0));
				//cout<<"dx = "<<dx<<"dy = "<<dy<endl;
				//cout<<"dis "<<i<<" "<<j<<":"<<dis<<endl;
				//cout<<"R"<<i<<"x - R"<<j<<"x :"<<R(i,0)-R(j,0)<<endl;
				//cout<<"R"<<i<<"y - R"<<j<<"y :"<<R(i,1)-R(j,1)<<endl;
				sum0+=(dx)*(2.0*pow(1.0/dis,14)-pow(1.0/dis,8));
				sum1+=(dy)*(2.0*pow(1.0/dis,14)-pow(1.0/dis,8));
			}
		}
		mat(i,0)=sum0*24.0;
		mat(i,1)=sum1*24.0;
	}
}	

void compute_E(Vector &vec,int n)
{
	double sum=0;
	double sumv=0;
	for(int i=0;i<ntot;i++)
	{
		sumv+=pow(V(i,0),2.0)+pow(V(i,1),2.0);
		for(int j=0;j<ntot;j++)
		{
			if(j!=i)
			{
				dx=(double)(R(i,0)-R(j,0));
				dy=(double)(R(i,1)-R(j,1));
				if(dx>(double)l/2.0)
				{
					dx-=(double)l;
				}
				if(dx<(double)l*(-1.0)/2.0)
				{
					dx+=(double)l;
				}
				if(dy>(double)l/2.0)
				{
					dy-=(double)l;
				}
				if(dy<(double)l*(-1.0)/2.0)
				{
					dy+=(double)l;
				}
				dis=sqrt(pow(dx,2.0)+pow(dy,2.0));
				sum+=(pow(1.0/dis,12)-pow(1.0/dis,6));
			}
		}
	}
	sum*=4.0;
	sumv/=2.0;
	KE(n)=sumv;
	//cout<<"sumv : "<<sumv<<" and sum : "<<sum<<endl;
	E(n)=sumv+sum;
}
	 

void compute_r(Matrix &mat,double t)
{
	mat=mat+(t*V)+(pow(t,2.0)*G_old)*0.5;
	for(int i=0;i<ntot;i++)
	{
		if(R(i,0)>l) { R(i,0)-=l;}
		if(R(i,0)<0) { R(i,0)+=l;}
		if(R(i,1)>l) { R(i,1)-=l;}
		if(R(i,1)<0) { R(i,1)+=l;}
	}
}	

void compute_v(Matrix &mat, double t)
{
	mat=mat+t*(G_old+G_new)*0.5;
}

void conserve_p(Matrix &mat)
{
	double sum0=0;
	double sum1=0;
	for(int i=0;i<ntot-1;i++)
	{
		sum0+=mat(i,0);
		sum1+=mat(i,1);
	}
	mat(ntot-1,0)=-sum0;
	mat(ntot-1,1)=-sum1;
}

double compute_temp(Vector &vec,int m,int n)
{
	double val=0;
	for(int i=m;i<n;i++)
	{
		val+=KE(i);
	}
	val/=n;
	val/=ntot;
	return val;
}

void compute_P(double &sump0,double &sump1)
{
	dx=0;dy=0;
	for(int i=0;i<ntot;i++)
	{
		for(int j=i+1;j<ntot;j++)
		{
				dx=(double)(R(i,0)-R(j,0));
				dy=(double)(R(i,1)-R(j,1));
				if(dx>(double)l/2.0)
				{
					dx-=(double)l;
				}
				if(dx<(double)l*(-1.0)/2.0)
				{
					dx+=(double)l;
				}
				if(dy>(double)l/2.0)
				{
					dy-=(double)l;
				}
				if(dy<(double)l*(-1.0)/2.0)
				{
					dy+=(double)l;
				}
				dis=sqrt(pow(dx,2.0)+pow(dy,2.0));
				//cout<<"dx = "<<dx<<"dy = "<<dy<endl;
				//cout<<"dis "<<i<<" "<<j<<":"<<dis<<endl;
				//cout<<"R"<<i<<"x - R"<<j<<"x :"<<R(i,0)-R(j,0)<<endl;
				//cout<<"R"<<i<<"y - R"<<j<<"y :"<<R(i,1)-R(j,1)<<endl;
				sump0+=(dx)*(dx)*(2.0*pow(1.0/dis,14)-pow(1.0/dis,8));
				sump1+=(dy)*(dy)*(2.0*pow(1.0/dis,14)-pow(1.0/dis,8));
		}
	}
}	

void compute_gr(Vector &vec)
{
	double dist=0;
	double dr=0.01;
	int index=0;
	for(int i=0;i<ntot;i++)
	{
			for(int j=0;j<ntot;j++)
			{
				dx=(double)(R(i,0)-R(j,0));
				dy=(double)(R(i,1)-R(j,1));
				dist=hypot(dx,dy);
				index=dist/dr;
				//cout<<"works "<<i<<" "<<j<<" "<<endl;	
				if(dist<4.5)
				{
					vec(index)+=1.0/(acos(-1.0)*dist*dr);
				}
			}
	}
}

void compute_prob(Vector &vec)
{
	double nbin=32;
	double dvel=0.16;
	double v;
	int index=0;
	for(int i=0;i<ntot;i++)
	{
		v=hypot(V(i,0),V(i,1));
		index=v/dvel;
		if(index>nbin)
		{
			index=nbin;
		}
		vec(index)+=1.0;
	}
}
