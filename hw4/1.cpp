#include <iostream>
#include <cmath>
#include <math.h>
#include <vector>

using namespace std;

//setting the dimensions of all vectors
const int n=4;

//declaring necessary vectors and variables
vector< vector<double> > a{{2,3,10,-1},{10,15,3,7},{-4,1,2,9},{0,0,0,0}};
vector< vector<double> > b{{1,0,0,0},{2,0,0,0},{3,0,0,0},{0,0,0,0}};
vector< vector<double> > aa(n, vector<double>(n));
vector< vector<double> > aat(n, vector<double>(n));
vector< vector<double> > rot(n, vector<double>(n));
vector< vector<double> > v(n, vector<double>(n));
vector< vector<double> > vt(n, vector<double>(n));
vector< vector<double> > w(n, vector<double>(n));
vector< vector<double> > u(n, vector<double>(n));
vector< vector<double> > ut(n, vector<double>(n));
vector< vector<double> > x(n, vector<double>(n));

double theta=0;
double t=0,c=0,s=0,tau=0;
int sweep_count=0;

//function to create ATA, AAT, initialize V and U
void initialize()
{
	for(int i=0;i<n;i++)
	{
		for(int k=0;k<n;k++)
		{ 
			for(int j=0;j<n;j++)
			{
				aa[i][k]+=a[j][i]*a[j][k];
				aat[i][k]+=a[i][j]*a[k][j];
			}
			//cout<<aa[i][k]<<" ";
			if(k==i)
			{
				v[i][k]=1;
				u[i][k]=1;
			}
		}
		//cout<<'\n';
	}
}

//function to reset the rotation matrix Ppq for every p,q
void initialize_rot()
{
	for(int i=0;i<n;i++)
	{
		for(int k=0;k<n;k++)
		{
			if(i==k)
				rot[i][k]=1;
			else
				rot[i][k]=0;
		}
	}
}

//function to print a matrix of size nxn
void print_matrix(vector < vector<double> > vec)
{
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			printf("\t%.8f\t",vec[i][j]);
		}
		cout<<'\n';
	}
	cout<<'\n';
}

//function to multiply two given matrices
vector < vector<double> >  multiply(vector < vector<double> > vec1,vector < vector<double> > vec2)
{
	vector <vector<double> > vec(n, vector<double>(n));
	for(int i=0;i<n;i++)
	{
		for(int k=0;k<n;k++)
		{ 
			for(int j=0;j<n;j++)
			{
				vec[i][k]+=vec1[i][j]*vec2[j][k];
			}
		}
		//cout<<'\n';
	}
	return vec;
}

//signum function			
double signum(double num)
{
	if(num<0.0)
		return -1.0;
	if(num>=0.0)
		return 1.0;
}

//function to update the sum of the off-diagonal elements for convergence
double update_s(vector < vector<double> > vec)
{
	double sum=0;
	for(int p=0;p<n;p++)
	{
		for(int q=0;q<n;q++)
		{
			if(q!=p)
				sum+=abs(vec[p][q]);
		}
	}
	return sum;
	//cout<<sum<<'\n';
}

//function to do the jacobi transformation
void jacobi(vector < vector<double> > &veca, vector < vector<double> > &vecb)
{
	double s0=update_s(veca);
	while(update_s(veca) >pow(10,-20))
	{
		sweep_count++;
		for(int p=0;p<n;p++)
		{
			for(int q=p+1;q<n;q++)
			{
				initialize_rot();
			//all variables created here
				theta=(veca[q][q]-veca[p][p])/(2.0*veca[p][q]);
				t = signum(theta)/ (abs(theta) + pow(pow(theta,2.0) + 1.0,0.5));
				c=1.0/sqrt(pow(t,2.0)+1.0);
				s=c*t;
				tau=s/(1.0+c);
			//updating diagonal elements of the vector to be diagonalized
				veca[p][p]=veca[p][p]-t*veca[p][q];
				veca[q][q]=veca[q][q]+t*veca[p][q];
			//updating non-diagonal elements of the vector
				for(int r=0;r<n;r++)
				{
					if(r!=p && r!=q)
					{
						double temp=veca[r][p];
						double temp1=veca[r][q];
						veca[r][p]=temp-(s*(temp1+(tau*temp)));
						veca[p][r]=veca[r][p];
						veca[r][q]=temp1+(s*(temp-(tau*temp1)));
						veca[q][r]=veca[r][q];
					}
				}
				veca[p][q]=0.0;
				veca[q][p]=0.0;
			//creating the rotation matrix for p,q
				rot[p][p]=c;
				rot[q][q]=c;
				rot[p][q]=s;
				rot[q][p]=-s;
			//multiplying to get the eigenvector matrix
				vecb=multiply(vecb,rot);
				//print_matrix(vecb);
				//cout<<"theta: "<<theta<<"\tt: "<<t<<"\tc: "<<c<<"\ts: "<<s<<'\t'<<'\n';
				//print_matrix(veca);
				//cout<<'\n';
				//print_matrix(rot);
				//cout<<'\n';
			}
		}
	}
}

//finding the largest element on the diagonal for sorting
int compare_y(int j, vector < vector<double> > veca)
{
	int largest_y=j;
	for(int i=j;i<n;i++)
	{
		if(abs(veca[i][i])>abs(veca[largest_y][largest_y]))
		{
			largest_y=i;
		}
	}
	return largest_y;
}

//sorting on k
void csort(vector < vector<double> > &veca, vector < vector<double> > &vecb)
{
	for(int k=0;k<n;k++)
	{
		vector <double> temp;
		int l=compare_y(k, veca);
		temp.push_back(veca[k][k]);
		veca[k][k]=veca[l][l];
		veca[l][l]=temp[0];
		temp.clear();
		for(int i=0;i<n;i++)
		{
			temp.push_back(vecb[i][k]);
			vecb[i][k]=vecb[i][l];
			vecb[i][l]=temp[i];
		}
		temp.clear();
	}
}

//function to do the SVD transformation
void svd()
{
	for(int i=0;i<n;i++)
	{
		w[i][i]=sqrt(aa[i][i]);
		if(std::isnan(w[i][i]) == true)
			w[i][i]=0;
	}
	csort(w,v);
	cout<<"W matrix:\n";
	print_matrix(w);
	csort(aat,u);
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			vt[i][j]=v[j][i];
			ut[i][j]=u[j][i];
		}
	}
	cout<<"UWVT:\n";
	print_matrix(multiply(multiply(u,w),vt));
	for(int i=0;i<n;i++)
	{
		if(std::isinf(1.0/w[i][i]) == false)
		{
			w[i][i]=1.0/w[i][i];
		}
	}
//getting X by multiplying V*W*UT
	x=multiply(v,w);
	x=multiply(x,ut);
	x=multiply(x,b);
	cout<<"X matrix:\n";
	print_matrix(x);
}	

int main()
{
	initialize();
	jacobi(aa, v);
	cout<<"\nJacobi done.\nNumber of sweeps: "<<sweep_count<<'\n';
	sweep_count=0;
	jacobi(aat, u);
	cout<<"Jacobi done.\nNumber of sweeps: "<<sweep_count<<'\n';
	svd();
	//print_matrix(aat);
	//print_matrix(v);
	cout<<"U matrix:\n";
	print_matrix(u);
	cout<<"VT matrix:\n";
	print_matrix(vt);
	return 0;
}
