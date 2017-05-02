#include <cmath>
#include <vector>
#include <iostream>

using std::cin;
using std::cout;
using std::vector;

//defining necessary variables
vector < vector<double> > mat { {1,-1,-2,1,1},{3,2,-1,2,4},{2,3,1,-3,2},{10,-4,3,2,3} };
double alpha[4][4];
double beta[4][4];
vector < vector<double> > newmat;
vector <double> y;
vector <double> x;
vector <double> temp1;
vector < vector<double> > mult(4, vector<double>(4,0));

//initializing pivot cout, setting pivoting to on
int pivot_count=0;
bool pivot_flag=true;

//comparing all values greater than row j in a column
int compare_y(int j)
{
	int largest_y=j;
	for(int i=j;i<4;i++)
	{
		//cout<<"mat"<<i<<" "<<j<<" = "<<mat[i][j]<<'\n';
		if(abs(newmat[i][j])>abs(newmat[largest_y][j]))
		{
			largest_y=i;
		}
	}
	return largest_y;
}

//pivoting on k
void pivot(int k)
{
	int l=compare_y(k);
	for(int i=0;i<5;i++)
	{
		temp1.push_back(newmat[k][i]);
		newmat[k][i]=newmat[l][i];
		newmat[l][i]=temp1[i];
	}
	temp1.clear();
	//cout<<"Pivot no: "<<pivot_count+1<<'\n';
	/*for(int i=0;i<5;i++)
	{
		cout<<"\nmat "<<l<<" "<<i<<" "<<mat[l][i]<<'\n';
	}
	for(int i=0;i<5;i++)
	{
		cout<<"\nmat 0 "<<i<<" "<<mat[0][i]<<'\n';
	}*/
	++pivot_count;
}

//function to print the matrix
void print_matrix(vector < vector<double> > matrix)
{
	for(int j=0;j<4;j++)
	{
		for(auto i:matrix[j])
		{
			cout<<i<<" ";
		}
		cout<<'\n';
	}
}

void print_all()
{
	cout<<"Alpha:\n";
	for(int i=0;i<4;i++)
	{
		for(int j=0;j<4;j++)
		{
			cout<<alpha[i][j]<<" ";
		}
		cout<<'\n';
	}
	cout<<"Beta:\n";
	for(int i=0;i<4;i++)
	{
		for(int j=0;j<4;j++)
		{
			cout<<beta[i][j]<<" ";
		}
		cout<<'\n';
	}
	for(int j=0;j<4;j++)
	{
		for(auto i:newmat[j])
		{
			cout<<i<<" ";
		}
		cout<<'\n';
	}
	cout<<"Matrix multiplication of alpha and beta:\n";
	for(int j=0;j<4;j++)
	{
		for(auto i:mult[j])
		{
			cout<<i<<" ";
		}
		cout<<'\n';
	}
}

//initializing alpha and beta matrices
void initialize()
{
	for(int i=0;i<4;i++)
	{
		alpha[i][i]=1;
	}
	for(int i=0;i<4;i++)
	{
		for(int j=i;j<4;j++)
		{
			if(j==i)
				alpha[i][j]=1;
			else
				alpha[i][j]=0;
		}
		for(int j=0;j<i+1;j++)
		{
			beta[i][j]=0;
		}
	newmat=mat;
	}
}

//using the formula to obtain non-zero elements of alpha and beta
void makeab()
{
	for(int j=0;j<4;j++)
	{
		if(pivot_flag==true && j<3)
		{
			pivot(j);
		}
		//cout<<"After pivoting:\n";
		//print_matrix(newmat);
		//print_all();
		for(int i=0;i<j+1;i++)
		{
			double sum2=0;
			for(int k=0;k<i;k++)
			{
				sum2+=newmat[i][k]*newmat[k][j];
			}
			newmat[i][j]=(newmat[i][j]-sum2);
		}
		for(int i=j+1;i<4;i++)
		{
			double sum1=0;
			for(int k=0;k<j;k++)
			{
				sum1+=newmat[i][k]*newmat[k][j];
			}
			newmat[i][j]=(newmat[i][j]-sum1)/newmat[j][j];
		}
		//cout<<"After updating:\n";
		//print_matrix(newmat);
		//print_all();
	}
}

//placing the elements of alpha and beta in the alpha and beta matrices
void fill()
{
	for(int i=0;i<4;i++)
	{
		for(int j=0;j<i;j++)
		{
			//alpha[i][j]=pow(-1,pivot_count)*newmat[i][j];
			alpha[i][j]=newmat[i][j];
		}
		for(int j=i;j<4;j++)
		{
			beta[i][j]=newmat[i][j];
		}
	}
}

//solving the 2 sets of linear equations with triangular matrices
void solve()
{
	y.push_back(mat[3][4]/alpha[0][0]);
	y.push_back((mat[2][4]-y[0]*alpha[1][0])/alpha[1][1]);
	y.push_back((mat[0][4]-y[1]*alpha[2][1]-y[0]*alpha[2][0])/alpha[2][2]);
	y.push_back((mat[1][4]-y[2]*alpha[3][2]-y[1]*alpha[3][1]-y[0]*alpha[3][0])/alpha[3][3]);
	x.push_back(y[3]/beta[3][3]);
	x.push_back((y[2]-beta[2][3]*x[0])/beta[2][2]);
	x.push_back((y[1]-beta[1][3]*x[0]-beta[1][2]*x[1])/beta[1][1]);
	x.push_back((y[0]-beta[0][3]*x[0]-beta[0][2]*x[1]-beta[0][1]*x[2])/beta[0][0]);
}

void matmultiply()
{
	for(int i=0;i<4;i++)
	{
		for(int j=0;j<4;j++)
		{
			mult[i][j]=0;
			for(int k=0;k<4;k++)
			{
				mult[i][j]+=alpha[i][k]*beta[k][j];
			}
		}
	}
}

double determinant()
{
	double det=1;
	for(int i=0;i<4;i++)
	{
		det*=beta[i][i];
	}
	return det;
}

int main()
{
	initialize();
	//print_all();
	//cout<<'\n';
	//print_matrix(mat);
	makeab();
	fill();
	solve();
	matmultiply();
	print_all();
	cout<<"Y matrix:\n";
	for(int i=0;i<4;i++)
	{
		cout<<y[i]<<" ";
	}
	cout<<'\n';
	cout<<"u = "<<x[0]<<" z = "<<x[1]<<" y = "<<x[2]<<" x = "<<x[3];
	cout<<'\n';
	cout<<"The determinant of A: "<<pow(-1,pivot_count)*determinant()<<'\n';
	return 0;
}
