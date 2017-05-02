#include <iostream>
#include <cmath>
#include <vector>

using std::cin;
using std::cout;
using std::vector;

//defining necessary variables
vector < vector<double> > mat { {1,-1,-2,1,1},{3,2,-1,2,4},{2,3,1,-3,2},{10,-4,3,2,3} };
vector <double> temp;
vector <double> var;

//Deciding whether to use pivoting or not
bool pivot_flag=true;

//finding the largest element in the row for normalization		
int compare_x(int j)
{
	int largest_x=0;
	for(int i=0;i<4;i++)
	{
		//cout<<"mat"<<i<<" "<<j<<" = "<<mat[i][j]<<'\n';
		if(abs(mat[j][i])>abs(mat[j][largest_x]))
		{
			largest_x=i;
		}
	}
	return largest_x;
}

//normalizing using the largest element in the row
void normalize()
{
	double temp2=0;
	for(int i=0;i<4;i++)
	{
		temp2=mat[i][compare_x(i)];
		for(int j=0;j<5;j++)
		{
			mat[i][j]/=temp2;
		}
	}
}

//finding the largest element in the column for pivoting
int compare_y(int j)
{
	int largest_y=j;
	for(int i=j;i<4;i++)
	{
		//cout<<"mat"<<i<<" "<<j<<" = "<<mat[i][j]<<'\n';
		if(abs(mat[i][j])>abs(mat[largest_y][j]))
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
		temp.push_back(mat[k][i]);
		mat[k][i]=mat[l][i];
		mat[l][i]=temp[i];
	}
	temp.clear();
	/*for(int i=0;i<5;i++)
	{
		cout<<"\nmat "<<l<<" "<<i<<" "<<mat[l][i]<<'\n';
	}
	for(int i=0;i<5;i++)
	{
		cout<<"\nmat 0 "<<i<<" "<<mat[0][i]<<'\n';
	}*/
}

//making elements below the pivot element zero
void clear(int m)
{
	double temp1;
	for(int i=m+1;i<4;i++)
	{
		temp1=mat[i][m];
		for(int j=0;j<5;j++)
		{
			mat[i][j]=mat[i][j]-(mat[m][j]*(temp1/mat[m][m]));
			//cout<<i<<j<<"factor :"<<mat[i][m]<<mat[m][m]<<" "<<mat[i][m]/mat[m][m]<<'\n';
		}
	}
}

//function to print the matrix
void print_matrix()
{
	for(int i=0;i<4;i++)
	{
		for(int j=0;j<5;j++)
		{
			cout<<mat[i][j]<<" ";
		}
		cout<<'\n';
	}
}

//function to reduce the matrix to row echelon form
void reduce()
{
	for(int i=0;i<3;i++)
	{
		if(pivot_flag==true)
			pivot(i);
	//cout<<"Pivot number "<<i<<'\n';
	//print_matrix();
	clear(i);
	//cout<<"Clear number "<<i<<'\n';
	//print_matrix();
	}
}

void solve()
{
	for(int i=3;i>=0;--i)
	{
	double diff=mat[i][4];
	//cout<<diff<<'\n';
	int count=0;
	if(i==3)
	{
		var.push_back(diff/mat[i][i]);
		cout<<var[0]<<'\n';
	}
	else
	{
		for(int j=4;j>i+1;--j)
		{
			diff=diff-(var[count]*mat[i][j-1]);
			//cout<<diff<<'\n';
			count++;
		}
		var.push_back(diff/mat[i][i]);
	}
	}
	//z=(mat[2][4]-(u*mat[2][3]))/mat[2][2];
	//y=(mat[1][4]-(z*mat[1][2])-(u*mat[1][3]))/mat[1][1];
	//x=(mat[0][4]-(y*mat[0][1])-(z*mat[0][2])-(u*mat[0][3]))/mat[0][0];
}

int main()
{
	cout<<"Initial matrix:\n";
	print_matrix();
	normalize();
	cout<<"Normalized matrix:\n";
	print_matrix();
	reduce();
	cout<<"Matrix in row echelon form:\n";
	print_matrix();
	solve();
	cout<<"\nThe values :x = "<<var[3]<<" y = "<<var[2]<<" z = "<<var[1]<<" u = "<<var[0]<<'\n';
	return 0;
}
