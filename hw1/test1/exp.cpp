#include <cmath>
#include <iostream>

using namespace std;

//the function's values must be specified here
const int nval=5;
double func[nval][2]={{0, 1.0},{0.5,1.64872},{1,2.71828},{2,7.38906},{4,54.59815}};

//defining other needed variables
double C[nval][nval];
double D[nval][nval];
double x=0;

//function to use the recursive formulas to get C and D arrays
void makeCD(double y)
{
	for(int m=0;m<(nval-1);++m)
	{
		for(int i=0;i<(nval-m);++i)
		{
			C[m+1][i]=((y-func[i][0])*(C[m][i+1]-D[m][i] ))/(func[i+m+1][0]-func[i][0]);
			D[m+1][i]=((y-func[i+m+1][0])*(D[m][i]-C[m][i+1] ))/(func[i][0]-func[i+m+1][0]);
		}
	}
}

//function to return the value needed
double value(double z)
{
	makeCD(z);
//checking if C and D values calculated are correct
	/*for(int m=0;m<5;m++)
	{
		for(int i=0;i<5-m;++i)
		{
		cout<<"C"<<m<<i<<" : "<<C[m][i]<<" D"<<m<<i<<": "<<D[m][i]<<'\n';
		}
	}*/
	double val=func[0][1]+C[1][0]+C[2][0]+C[3][0]+C[4][0];
	return val;
}

//function for linear interpolation
double linear(double w)
{
	double linval=0;
	if(w==0.16)
	{
		linval=(w-func[0][0])*((func[1][1]-func[0][1])/(func[1][0]-func[0][0]));
		linval=func[0][1]+linval;
	}
	if(w==0.46)
	{
		linval=(w-func[2][0])*((func[3][1]-func[2][1])/(func[3][0]-func[2][0]));
		linval=func[2][1]+linval;
	}
	return linval;
}

int main()
{
//initializing C0i's and D0i's to be fi's
	for(int i=0;i<nval;++i)
	{
		C[0][i]=func[i][1];
		D[0][i]=func[i][1];
	}
//checking the initialization
	/*for(int i=0;i<5;++i)
	{
		cout<<"C0"<<i<<" : "<<C[0][i]<<" D0"<<i<<": "<<D[0][i]<<'\n';
	}*/
//Calling the value() function and outputting the result
	cout<<"The function's value at 0.16 is "<<value(0.16)<<'\n';
	cout<<"The function's value at 0.46 is "<<value(0.46)<<'\n';
	/*cout<<"The function's value by linear interpolation  at 0.46 is "<<linear(0.46)<<'\n';
	cout<<"The function's value by linear interpolation at 0.16 is "<<linear(0.16)<<'\n';*/
	cout<<"The numerical uncertainty is: "<<0.5*abs((C[4][0])+abs(D[4][0]))<<'\n';
	return 0;
}
