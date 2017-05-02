#include <cmath>
#include <iostream>
#include <vector>

using std::cin;
using std::cout;
using std::vector;

//defining needed variables and arrays
vector <double> func;
vector <double> pt;
double maxm=3*acos(0);
double minm=acos(0);
double guess=maxm;
double h=(maxm-minm)/100;
double guess1=maxm-h;
double a=5.0;
double epsilon=0.001;
int i;

//function to find the next point close to the root using the secant method
double find_pt(int n)
{
	double val=(pt[n-1]-pt[n-2])/(func[n-1]-func[n-2]);
	val=func[n-1]*val;
	val=pt[n-1]-val;
	return val;
}

//obtaining the value of the function at a point
double function(double y)
{
	double value=tan(y)-(a/y);
	return value;
}

int main()
{	
	//cout<<maxm<<" "<<minm<<" "<<h<<'\n';
//stepping through the range to find a range where the function crosses the axis
	while(function(guess)*function(guess1) > 0)
	{
		guess1=guess1-h;
		guess=guess-h;
		//cout<<guess<<" "<<guess1<<'\n';
	}
//storing the optimized guess values
	pt.push_back(guess);
	func.push_back(function(guess));
	pt.push_back(guess1);
	func.push_back(function(guess1));
//checking the guess values
	//cout<<pt[0]<<" "<<pt[1]<<'\n';
	//cout<<func[0]<<" "<<func[1]<<'\n';
//using the secant method to find closer values to the root and storing them
	for(i=2;i<101;i++)
	{
		pt.push_back(find_pt(i));
		func.push_back(function(pt[i]));
//convergence criterion
		if(pt[i]-pt[i-1]<epsilon)
			break;
	}
	cout<<"The value of the root is: "<<pt[i]<<'\n';
	return 0;
}
