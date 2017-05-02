#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

//defining needed variables
vector <double> func;
vector <double> derv1;
vector <double> derv2;
double maxm = acos(-1);
double minm = asin(1);
double h=(maxm-minm)/100;

//function to calculate value of the function at a point
double function(double y)
{
	double val=y*y*y*cos(y);
	return val;
}

//implementing the 3 point formula for first derivatives
double derive1(int j)
{
	double der=(func[j+1]-func[j-1])/(2*h);
	//cout<<func[j+1]-func[j-1]<<'\n';
	return der;
}

//implementing the 3 point formula for second derivatives
double derive2(int k)
{
	double dder=(func[k+1]+func[k-1]-(2*func[k]))/(h*h);
	//cout<<func[k+1]-func[k-1]-(2*func[k])<<'\n';
	return dder;
}

//function to obtain values of the analytical first derivative at a point
double der_function(double z)
{
	double derf=(z*z)*(3*cos(z)-z*sin(z));
	return derf;
}

//function to obtain values of the analytical second derivative at a point
double dder_function(double w)
{
	double dderf=-(w*((-6+w*w)*cos(w)+6*w*sin(w)));
	return dderf;
}

int main()
{
	double x=minm;
	//cout<<x<<'\n';
	//cout<<function(x)<<'\n';
	func.push_back(function(x));
	//cout<<func[0];
//storing function values in a vector called func
	for(int i=0;i<100;++i)
	{
		x=x+h;
		func.push_back(function(x));
		//cout<<x<<" "<<func[i+1]<<'\n';
	}
//Calling the derive functions and storing the derivatives	
	for(int i=1;i<100;i++)
	{
		derv1.push_back(derive1(i));
		derv2.push_back(derive2(i));
	}
//Finding the boundary values using the 3 point boundary formula
	double end1=(-3*func[0]+4*func[1]-func[2])/(2*h);
	double end2=(-3*func[98]+4*func[99]-func[100])/(2*h);
	double ddend1=(2*func[0]-5*func[1]+4*func[2]-func[3])/(h*h);
	double ddend2=(2*func[100]-5*func[99]+4*func[98]-func[97])/(h*h);
//Adding the boundary values to the derivative vectors
	derv1.insert(derv1.begin(),end1);
	//cout<<derv1[0]<<'\n';
	derv1.push_back(end2);
	//cout<<end1<<" "<<end2<<'\n';
	derv2.insert(derv2.begin(),ddend1);
	//cout<<ddend1<<" "<<ddend2<<'\n';
	derv2.push_back(ddend2);

//Checking that the derivative values are being calculated correctly		
	/*cout<<"The derivatives at all 100 points:"<<'\n';
	cout<<"First      Second"<<'\n';*/
//printing first and second derivatives to files
	ofstream data1("first_derivative.dat");
	ofstream data2("second_derivative.dat");
	x=minm;
	for(int i=0;i<101;i++)
	{
		data1.precision(6);
		data2.precision(6);
		data1<<x<<" "<<derv1[i]<<'\n';
		data2<<x<<" "<<derv2[i]<<'\n';
		x+=h;
		//cout<<derv1[i]<<"  "<<derv2[i]<<'\n';
	}
	data1.close();
	data2.close();
	
	ofstream data11("first_analytical.dat");
	ofstream data22("second_analytical.dat");
	x=minm;
	for(int i=0;i<101;i++)
	{
		data11.precision(6);
		data22.precision(6);
		data11<<x<<" "<<derv1[i]-der_function(x)<<'\n';
		data22<<x<<" "<<derv2[i]-dder_function(x)<<'\n';
		x+=h;
		//cout<<derv1[i]<<"  "<<derv2[i]<<'\n';
	}
	data11.close();
	data22.close();
	cout<<"The numerical derivatives and differences from analytical values have printed to files *_derivatives.dat and *_analytical.dat."<<'\n';	
//Checking if the vector sizes are correct
	//cout<<derv1[100]<<'\n';
	//cout<<func.size()<<'\n';
	//cout<<derv1.size()<<'\n';
	//cout<<derv2.size()<<'\n';
	return 0;
}
