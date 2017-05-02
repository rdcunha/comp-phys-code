#include <vector>
#include <cmath>

using namespace std;

class Function
{
	public:
	const int nval=5;
	double func1[nval][2];
	double C[nval][nval];
	double D[nval][nval];
	double x=0;
	vector <double> func2;
	vector <double> derv1;
	vector <double> derv2;
	

	double function(double w);
	void makeCD(double y);
	double value(double z);
	double derive1(int j);
	double derive2(int k);
	
	Function();
	~Function();
};
