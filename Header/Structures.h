// Data structures

#ifndef Structures_H
#define Structures_H

//#include <vector>
using namespace std;

struct RMSEinputs{
	double S;
	vector<double> K;
	vector<double> MPrice;
	double r;
	double T;
};

struct NMsettings{
	int N;
	double MaxIters;
	double tolerance;
	RMSEinputs RMSEinp;
};

#endif // Structures_H 
