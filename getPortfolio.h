#define GET_PORTFOLIO
#ifndef M_PI
#define M_PI 3.14159265358979323846


#include "csv.h"
#include "varCal.h"
#include "EMtry.h"
#include <iostream>
#include <string>
#include <vector>
#include <map>

using namespace std;

struct Future {
public:
	string title;
	double position;
	double price;
	double gamma;
	double priceChange;
	double alpha;
	double sigma;
	double v;
	
};

struct option{
public:
	string underlying;
	double underlyingPrice;

	double price;
	double position;
	double maturity;
	double rd;
	double rf;
	double impliedVol;
	string optionType;
	double delta;
	double d1;
	double gamma;
	double strike;

	double MCunderlyingPrice;
	double impliedVol_LP;
	double priceEvaluatedByBSformula;
	double priceEvaluatedByVG;

	double imp_alpha;
	double imp_sigma;
	double imp_v;
};

vector<Future> getFuturesPortfolio();
vector<vector<option>> getOptionsPortfolio(map<string, map<double, double> >& impliedVol_grid);
CovCal getHistoricalData();
double** getCovMat(int& col);
double* getMean();
string* getTitle();

#endif