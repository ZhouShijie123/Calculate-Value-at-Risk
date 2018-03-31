#ifndef VAR_METHODS
#define VAR_METHODS
#include<vector>
#include<string>
#include"getPortfolio.h"
using namespace std;
void deltaNormal(vector<Future> futuresPortfolio, vector<vector<option>> optionsPortfolio, double** covMat);
void deltaGamma(vector<Future> futuresPortfolio, vector<vector<option>> optionsPortfolio, double** covMat, double* meanPtr, string* title, int nCol);
void historical(vector<Future> futuresPortfolio, vector<vector<option>> optionsPortfolio, map<string, map<double, double> > impliedVol_grid);
void MonteCarloVaR(vector<Future> futuresPortfolio, vector<vector<option>> optionsPortfolio, map<string, map<double, double> > impliedVol_grid,
	double** covMat, string* title, int nCol);
void VG_VaR(vector<Future> futuresPortfolio, vector<vector<option>> optionsPortfolio, double ** CovMat, string* title, int nCol, double *mean);
#endif