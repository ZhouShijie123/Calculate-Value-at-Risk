#include "stdafx.h"
#include <iostream>
#include "matrixops.h"
#include "normalDist.h"
#include "varMethods.h"
#include "getPortfolio.h"
#include <fstream>

using namespace std;

int main(){
	
	//Getting Data
	cout << "Processing Futures Portfolio ..." << endl;
	vector<Future> futuresPortfolio = getFuturesPortfolio();
	cout << "done" << endl;
	cout << "Processing Options Portfolio ..." << endl;
	//impliedVol added, for further linear interpolation 
	map<string, map<double, double> > impliedVol_grid;
	vector<vector<option>> optionsPortfolio = getOptionsPortfolio(impliedVol_grid);

	cout << "Getting Covariance Matrix ..." << endl;
	//Check if cov and mean are already computed and stored, if not fix/fill data with EM algorithm and compute them
	ifstream meanfile("mean.csv");
	ifstream covfile("cov.csv");
	ifstream contractfile("title.csv");
	if (!meanfile | !covfile ){
		cout << "Computing Covariance Matrix ..." << endl;
		CovCal historicalData = getHistoricalData();
	}


	covfile.close();
	meanfile.close();
	contractfile.close();
	int nCol=0; //number of underlying assets(futures)
	double ** covMat = getCovMat(nCol);
	double* meanPtr = getMean();
	string* title = getTitle();
	VG_VaR(futuresPortfolio, optionsPortfolio, covMat, title, nCol, meanPtr);


	//Delta Normal Value at Risk
	cout << "================================" << endl;
	cout << "Delta Normal Value at Risk: " << endl;
	deltaNormal(futuresPortfolio,optionsPortfolio,covMat);
	cout << "================================" << endl;
	
	//Delta Gamma Value at Risk
	cout << "================================" << endl;
	cout << "Delta Gamma Value at Risk: " << endl;
	//deltaGamma(futuresPortfolio,optionsPortfolio,covMat,meanPtr,title,col);
	cout << "================================" << endl;

	
	//Historical Value at Risk
	cout << "================================" << endl;
	cout << "Historical Value at Risk: " << endl;
	historical(futuresPortfolio, optionsPortfolio, impliedVol_grid);
	cout << "================================" << endl;
	
	//MonteCarlo Value at Risk
	cout << "================================" << endl;
	cout << "MonteCarlo Value at Risk: " << endl;
	MonteCarloVaR(futuresPortfolio, optionsPortfolio, impliedVol_grid, covMat, title, nCol);
	cout << "================================" << endl;


	cin.get();
	return 0;
}

