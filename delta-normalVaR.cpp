#include "stdafx.h"
#include "varCal.h"
#include "csv.h"
#include "matrixops.h"
#include "normalDist.h"
#include "EMtry.h"
#include "varMethods.h"
#include "getPortfolio.h"
#include <map>
#include <iostream>


void deltaNormal(vector<Future> futuresPortfolio, vector<vector<option>> optionsPortfolio,double **covMat){
	//Computing delta adjusted position
	//Optimized O(M+V) search where M: number of underlying V: number of options 
	map<string, int> position;
	for (int i1 = 0; i1 < futuresPortfolio.size(); ++i1) {
		position[futuresPortfolio[i1].title] = i1;
	}
	for (int i2 = 0; i2 < optionsPortfolio.size(); ++i2) {
		for (int i3 = 0; i3 < optionsPortfolio[i2].size(); ++i3) {
				futuresPortfolio[position[optionsPortfolio[i2][i3].underlying]].position += optionsPortfolio[i2][i3].delta * optionsPortfolio[i2][i3].position;
		}
	}

	double **delta;
	double **deltaT; //Transpose of delta
	delta = new double*[1]; //1xn vector
	deltaT = new double*[futuresPortfolio.size()]; //nx1 vector
	delta[0] = new double[futuresPortfolio.size()];
	for (int i1 = 0; i1 < futuresPortfolio.size(); ++i1) {
		delta[0][i1] = futuresPortfolio[i1].position * futuresPortfolio[i1].price;
		deltaT[i1] = new double[1];
		deltaT[i1][0] = futuresPortfolio[i1].position * futuresPortfolio[i1].price;
	}
	
	//Computing delta*CovMat*deltaT
	double **result1;
	result1 = new double *[1];
	result1[0] = new double[futuresPortfolio.size()];
	double **result2;
	result2 = new double *[1];
	result2[0] = new double[1];
	matrixMultiply(delta, covMat, result1, 1, futuresPortfolio.size(), futuresPortfolio.size());
	matrixMultiply(result1, deltaT, result2, 1, futuresPortfolio.size(), 1);
	cout << "Value at Risk: " << 2.33 * sqrt(result2[0][0]) << endl;

	//Deallocating memory
	delete delta[0];
	delete[] delta;

	for (int i1 = 0; i1 < futuresPortfolio.size(); ++i1) {
		delete deltaT[i1];
	}
	delete[] deltaT;

	delete result1[0];
	delete[] result1;

	delete result2[0];
	delete[] result2;
}




