#include "stdafx.h"
#include "varCal.h"
#include "csv.h"
#include <cmath>
#include "getPortfolio.h"
#include "matrixops.h"

void deltaGamma(vector<Future> futuresPortfolio, vector<vector<option>> optionsPortfolio,double** covMat,double* meanPtr,string* title, int nCol){

	map<string, map<string, double>> cov;
	for (int i1 = 0; i1 < nCol; ++i1) {
		for (int i2 = 0; i2 < nCol; ++i2) {
			cov[title[i1]][title[i2]] = covMat[i1][i2];
		}
	}

	//Delta Gamma VaR computation

	//Optimized search
	map<string, int> position;
	for (int i1 = 0; i1 < futuresPortfolio.size(); ++i1) {
		position[futuresPortfolio[i1].title] = i1;
	}
	for (int i2 = 0; i2 < optionsPortfolio.size(); ++i2){
		for (int i3 = 0; i3 < optionsPortfolio[i2].size(); ++i3){
			futuresPortfolio[position[optionsPortfolio[i2][i3].underlying]].position += optionsPortfolio[i2][i3].delta * optionsPortfolio[i2][i3].position;
			futuresPortfolio[position[optionsPortfolio[i2][i3].underlying]].gamma += optionsPortfolio[i2][i3].gamma * optionsPortfolio[i2][i3].position;
		}
	}

	map<string, double> mean;
	for (int i1 = 0; i1 < nCol; ++i1) {
		mean[title[i1]] = meanPtr[i1];
	}

	double **covMatrix12; // i,j element is Cov(ri,rj^2) where ri is the return of the i-th underlying
	covMatrix12 = new double *[futuresPortfolio.size()];
	for (int i1 = 0; i1 < futuresPortfolio.size(); ++i1) {
		covMatrix12[i1] = new double[futuresPortfolio.size()];
	}

	double **covMatrix22; //i,j element is Cov(ri^2,rj^2) where ri is the return of the i-th underlying
	covMatrix22 = new double *[futuresPortfolio.size()];
	for (int i1 = 0; i1 < futuresPortfolio.size(); ++i1){
		covMatrix22[i1] = new double[futuresPortfolio.size()];
	}

	string element1;
	string element2;

	cout << "Computing Covariance matrices" << endl;
	for (int i1 = 0; i1 < futuresPortfolio.size(); ++i1) {
		for (int i2 = 0; i2 < futuresPortfolio.size(); ++i2) {
			element1 = futuresPortfolio[i1].title;
			element2 = futuresPortfolio[i2].title;
			covMatrix12[i1][i2] = momentGenerate(element1, element2, mean, cov, 1, 2) - momentGenerate(element1, element2, mean, cov, 1, 0) * momentGenerate(element1, element2, mean, cov, 0, 2);
			covMatrix22[i1][i2] = momentGenerate(element1, element2, mean, cov, 2, 2) - momentGenerate(element1, element2, mean, cov, 2, 0) * momentGenerate(element1, element2, mean, cov, 0, 2);
		}
	}

	cout << "done" << endl;

	double **delta;
	delta = new double *[1];
	delta[0] = new double[futuresPortfolio.size()];

	double **gamma;
	gamma = new double *[1];

	gamma[0] = new double[futuresPortfolio.size()];
	for (int i1 = 0; i1 < futuresPortfolio.size(); ++i1) {
		delta[0][i1] = futuresPortfolio[i1].position * futuresPortfolio[i1].price;
		gamma[0][i1] = futuresPortfolio[i1].gamma* pow(futuresPortfolio[i1].price, 2);
	}

	double **deltaT;
	deltaT = new double *[futuresPortfolio.size()];
	double **gammaT;
	gammaT = new double *[futuresPortfolio.size()];
	for (int i1 = 0; i1 < futuresPortfolio.size(); ++i1) {
		deltaT[i1] = new double[1];
		gammaT[i1] = new double[1];
		deltaT[i1][0] = futuresPortfolio[i1].position * futuresPortfolio[i1].price;
		gammaT[i1][0] = futuresPortfolio[i1].gamma* pow(futuresPortfolio[i1].price, 2);
	}

	double deltaDelta = 0.0;
	double deltaGamma = 0.0;
	double gammaGamma = 0.0;

	double **result1;
	result1 = new double *[1];
	result1[0] = new double[futuresPortfolio.size()];

	double **result2;
	result2 = new double *[1];
	result2[0] = new double[1];

	//Computing delta*CovMat*deltaT
	matrixMultiply(delta, covMat, result1, 1, futuresPortfolio.size(), futuresPortfolio.size());
	matrixMultiply(result1, deltaT, result2, 1, futuresPortfolio.size(), 1);
	deltaDelta = result2[0][0];

	//Computing delta*Cov12*gammaT
	matrixMultiply(delta, covMatrix12, result1, 1, futuresPortfolio.size(), futuresPortfolio.size());
	matrixMultiply(result1, gammaT, result2, 1, futuresPortfolio.size(), 1);
	deltaGamma = result2[0][0];

	//Computing gamma*Cov22*gammaT
	matrixMultiply(gamma, covMatrix22, result1, 1, futuresPortfolio.size(), futuresPortfolio.size());
	matrixMultiply(result1, gammaT, result2, 1, futuresPortfolio.size(), 1);
	gammaGamma = result2[0][0];

	//Output VaR
	cout << "Value at risk is : " << 2.33 * sqrt(deltaDelta + deltaGamma + 0.25 * gammaGamma) << endl; //0.5 factor missing in deltaGamma term ?


	//Deallocating memory
	delete delta[0];
	delete[] delta;

	delete gamma[0];
	delete[] gamma;

	for (int i1 = 0; i1 < futuresPortfolio.size(); ++i1) {
		delete deltaT[i1];
		delete gammaT[i1];
	}
	delete[] deltaT;
	delete[] gammaT;

	for (int i1 = 0; i1 < futuresPortfolio.size(); ++i1) {
		delete covMatrix12[i1];
	}
	delete[] covMatrix12;
	for (int i1 = 0; i1 < futuresPortfolio.size(); ++i1) {
		delete covMatrix22[i1];
	}
	delete[] covMatrix22;

	delete result1[0];
	delete[] result1;

	delete result2[0];
	delete[] result2;
}