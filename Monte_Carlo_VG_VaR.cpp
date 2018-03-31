#include "stdafx.h"
#include <iostream>
#include "matrixops.h"
#include "normalDist.h"
#include "varMethods.h"
#include "getPortfolio.h"
#include <fstream>
#include <string>
#include <math.h>
#include <ctime>
#include <random>
#include "BlackScholes.h"
#include "NelderMead.h"
#include "VGFiles.h"
#include "RMSE.h"
using namespace std;
double normcdf(double z) {
	double pi = 3.141592653589793;
	double b1 = -0.0004406;
	double b2 = 0.0418198;
	double b3 = 0.9;
	if (z < -8)
		return 0.0;
	else if (z > 8)
		return 1.0;
	else
		return 1.0 / (1.0 + exp(-sqrt(pi)*(b1*pow(z, 5.0) + b2*pow(z, 3.0) + b3*z)));
}

double VGCall(double alpha, double sigma, double v, double S, double K, double r, double T, string type) {
	double a = pow(alpha + sigma, 2.0);
	double num = 1.0 - v*a / 2.0;
	double den = 1.0 - v*alpha*alpha / 2.0;
	double d1 = log(S / K) / sigma / T + ((r + 1.0 / v*log(num / den)) / sigma + alpha + sigma)*sqrt(T);
	double d2 = d1 - sigma*sqrt(T);
	double temp = S*exp(a*T / 2.0) * pow((1.0 - v*a / 2.0), T / v)*normcdf(d1)
		- K*exp(-r*T + alpha*alpha*T / 2.0) * pow((1.0 - v*alpha*alpha / 2.0), T / v)*normcdf(d2);
	if (type == "call")
		return temp;
	else if (type == "put")
		return temp + K * exp(-r * T) - S;
	else {
		cout << "Should be call or put";
		return 0;
	}
}
double RMSE(vector<double> B, RMSEinputs rmsein) {
	CVGFiles VG;
	int N = rmsein.K.size();
	vector<double> VGPrice(N, 0.0);
	double Sum = 0.0;
	// Penalty for negative parameter values
	if ((B[1]<0.0) || (B[2]<0.0))
		return 1.0e100;
	else {
		// Create the loss function
		for (int i = 0; i <= N - 1; i++) {
			VGPrice[i] = VG.VGCall(B[0], B[1], B[2], rmsein.S, rmsein.K[i], rmsein.r, rmsein.T);
			// Sum +=  pow( (rmsein.MPrice[i] - VGPrice[i])/rmsein.MPrice[i], 2.0);
			// Sum +=  pow((log(rmsein.MPrice[i]) - log(VGPrice[i])), 2.0);
			Sum += pow(rmsein.MPrice[i] - VGPrice[i], 2.0);
		}
		return Sum;
	}
}

double generator(double theta, double sigma, double v,int n) {
	double delta_t = 1.0; 
	
	default_random_engine generator_gamma;
	generator_gamma.seed(n);
	gamma_distribution<double> distribution_gamma(delta_t / v, v);
	double delta_G = distribution_gamma(generator_gamma);

	default_random_engine generator_normal;
	generator_normal.seed(n);
	normal_distribution<double> distribution_normal(0.0, 1.0);
	double z = distribution_normal(generator_normal);

	double res = delta_G * theta + sigma * sqrt(delta_G) * z;
	return res;
}

void VG_VaR(vector<Future> futuresPortfolio,vector<vector<option>> optionsPortfolio, double ** CovMat, string* title, int nCol, double *mean) {
	
	int numOfPath = 1000;
	double rf = 0.0;
	
	for (int i = 0; i < optionsPortfolio[0].size(); i += 2) {
		string m = optionsPortfolio[0][i].underlying;
		m = m + ".csv";
		ofstream output(m);
		cout << optionsPortfolio[0][i].underlying << endl;
		output << optionsPortfolio[0][i].underlying << endl;
		double RMSE(vector<double>, RMSEinputs);
		CNelderMead NM;
		CVGFiles VG;
		vector<double> strike;
		vector<double> call_price;
		double T = 0;
		double spot = 0;
		for (int j = 0; j < optionsPortfolio.size(); j++) {
			
			if (optionsPortfolio[j][i].optionType == "call") {
				call_price.push_back(optionsPortfolio[j][i].price);
				strike.push_back(optionsPortfolio[j][i].strike);
				T = optionsPortfolio[j][i].maturity;
				spot = optionsPortfolio[j][i].underlyingPrice;
			}
		}
		RMSEinputs rmsein;
		rmsein.K = strike;
		rmsein.MPrice = call_price;
		rmsein.S = spot;				// Spot Price
		rmsein.T = T;			// Maturity in Years
		rmsein.r = 0.0179;				// Interest Rate
		double q = 0.0;					// Dividend Yield

		NMsettings nmset;
		nmset.N = 3;						// Number of Variance Gamma parameters
		nmset.MaxIters = 2000;				// Maximum number of iterations
		nmset.tolerance = 1.0e-12;			// Tolerance on best and worst function values
		nmset.RMSEinp = rmsein;				// Settings for the RMSE function

		double alpha = 0.1;
		double sigma = 0.1;
		double v = 0.1;

		// VG parameter starting values (vertices) in vector form
		int N = nmset.N;
		double a = -0.1;
		double b = 0.1;
		vector<vector<double> > start(N, vector<double>(N + 1));
		start[0][0] = alpha + VG.U(a, b);	start[0][1] = alpha + VG.U(a, b);	start[0][2] = alpha + VG.U(a, b);	start[0][3] = alpha + VG.U(a, b);	    // alpha
		start[1][0] = sigma + VG.U(a, b);	start[1][1] = sigma + VG.U(a, b);	start[1][2] = sigma + VG.U(a, b);	start[1][3] = sigma + VG.U(a, b);		// gamma
		start[2][0] = v + VG.U(a, b);		start[2][1] = v + VG.U(a, b);		start[2][2] = v + VG.U(a, b);		start[2][3] = v + VG.U(a, b);		// v

		CRMSE crm;
		vector<double> B = NM.NelderMeadFunction(RMSE, nmset, start);		
		output << "Strike Price,BS Price,VG Price,Market Price" << endl;
		for (int j = 0; j < optionsPortfolio.size(); j++) {
			if (optionsPortfolio[j][i].optionType == "call") {
				optionsPortfolio[j][i].imp_alpha = B[0];
				optionsPortfolio[j][i].imp_sigma = B[1];
				optionsPortfolio[j][i].imp_v = B[2];

				optionsPortfolio[j][i + 1].imp_alpha = B[0];
				optionsPortfolio[j][i + 1].imp_sigma = B[1];
				optionsPortfolio[j][i + 1].imp_v = B[2];
				double temp = VGCall(optionsPortfolio[j][i].imp_alpha, optionsPortfolio[j][i].imp_sigma, optionsPortfolio[j][i].imp_v, optionsPortfolio[j][i].underlyingPrice, optionsPortfolio[j][i].strike, rf, optionsPortfolio[j][i].maturity, optionsPortfolio[j][i].optionType);
				//cout << temp << endl;
				double temp2 = optionPriceBS(optionsPortfolio[j][i].underlyingPrice, optionsPortfolio[j][i].strike, optionsPortfolio[j][i].rd, optionsPortfolio[j][i].rf, optionsPortfolio[j][i].maturity, optionsPortfolio[j][i].impliedVol, optionsPortfolio[j][i].optionType);
				output << optionsPortfolio[j][i].strike << "," << temp2 << "," << temp << "," << optionsPortfolio[j][i].price << endl;
			}

		}
		output.close();
	
		/*
	for (int i2 = 0; i2 < optionsPortfolio.size(); i2++) {
		for (int i3 = 0; i3 < optionsPortfolio[i2].size(); i3++) {
			cout << optionsPortfolio[i2][i3].underlying;
		}
		cout << endl;
	}
	*/
	}
	
	/*
	for (int i2 = 0; i2 < optionsPortfolio.size(); ++i2) {
		for (int i3 = 0; i3 < optionsPortfolio[i2].size(); ++i3) {
			double temp = VGCall(optionsPortfolio[i2][i3].imp_alpha, optionsPortfolio[i2][i3].imp_sigma, optionsPortfolio[i2][i3].imp_v, optionsPortfolio[i2][i3].underlyingPrice, optionsPortfolio[i2][i3].strike, rf, optionsPortfolio[i2][i3].maturity,optionsPortfolio[i2][i3].optionType);
			cout << "Market Price : " << optionsPortfolio[i2][i3].price << "   Cal_Price: " << temp << endl;
			cout << optionsPortfolio[i2][i3].underlyingPrice << "  " << optionsPortfolio[i2][i3].strike << "   " << optionsPortfolio[i2][i3].optionType << endl;
			cout << "***************************************************" << endl;
		}
	}
	*/



	ifstream inFile;
	inFile.open("data_est.txt");
	double m1, m2;
	
	for (int j = 0; j < 40; j++) {
		inFile >> m1;
		futuresPortfolio[j].sigma = m1;
	}
	for (int j = 0; j < 40; j++) {
		inFile >> m1;
		futuresPortfolio[j].alpha = 0;// m1;
	}
	for (int j = 0; j < 40; j++) {
		inFile >> m1;
		futuresPortfolio[j].v = m1;
	}
	double ** simuDataMatrix1 = new double*[numOfPath];
	for (int i = 0; i < numOfPath; i++) {
		simuDataMatrix1[i] = new double[futuresPortfolio.size()];
	}

	for (int i = 0; i < 1000; i++) {
		for (int j = 0; j < futuresPortfolio.size(); j++) {
			simuDataMatrix1[i][j] = generator(futuresPortfolio[j].alpha, futuresPortfolio[j].sigma, futuresPortfolio[j].v,i) - 0.1;
		}
	}

	double** squareRootOfCov;
	squareRootOfCov = new double*[nCol];
	for (int i1 = 0; i1 < nCol; ++i1) {
		squareRootOfCov[i1] = new double[nCol];
	}

	getSquareRootofSymMat(CovMat, squareRootOfCov, nCol);
	double** transposeOfSquareRoot;
	transposeOfSquareRoot = new double*[nCol];
	for (int i1 = 0; i1 < nCol; ++i1) {
		transposeOfSquareRoot[i1] = new double[nCol];
	}
	for (int i1 = 0; i1 < nCol; ++i1) {
		for (int i2 = 0; i2 < nCol; ++i2) {
			transposeOfSquareRoot[i1][i2] = squareRootOfCov[i2][i1];
		}
	}
	double **simuDataMatrix2;
	simuDataMatrix2 = new double *[numOfPath];
	for (int i1 = 0; i1 < numOfPath; ++i1) {
		simuDataMatrix2[i1] = new double[nCol];
	}

	matrixMultiply(simuDataMatrix1, transposeOfSquareRoot, simuDataMatrix2, numOfPath, nCol, nCol);
	for (int i1 = 0; i1 < numOfPath; ++i1) {
		for (int i2 = 0; i2 < nCol; ++i2) {
			simuDataMatrix2[i1][i2] += mean[i2];
		}
	}
	ofstream outputFile;
	outputFile.open("SimulationDataVG_theta0.csv");
	outputFile << "date,";
	for (int i1 = 0; i1 < nCol; ++i1) {
		outputFile << title[i1] << ",";
	}
	outputFile << 0 << endl;
	for (int i1 = 0; i1 < numOfPath; ++i1) {
		outputFile << "time,";
		for (int i2 = 0; i2 < nCol; ++i2) {
			outputFile << simuDataMatrix2[i1][i2] << ",";
		}
		outputFile << 0 << endl;
	}
	outputFile.close();

	double * diff = new double[numOfPath];
	double * futurediff = new double[numOfPath];
	/*
	
	for (int i = 0; i < numOfPath; i++) {
		double total_change = 0;
		for (int j = 0; j < futuresPortfolio.size(); j++) {
			double simulated_price = futuresPortfolio[j].price * (1 + simuDataMatrix2[i][j]);
			total_change += futuresPortfolio[j].price * simuDataMatrix2[i][j] * futuresPortfolio[j].position;
			for (int k = 0; k < optionsPortfolio.size(); k++) {
				for (int q = 0; q < optionsPortfolio[k].size(); q++) {
					if (futuresPortfolio[j].title == optionsPortfolio[k][q].underlying && optionsPortfolio[k][q].position != 0) {
						optionsPortfolio[k][q].MCunderlyingPrice = simulated_price;
						double temp = VGCall(futuresPortfolio[j].implied_alpha, futuresPortfolio[j].implied_sigma, futuresPortfolio[j].implied_v, simulated_price,optionsPortfolio[k][q].strike,rf,optionsPortfolio[k][q].maturity - 1.0/252.0);
						if (optionsPortfolio[k][q].optionType == "call")
							optionsPortfolio[k][q].priceEvaluatedByVG = temp;
						else
							optionsPortfolio[k][q].priceEvaluatedByVG = temp - simulated_price + optionsPortfolio[k][q].strike * exp(-rf * optionsPortfolio[k][q].maturity);
						total_change += (optionsPortfolio[k][q].priceEvaluatedByVG - optionsPortfolio[k][q].price) * optionsPortfolio[k][q].position;
					}
				}
			}
		}
		diff[i] = total_change;
		
		//cout << i << " : " << total_change << endl;
	}
	*/
	double * simulatedPortPriceDif = new double[numOfPath];
	for (int i = 0; i < numOfPath; i++) {
		simulatedPortPriceDif[i] = 0.0;
		futurediff[i] = 0.0;
	}
	CovCal simulationData("SimulationDataVG.csv");
	for (int i1 = 0; i1 < numOfPath; ++i1) {
		for (int i2 = 0; i2 < futuresPortfolio.size(); ++i2) {
			for (int i3 = 0; i3 < simulationData.nCol; ++i3) {
				if (simulationData.title[i3] == futuresPortfolio[i2].title) {					
					futuresPortfolio[i2].priceChange = futuresPortfolio[i2].price * simulationData.OriginalData[i1][i3];
				}
			}
			simulatedPortPriceDif[i1] += futuresPortfolio[i2].priceChange * futuresPortfolio[i2].position;
			futurediff[i1] += futuresPortfolio[i2].priceChange * futuresPortfolio[i2].position;
		}
		double total = 0;
		for (int i2 = 0; i2 < optionsPortfolio.size(); ++i2) {
			for (int i3 = 0; i3 < optionsPortfolio[i2].size(); ++i3) {
				for (int i4 = 0; i4 < simulationData.nCol; ++i4) {
					if ((optionsPortfolio[i2][i3].underlying == simulationData.title[i4]) && (optionsPortfolio[i2][i3].position != 0.0)) {
						optionsPortfolio[i2][i3].MCunderlyingPrice = optionsPortfolio[i2][i3].underlyingPrice * (1.0 + simulationData.OriginalData[i1][i4]);
						
						//optionsPortfolio[i2][i3].priceEvaluatedByVG = optionPriceBS(optionsPortfolio[i2][i3].MCunderlyingPrice, optionsPortfolio[i2][i3].strike, optionsPortfolio[i2][i3].rd, optionsPortfolio[i2][i3].rf, optionsPortfolio[i2][i3].maturity - (1.0 / 252.0), optionsPortfolio[i2][i3].impliedVol_LP, optionsPortfolio[i2][i3].optionType);
						optionsPortfolio[i2][i3].priceEvaluatedByVG = VGCall(optionsPortfolio[i2][i3].imp_alpha, optionsPortfolio[i2][i3].imp_sigma, optionsPortfolio[i2][i3].imp_v, optionsPortfolio[i2][i3].MCunderlyingPrice, optionsPortfolio[i2][i3].strike, rf, optionsPortfolio[i2][i3].maturity - (1.0 / 252.0), optionsPortfolio[i2][i3].optionType);
						if (isnan(optionsPortfolio[i2][i3].priceEvaluatedByVG - optionsPortfolio[i2][i3].price)) {
							simulatedPortPriceDif[i1] += 0;
						}
						else {
							simulatedPortPriceDif[i1] += optionsPortfolio[i2][i3].position * (optionsPortfolio[i2][i3].priceEvaluatedByVG - optionsPortfolio[i2][i3].price);
							/*
							if (i1 < 1) {
								double te = 0;
								
								cout << optionsPortfolio[i2][i3].priceEvaluatedByVG << "   " << optionsPortfolio[i2][i3].price << endl;
								cout << optionsPortfolio[i2][i3].underlyingPrice << "   " << optionsPortfolio[i2][i3].MCunderlyingPrice << "   " << optionsPortfolio[i2][i3].strike << "  " << optionsPortfolio[i2][i3].optionType << endl;
								cout << "DIF: " << (optionsPortfolio[i2][i3].priceEvaluatedByVG - optionsPortfolio[i2][i3].price) << "  Position: " << optionsPortfolio[i2][i3].position << endl;
								te = optionsPortfolio[i2][i3].position * (optionsPortfolio[i2][i3].priceEvaluatedByVG - optionsPortfolio[i2][i3].price);
								total = total + te;
								cout << "Change: " << te << "   Total Change: " << total << endl;
								cout << "**************************************" << endl;
							}
							*/
						}
					}
				}
			}
		}
	}


	sort(simulatedPortPriceDif, simulatedPortPriceDif + numOfPath);
	string outputFileName_PL = "unrealistic_P&L_MonteCarloVaR_VG.csv";

	ofstream output_PL(outputFileName_PL);
	output_PL << "MC simulation," << "P&L of the portfolio," << "P&L of the future" << endl;
	for (int i1 = 0; i1 < numOfPath; ++i1) {
		output_PL << "path" + to_string(i1 + 1) << "," << simulatedPortPriceDif[i1] << "," << futurediff[i1] << endl;
		cout << i1 + 1 << " : " << simulatedPortPriceDif[i1] << endl;
	}
	output_PL.close();
	cout << "Value at Risk: " << -simulatedPortPriceDif[int(numOfPath * 0.01)];
}