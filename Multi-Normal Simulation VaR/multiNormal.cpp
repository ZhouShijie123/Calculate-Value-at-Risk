#include "stdafx.h"
#include "time.h"
#include "matrixops.h"
#include "normalDist.h"
#include <time.h>
#include <io.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdlib.h>

using namespace std;

void simulationOfMultinormal(double *mean, double **CovMat, string outputFileName, int nCol, int numOfPath, string *title) {
	srand(time(NULL));
	double **simuDataMatrix1;
	simuDataMatrix1 = new double *[numOfPath];
	for (int i1 = 0; i1 < numOfPath; ++i1) {
		simuDataMatrix1[i1] = new double[nCol];
	}
	for (int i1 = 0; i1 < numOfPath; ++i1) {
		for (int i2 = 0; i2 < nCol; ++i2) {
			simuDataMatrix1[i1][i2] = invertingStandardNormalDistribution((double)rand() / (double)(RAND_MAX + 1.0));
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
	outputFile.open(outputFileName);
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

	for (int i1 = 0; i1 < numOfPath; ++i1) {
		delete simuDataMatrix1[i1];
		delete simuDataMatrix2[i1];
	}
	for (int i1 = 0; i1 < nCol; ++i1) {
		delete squareRootOfCov[i1];
		delete transposeOfSquareRoot[i1];
	}

	delete[] simuDataMatrix1;
	delete[] simuDataMatrix2;
	delete[] squareRootOfCov;
	delete[] transposeOfSquareRoot;

}

void test_multiNormal() {

	int ncol = 5;
	double* mean = new double[5]{ 1,2,3,4,5 };
	double** covMat = new double*[ncol];
	int counter = 1;
	for (int i = 0; i < ncol; i++) {
		covMat[i] = new double[ncol];
		for (int j = 0; j < ncol; j++) {
			if (i == j) {
				covMat[i][j] = i + 1;
			}
			else {
				covMat[i][j] = 0;
			}
		}
	}
	for (int i = 0; i < ncol; i++) {
		for (int j = 0; j < ncol; j++) {
			cout << covMat[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	double** S = new double*[ncol];
	for (int i = 0; i < ncol; i++) {
		S[i] = new double[ncol];
		for (int j = 0; j < ncol; j++) {
			if (i == j) {
				S[i][j] = 1;
			}
			else {
				S[i][j] = 0;
			}
		}
	}
	for (int i = 0; i < ncol; i++) {
		for (int j = 0; j < ncol; j++) {
			cout << S[i][j] << " ";
		}
		cout << endl;
	}


	string outputFileName = "output222017.csv";
	int numPath = 100;
	string* title = new string[5]{ "1","2","3","4","5" };

	simulationOfMultinormal(mean, covMat, outputFileName, ncol, numPath, title);

}
