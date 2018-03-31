#include "stdafx.h"
#include <iostream>
#include <string>
#include "matrixops.h"
using namespace std;

void EMfixData(double **myData, int **ifExist, int nRow, int nCol, double **fixedData, double **estiCovMat, double *estiMean) {
	double ***covMatrix;
	double **mean;
	int maxIter = 500;
	int iter = 0;
	double valueChange = 1.0;
	double threshold = 0.0001;

	covMatrix = new double **[maxIter + 1];
	mean = new double *[maxIter + 1];

	for (int i1 = 0; i1 < maxIter + 1; ++i1) {
		covMatrix[i1] = new double *[nCol];
	}
	for (int i1 = 0; i1 < maxIter + 1; ++i1) {
		for (int i2 = 0; i2 < nCol; ++i2) {
			covMatrix[i1][i2] = new double[nCol];
		}
	}
	for (int i1 = 0; i1 < maxIter + 1; ++i1) {
		mean[i1] = new double[nCol];
	}
	for (int i1 = 0; i1 < maxIter + 1; ++i1) {
		for (int i2 = 0; i2 < nCol; ++i2) {
			mean[i1][i2] = 0.0;
		}
	}

	for (int i1 = 0; i1 < nCol; ++i1) {
		for (int i2 = 0; i2 < nCol; ++i2) {
			if (i1 == i2) {
				covMatrix[0][i1][i2] = 1.0;
			}
			else {
				covMatrix[0][i1][i2] = 0.0;
			}
		}
	}
	for (int i1 = 1; i1 < maxIter + 1; ++i1) {
		for (int i2 = 0; i2 < nCol; ++i2) {
			for (int i3 = 0; i3 < nCol; ++i3) {
				covMatrix[i1][i2][i3] = 0.0;
			}
		}
	}
	/*

	EM algorithm

	*/

	while ((iter < maxIter) && (valueChange > threshold)) {
		valueChange = 0.0;
		int **locOfValue;
		double **condiExpectValue;
		double ***condiCovMat;
		locOfValue = new int *[nRow];
		for (int i1 = 0; i1 < nRow; ++i1) {
			int numOfValue = nCol;
			for (int i2 = 0; i2 < nCol; ++i2) {
				if (ifExist[i1][i2] == 1) {
					numOfValue -= 1;
				}
			}
			locOfValue[i1] = new int[numOfValue];
			int loc = 0;
			for (int i2 = 0; i2 < nCol; ++i2) {
				if (ifExist[i1][i2] == 0) {
					locOfValue[i1][loc] = i2;
					++loc;
				}
			}
		}
		condiExpectValue = new double *[nRow];
		for (int i1 = 0; i1 < nRow; ++i1) {
			condiExpectValue[i1] = new double[nCol];
		}
		condiCovMat = new double **[nRow];
		for (int i1 = 0; i1 < nRow; ++i1) {
			condiCovMat[i1] = new double *[nCol];
		}
		for (int i1 = 0; i1 < nRow; ++i1) {
			for (int i2 = 0; i2 < nCol; ++i2) {
				condiCovMat[i1][i2] = new double[nCol];
			}
		}
		for (int i1 = 0; i1 < nRow; ++i1) {
			double **covMatrix21;
			double **covMatrix12;
			double **covMatrix22;
			double *existValue;
			double *existAveValue;

			int numOfValue = nCol;
			for (int i2 = 0; i2 < nCol; ++i2) {
				if (ifExist[i1][i2] == 1) {
					numOfValue = numOfValue - 1;
				}
			}
			for (int i2 = 0; i2 < nCol; ++i2) {
				condiExpectValue[i1][i2] = mean[iter][i2];
			}
			covMatrix21 = new double *[numOfValue];
			for (int i2 = 0; i2 < numOfValue; i2++) {
				covMatrix21[i2] = new double[nCol];
			}
			covMatrix12 = new double *[nCol];
			for (int i2 = 0; i2 < nCol; ++i2) {
				covMatrix12[i2] = new double[numOfValue];
			}
			covMatrix22 = new double *[numOfValue];
			for (int i2 = 0; i2 < numOfValue; ++i2) {
				covMatrix22[i2] = new double[numOfValue];
			}
			for (int i2 = 0; i2 < numOfValue; ++i2) {
				for (int i3 = 0; i3 < nCol; ++i3) {
					covMatrix21[i2][i3] = covMatrix[iter][locOfValue[i1][i2]][i3];
				}
			}
			for (int i2 = 0; i2 < nCol; ++i2) {
				for (int i3 = 0; i3 < numOfValue; ++i3) {
					covMatrix12[i2][i3] = covMatrix[iter][i2][locOfValue[i1][i3]];
				}
			}
			for (int i2 = 0; i2 < numOfValue; ++i2) {
				for (int i3 = 0; i3 < numOfValue; ++i3) {
					covMatrix22[i2][i3] = covMatrix[iter][locOfValue[i1][i2]][locOfValue[i1][i3]];
				}
			}

			existValue = new double[numOfValue];
			for (int i2 = 0; i2 < numOfValue; ++i2) {
				existValue[i2] = myData[i1][locOfValue[i1][i2]];
			}
			existAveValue = new double[numOfValue];
			for (int i2 = 0; i2 < numOfValue; ++i2) {
				existAveValue[i2] = mean[iter][locOfValue[i1][i2]];
			}

			double **inverseOfCov;
			inverseOfCov = new double *[numOfValue];
			for (int i2 = 0; i2 < numOfValue; ++i2) {
				inverseOfCov[i2] = new double[numOfValue];
			}

			calInverseMatrix(covMatrix22, inverseOfCov, numOfValue);

			double **multiMatrix1;

			multiMatrix1 = new double *[nCol];
			for (int i2 = 0; i2 < nCol; ++i2) {
				multiMatrix1[i2] = new double[numOfValue];
			}
			matrixMultiply(covMatrix12, inverseOfCov, multiMatrix1, nCol, numOfValue, numOfValue);
			for (int i2 = 0; i2 < nCol; ++i2) {
				for (int i3 = 0; i3 < numOfValue; ++i3) {
					condiExpectValue[i1][i2] += multiMatrix1[i2][i3] * (existValue[i3] - existAveValue[i3]);
				}
			}

			double **multiMatrix2;
			multiMatrix2 = new double *[nCol];
			for (int i2 = 0; i2 < nCol; ++i2) {
				multiMatrix2[i2] = new double[nCol];
			}

			matrixMultiply(multiMatrix1, covMatrix21, multiMatrix2, nCol, numOfValue, nCol);
			for (int i2 = 0; i2 < nCol; ++i2) {
				for (int i3 = 0; i3 < nCol; ++i3) {
					condiCovMat[i1][i2][i3] = covMatrix[iter][i2][i3] - multiMatrix2[i2][i3];
				}
			}

			for (int i2 = 0; i2 < nCol; ++i2) {
				for (int i3 = 0; i3 < nCol; ++i3) {
					condiCovMat[i1][i2][i3] += (condiExpectValue[i1][i2]) * (condiExpectValue[i1][i3]);
				}
			}

			for (int i2 = 0; i2 < numOfValue; ++i2) {
				delete covMatrix21[i2];
			}
			delete[] covMatrix21;

			for (int i2 = 0; i2 < nCol; ++i2) {
				delete covMatrix12[i2];
			}
			delete[] covMatrix12;

			for (int i2 = 0; i2 < numOfValue; ++i2) {
				delete covMatrix22[i2];
			}
			delete[] covMatrix22;

			for (int i2 = 0; i2 < numOfValue; ++i2) {
				delete inverseOfCov[i2];
			}
			delete[] inverseOfCov;

			for (int i2 = 0; i2 < nCol; ++i2) {
				delete multiMatrix1[i2];
			}
			delete[] multiMatrix1;

			for (int i2 = 0; i2 < nCol; ++i2) {
				delete multiMatrix2[i2];
			}
			delete[] multiMatrix2;

			delete existValue;
			delete existAveValue;

		}

		for (int i1 = 0; i1 < nRow; ++i1) {
			for (int i2 = 0; i2 < nCol; ++i2) {
				fixedData[i1][i2] = condiExpectValue[i1][i2];
			}
		}
		for (int i1 = 0; i1 < nCol; ++i1) {
			for (int i2 = 0; i2 < nRow; ++i2) {
				mean[iter + 1][i1] += 1.0 / nRow * condiExpectValue[i2][i1];
			}
		}
		for (int i1 = 0; i1 < nCol; ++i1) {
			for (int i2 = 0; i2 < nCol; ++i2) {
				for (int i3 = 0; i3 < nRow; ++i3) {
					covMatrix[iter + 1][i1][i2] += 1.0 / nRow * condiCovMat[i3][i1][i2];
				}
			}
		}

		for (int i1 = 0; i1 < nCol; ++i1) {
			for (int i2 = 0; i2 < nCol; ++i2) {
				valueChange += fabs(covMatrix[iter + 1][i1][i2] - covMatrix[iter][i1][i2]);
			}
		}

		iter += 1;

		for (int i1 = 0; i1 < nRow; ++i1) {
			delete locOfValue[i1];
		}
		for (int i1 = 0; i1 < nRow; ++i1) {
			delete condiExpectValue[i1];
		}

		for (int i1 = 0; i1 < nRow; ++i1) {
			for (int i2 = 0; i2 < nCol; ++i2) {
				delete condiCovMat[i1][i2];
			}
		}
		for (int i1 = 0; i1 < nRow; ++i1) {
			delete[] condiCovMat[i1];
		}



		delete[] locOfValue;
		delete[] condiExpectValue;
		delete[] condiCovMat;

	}

	for (int i1 = 0; i1 < nCol; ++i1) {
		for (int i2 = 0; i2 < nCol; ++i2) {
			covMatrix[iter][i1][i2] -= mean[iter][i1] * mean[iter][i2];
		}
	}
	/*
	cout << "iter" << endl;
	cout << iter << endl;

	cout << " mean " << endl;*/
	for (int i1 = 0; i1 < nCol; ++i1) {
		estiMean[i1] = mean[iter][i1];
	}
	/*
	for(int i1=0; i1 < nCol; ++i1){
	cout << mean[iter][i1] << " , ";
	}
	cout << " " << endl;

	cout << "covariance matrix" << endl;
	*/

	for (int i1 = 0; i1 < nCol; ++i1) {
		for (int i2 = 0; i2 < nCol; ++i2) {
			estiCovMat[i1][i2] = covMatrix[iter][i1][i2];
		}
	}
	/*
	for(int i1=0; i1 < nCol; ++i1){
	for(int i2=0; i2 < nCol; ++i2){
	cout << estiCovMat[i1][i2] << " , ";
	}
	cout << " " << endl;
	}
	*/
	for (int i1 = 0; i1 < maxIter + 1; ++i1) {
		for (int i2 = 0; i2 < nCol; ++i2) {
			delete covMatrix[i1][i2];
		}
	}
	for (int i1 = 0; i1 < maxIter + 1; ++i1) {
		delete[] covMatrix[i1];
	}


	for (int i1 = 0; i1 < maxIter + 1; ++i1) {
		delete mean[i1];
	}
	delete[] mean;
	delete[] covMatrix;
}
