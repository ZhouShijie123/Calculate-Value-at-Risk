#include "stdafx.h"
#include "matrixops.h"
#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

void EMfixData(double **myData, int **ifExist, int T, int K, double **fixedData, double **estiCovMat, double *estiMean,string* title,string outputMeanFile,string outputCovFile){
	cout << T << endl;
	//Parameters
	int maxIter = 500;
	double valueChange = 1.0;
	double threshold = 0.0001;
	//Declarations
	int iter = 0;
	double ***covMatrix; //iterx(ncolxncol) tensor
	double **mean; //iterxncol 
	covMatrix = new double**[maxIter + 1]; //500*40*40 = 3.2 mb data size, can largely be optimized
	mean = new double*[maxIter + 1];
	//CovMatrix alloc
	for (int i1 = 0; i1 < maxIter + 1; ++i1){
		covMatrix[i1] = new double*[K];
	}
	for (int i1 = 0; i1 < maxIter + 1; ++i1){
		for (int i2 = 0; i2 < K; ++i2){
			covMatrix[i1][i2] = new double[K];
		}
	}
	//Mean alloc
	for (int i1 = 0; i1 < maxIter + 1; ++i1){
		mean[i1] = new double[K];
	}
	//Mean init
	for (int i1 = 0; i1 < maxIter + 1; ++i1){
		for (int i2 = 0; i2 < K; ++i2){
			mean[i1][i2] = 0.0;
		}
	}
	//CovMatrix init for first itertion as Identity matrix -> bad idea as it will slow down convergence A LOT
	//Better idea: Initialize cov matrix using sample cov on complete lines (from 1 to 249)
	for (int i1 = 0; i1 < K; ++i1){
		for (int i2 = 0; i2 < K; ++i2) {
			if (i1 == i2) {
				covMatrix[0][i1][i2] = 1.0;
			}
			else {
				covMatrix[0][i1][i2] = 0.0;
			}
		}
	}
	//CovMatrix init for following iterations
	for (int i1 = 1; i1 < maxIter + 1; ++i1){
		for (int i2 = 0; i2 < K; ++i2) {
			for (int i3 = 0; i3 < K; ++i3) {
				covMatrix[i1][i2][i3] = 0.0;
			}
		}
	}
	//EM algorithm
	//For each iteration
	while ((iter < maxIter) && (valueChange > threshold)){
		valueChange = 0.0;
		int** ZLoc;
		int** RLoc;
		double** condiExpectValue;
		double*** condiCovMat;
		ZLoc = new int*[T];//array containing location of non missing values
		RLoc = new int*[T];//array containing location of missing values
		// Get location of missing and existent
		for (int t = 0; t < T; ++t){
			int nExistent = K; //Number of non missing values in current row (t)
			for (int k = 0; k < K; ++k) {
				nExistent -= ifExist[t][k];
			}
			int nMissing = K - nExistent;
			ZLoc[t] = new int[nExistent];
			RLoc[t] = new int[nMissing];
			int jexistent = 0;
			int jmissing = 0;
			for (int k = 0; k < K; ++k) {
				if (ifExist[t][k] == 0) {
					ZLoc[t][jexistent] = k;
					++jexistent;
				}
				else {
					RLoc[t][jmissing] = k;
					++jmissing;
				}
			}
		}
		//Allocation of E[rt|zt,theta_iter] (TxK matrix)
		condiExpectValue = new double *[T]; 
		for (int t = 0; t < T; ++t) {
			condiExpectValue[t] = new double[K];
		}
		//Allocation of Cov[rt|zt,theta_iter] (Tx (KxK) tensor)
		condiCovMat = new double **[T];
		for (int i1 = 0; i1 < T; ++i1) {
			condiCovMat[i1] = new double *[K];
		}
		for (int i1 = 0; i1 < T; ++i1) {
			for (int i2 = 0; i2 < K; ++i2) {
				condiCovMat[i1][i2] = new double[K];
			}
		}
		//For each row/time
		for (int t = 0; t < T; ++t){
			double **Covzr; //Cov(rz)
			double **Covrz; //Cov(zr)
			double **Covzz; //Cov(zz)
			double *existValue;
			double *existAveValue;
			//Count non missing values(again ?)
			int nExistent = K;
			for (int i2 = 0; i2 < K; ++i2) {
					nExistent-=ifExist[t][i2];
			}
			for (int i2 = 0; i2 < K; ++i2) {
				condiExpectValue[t][i2] = mean[iter][i2];
			}
			//Alloc Covzr
			Covrz = new double*[K];
			for (int i2 = 0; i2 < K; ++i2) {
				Covrz[i2] = new double[nExistent];
			}
			//Alloc Covrz
			Covzr = new double *[nExistent];
			for (int i2 = 0; i2 < nExistent; i2++) {
				Covzr[i2] = new double[K];
			}
			//Alloc Covzz
			Covzz = new double *[nExistent];
			for (int i2 = 0; i2 < nExistent; ++i2) {
				Covzz[i2] = new double[nExistent];
			}
			//Fill Covrz
			for (int i2 = 0; i2 < K; ++i2) {
				for (int i3 = 0; i3 < nExistent; ++i3) {
					Covrz[i2][i3] = covMatrix[iter][i2][ZLoc[t][i3]];
				}
			}
			//Fill Covzr
			for (int i2 = 0; i2 < nExistent; ++i2) {
				for (int i3 = 0; i3 < K; ++i3) {
					Covzr[i2][i3] = covMatrix[iter][ZLoc[t][i2]][i3];
				}
			}
			//Fill Covzz
			for (int i2 = 0; i2 < nExistent; ++i2) {
				for (int i3 = 0; i3 < nExistent; ++i3){
					Covzz[i2][i3] = covMatrix[iter][ZLoc[t][i2]][ZLoc[t][i3]];
				}
			}
			//M-Step of EM algorithm, compute new estimate
			existValue = new double[nExistent];
			for (int i2 = 0; i2 < nExistent; ++i2) {
				existValue[i2] = myData[t][ZLoc[t][i2]];
			}
			existAveValue = new double[nExistent];
			for (int i2 = 0; i2 < nExistent; ++i2) {
				existAveValue[i2] = mean[iter][ZLoc[t][i2]];
			}
			//Compute covzz^-1
			double **inverseOfCov;
			inverseOfCov = new double *[nExistent];
			for (int i2 = 0; i2 < nExistent; ++i2) {
				inverseOfCov[i2] = new double[nExistent];
			}
			calInverseMatrix(Covzz, inverseOfCov, nExistent);
			//Compute Covrz*Covzz^-1 = CovrzXCovzz_inv
			double **CovrzXCovzz_inv;
			CovrzXCovzz_inv = new double *[K];
			for (int i2 = 0; i2 < K; ++i2) {
				CovrzXCovzz_inv[i2] = new double[nExistent];
			}
			matrixMultiply(Covrz, inverseOfCov, CovrzXCovzz_inv, K, nExistent, nExistent);
			//Conditional expected value E[rt|zt,thetaiter] = mur + Covrz*Covzz_inv*(zt-uz)
			for (int k = 0; k < K; ++k) {
				for (int j = 0; j < nExistent; ++j) {
					condiExpectValue[t][k] += CovrzXCovzz_inv[k][j] * (existValue[j] - existAveValue[j]);
				}
			}
			//Compute Covrz*Covzz^-1*Covzr
			double **CovrzXCovzz_invXCovzr;
			CovrzXCovzz_invXCovzr = new double *[K];
			for (int i2 = 0; i2 < K; ++i2) {
				CovrzXCovzz_invXCovzr[i2] = new double[K];
			}
			matrixMultiply(CovrzXCovzz_inv, Covzr, CovrzXCovzz_invXCovzr, K, nExistent, K);
			//Cov[rtT|zt,theta] = Covrr - Covrz*Covzz^-1*Covzr
			for (int i2 = 0; i2 < K; ++i2) {
				for (int i3 = 0; i3 < K; ++i3) {
					condiCovMat[t][i2][i3] = covMatrix[iter][i2][i3] - CovrzXCovzz_invXCovzr[i2][i3]; //Why not Covrr -  Covrz*Covzz^-1*Covzr, any trick here ?
				}
			}
			//E[rtrtT|z_t,theta_iter] = Cov[rt_T|Zt,theta] + E[rt|zt]*E[rt|zt]T
			for (int i2 = 0; i2 < K; ++i2) {
				for (int i3 = 0; i3 < K; ++i3) {
					condiCovMat[t][i2][i3] += (condiExpectValue[t][i2]) * (condiExpectValue[t][i3]);
				}
			}
			//Deallocation
			for (int i2 = 0; i2 < nExistent; ++i2) {
				delete Covzr[i2];
			}
			delete[] Covzr;

			for (int i2 = 0; i2 < K; ++i2) {
				delete Covrz[i2];
			}
			delete[] Covrz;

			for (int i2 = 0; i2 < nExistent; ++i2) {
				delete Covzz[i2];
			}
			delete[] Covzz;

			for (int i2 = 0; i2 < nExistent; ++i2) {
				delete inverseOfCov[i2];
			}
			delete[] inverseOfCov;

			for (int i2 = 0; i2 < K; ++i2) {
				delete CovrzXCovzz_inv[i2];
			}
			delete[] CovrzXCovzz_inv;

			for (int i2 = 0; i2 < K; ++i2) {
				delete CovrzXCovzz_invXCovzr[i2];
			}
			delete[] CovrzXCovzz_invXCovzr;

			delete existValue;
			delete existAveValue;

		}
		//Compute 
		for (int i1 = 0; i1 < T; ++i1){
			for (int i2 = 0; i2 < K; ++i2){
				fixedData[i1][i2] = condiExpectValue[i1][i2];
			}
		}
		//new estimate of mean
		for (int k = 0; k < K; ++k){
			for (int t = 0; t < T; ++t){
				mean[iter + 1][k] += 1.0 / T * condiExpectValue[t][k];
			}
		}
		//new estimation of CovMatrix
		for (int i1 = 0; i1 < K; ++i1) {
			for (int i2 = 0; i2 < K; ++i2) {
				for (int i3 = 0; i3 < T; ++i3) {
					covMatrix[iter + 1][i1][i2] += 1.0 / T * condiCovMat[i3][i1][i2];
				}
			}
		}
		//Compute value change (Norm 1)
		for (int i1 = 0; i1 < K; ++i1) {
			for (int i2 = 0; i2 < K; ++i2) {
				valueChange += fabs(covMatrix[iter + 1][i1][i2] - covMatrix[iter][i1][i2]);
			}
		}
		iter += 1;
		//Deallocation
		for (int t = 0; t < T; ++t) {
			delete ZLoc[t];
		}
		for (int t = 0; t < T; ++t) {
			delete condiExpectValue[t];
		}

		for (int t = 0;t < T; ++t) {
			for (int k = 0; k < K; ++k) {
				delete condiCovMat[t][k];
			}
		}
		for (int t = 0; t < T; ++t) {
			delete[] condiCovMat[t];
		}
		delete[] ZLoc;
		delete[] condiExpectValue;
		delete[] condiCovMat;

	}
	//End of iterations
	//Cov[rtT|zt,theta_iter] = E[rtrtT|Zt,theta_iter] - E[rt|Zt]*E[rt|zt]T but what has this to do with the below ???
	/*for (int i1 = 0; i1 < nCol; ++i1) {
		for (int i2 = 0; i2 < nCol; ++i2) {
			covMatrix[iter][i1][i2] -= mean[iter][i1] * mean[iter][i2];
		}
	}*/
	//Save mean
	ofstream meanfile;
	meanfile.open(outputMeanFile);
	for (int i1 = 0; i1 < K; ++i1) {
		estiMean[i1] = mean[iter][i1];
	}

	for(int i1=0; i1 < K; ++i1){
		meanfile << mean[iter][i1] << ",";
	}
	meanfile.close();
	//Save fixedData
	ofstream fixedfile;
	fixedfile.open("fixedData.csv");
	//header
	fixedfile << "day" << ',';
	for (int i2 = 0; i2 < K; ++i2){
		fixedfile << title[i2] << ',';
	}
	fixedfile << "\n";
	//core file
	int d = 0; // day counter
	for (int i1 = 0; i1 < T; ++i1) {
		d++;
		fixedfile << d <<',';
		for (int i2 = 0; i2 < K; ++i2) {
			fixedfile << fixedData[i1][i2] << ',';
		}
		fixedfile << "\n";
	}
	//Save covariance
	ofstream covfile;
	covfile.open(outputCovFile);
	for (int i1 = 0; i1 < K; ++i1) {
		for (int i2 = 0; i2 < K; ++i2) {
			estiCovMat[i1][i2] = covMatrix[iter][i1][i2];
		}
	}
	
	for(int i1=0; i1 < K; ++i1){
	for(int i2=0; i2 < K; ++i2){
	covfile << estiCovMat[i1][i2] << ",";
	}
	covfile << "\n";
	}
	covfile.close();

	//Deallocating Memory
	for (int i1 = 0; i1 < maxIter + 1; ++i1) {
		for (int i2 = 0; i2 < K; ++i2) {
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
