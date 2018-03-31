#include "stdafx.h"
#include <string>
#include <map>
#include <iostream>
using namespace std;

//Matrix multiplication
void matrixMultiply(double **matrix1, double **matrix2, double **resultMatrix, int nRow1, int nRow2, int nCol) {

	for (int i1 = 0; i1 < nRow1; ++i1) {
		for (int i2 = 0; i2 < nCol; ++i2) {
			resultMatrix[i1][i2] = 0;
		}
	}
	for (int i1 = 0; i1 < nRow1; ++i1) {
		for (int i2 = 0; i2 < nCol; ++i2) {
			for (int i3 = 0; i3 < nRow2; ++i3) {
				resultMatrix[i1][i2] += matrix1[i1][i3] * matrix2[i3][i2];
			}
		}
	}
}

//Calculate the inverse of a matrix
void calInverseMatrix(double **matrix1, double **inverse, int n) {

	double **matrix;
	matrix = new double *[n];
	for (int i1 = 0; i1 < n; ++i1) {
		matrix[i1] = new double[n];
	}

	for (int i1 = 0; i1 < n; ++i1) {
		for (int i2 = 0; i2 < n; ++i2) {
			matrix[i1][i2] = matrix1[i1][i2];
		}
	}
	for (int i1 = 0; i1 < n; ++i1) {
		for (int i2 = 0; i2 < n; ++i2) {
			if (i1 == i2) {
				inverse[i1][i2] = 1.0;
			}
			else {
				inverse[i1][i2] = 0.0;
			}
		}
	}

	double divider;
	double deno_coef;
	for (int i1 = 0; i1 < n; ++i1) {
		divider = matrix[i1][i1];
		for (int i2 = 0; i2 < n; ++i2) {
			matrix[i1][i2] = matrix[i1][i2] / divider;
			inverse[i1][i2] = inverse[i1][i2] / divider;
		}
		for (int i3 = i1 + 1; i3 < n; ++i3) {
			deno_coef = matrix[i3][i1];
			for (int i4 = 0; i4 < n; ++i4) {
				matrix[i3][i4] = matrix[i3][i4] - deno_coef * matrix[i1][i4];
				inverse[i3][i4] = inverse[i3][i4] - deno_coef * inverse[i1][i4];
			}
		}
	}//row downward

	for (int i1 = n - 1; i1 > -1; --i1) {
		for (int i2 = i1 - 1; i2 > -1; --i2) {
			deno_coef = matrix[i2][i1];
			for (int i3 = 0; i3 < n; ++i3) {
				matrix[i2][i3] = matrix[i2][i3] - deno_coef * matrix[i1][i3];
				inverse[i2][i3] = inverse[i2][i3] - deno_coef * inverse[i1][i3];
			}
		}
	}  //row upward
	for (int i1 = 0; i1 < n; ++i1) {
		delete matrix[i1];
	}
	delete[] matrix;
}


//Compute the Cholesky decomposition of a a positive semi-definite matrix
void getSquareRootofSymMat(double **matrix, double **squareRoot, int nRow) {
	double **Lmatrix;
	double **LinverseMatrix;
	double **LtransposeOfInverse;
	double **copyMatrix;
	Lmatrix = new double *[nRow];
	LinverseMatrix = new double *[nRow];
	LtransposeOfInverse = new double *[nRow];
	copyMatrix = new double *[nRow];
	for (int i1 = 0; i1 < nRow; ++i1) {
		Lmatrix[i1] = new double[nRow];
		LinverseMatrix[i1] = new double[nRow];
		LtransposeOfInverse[i1] = new double[nRow];
		copyMatrix[i1] = new double[nRow];
	}
	for (int i1 = 0; i1 < nRow; ++i1) {
		for (int i2 = 0; i2 < nRow; ++i2) {
			if (i1 == i2) {
				LinverseMatrix[i1][i2] = 1.0;
			}
			else {
				LinverseMatrix[i1][i2] = 0.0;
			}
		}
	}
	for (int i1 = 0; i1 < nRow; ++i1) {
		for (int i2 = 0; i2 < nRow; ++i2) {
			copyMatrix[i1][i2] = matrix[i1][i2];
		}
	}
	double divider;
	double deno_coef;
	for (int i1 = 0; i1 < nRow; ++i1) {
		divider = copyMatrix[i1][i1];
		for (int i3 = i1 + 1; i3 < nRow; ++i3) {
			deno_coef = copyMatrix[i3][i1];
			for (int i4 = 0; i4 < nRow; ++i4) {
				copyMatrix[i3][i4] = copyMatrix[i3][i4] - deno_coef / divider * copyMatrix[i1][i4];
				LinverseMatrix[i3][i4] = LinverseMatrix[i3][i4] - deno_coef / divider * LinverseMatrix[i1][i4];
			}
		}
	}

	for (int i1 = 0; i1 < nRow; ++i1) {
		for (int i2 = 0; i2 < nRow; ++i2) {
			LtransposeOfInverse[i1][i2] = LinverseMatrix[i2][i1];
		}
	}

	double **multiMatrix1;
	double **multiMatrix2;

	multiMatrix1 = new double *[nRow];
	multiMatrix2 = new double *[nRow];

	for (int i1 = 0; i1 < nRow; ++i1) {
		multiMatrix1[i1] = new double[nRow];
		multiMatrix2[i1] = new double[nRow];
	}
	calInverseMatrix(LinverseMatrix, Lmatrix, nRow);

	matrixMultiply(LinverseMatrix, matrix, multiMatrix1, nRow, nRow, nRow);
	matrixMultiply(multiMatrix1, LtransposeOfInverse, multiMatrix2, nRow, nRow, nRow);

	for (int i1 = 0; i1 < nRow; ++i1) {
		multiMatrix2[i1][i1] = sqrt(multiMatrix2[i1][i1]);
	}

	matrixMultiply(Lmatrix, multiMatrix2, squareRoot, nRow, nRow, nRow);

	for (int i1 = 0; i1 < nRow; ++i1) {
		delete Lmatrix[i1];
		delete LinverseMatrix[i1];
		delete LtransposeOfInverse[i1];
		delete copyMatrix[i1];
		delete multiMatrix1[i1];
		delete multiMatrix2[i1];
	}
	delete[] Lmatrix;
	delete[] LinverseMatrix;
	delete[] LtransposeOfInverse;
	delete[] copyMatrix;
	delete[] multiMatrix1;
	delete[] multiMatrix2;
}

double momentGenerate(string element1, string element2, map<string, double> mean, map<string, map<string, double> > cov, int order1, int order2){
	//The following function computer E[Xi*Xj]
	double *meanVector;
	meanVector = new double[1];
	double **covMatrix;
	covMatrix = new double *[2];
	for (int i1 = 0; i1 < 2; ++i1){
		covMatrix[i1] = new double[2];
	}
	meanVector[0] = mean[element1];
	meanVector[1] = mean[element2];
	covMatrix[0][0] = cov[element1][element1];
	covMatrix[0][1] = cov[element1][element2];
	covMatrix[1][0] = cov[element2][element1];
	covMatrix[1][1] = cov[element2][element2];

	double moment;
	//Deallocation of memory
	if (order1 == 0) {
		if (order2 == 0) {

			delete meanVector;
			for (int i1 = 0; i1 < 2; ++i1) {
				delete covMatrix[i1];
			}
			delete[] covMatrix;

			return 1;
		}
		else {
			if (order2 == 1) {
				moment = meanVector[1];

				delete meanVector;
				for (int i1 = 0; i1 < 2; ++i1) {
					delete covMatrix[i1];
				}
				delete[] covMatrix;

				return moment;
			}
			else {
				if (order2 == 2) {
					moment = covMatrix[1][1] + pow(meanVector[1], 2);

					delete meanVector;
					for (int i1 = 0; i1 < 2; ++i1) {
						delete covMatrix[i1];
					}
					delete[] covMatrix;

					return moment;
				}
				else {
					return 0.0;
				}
			}
		}
	}
	else {
		if (order2 == 0) {
			if (order1 == 1) {
				moment = meanVector[0];

				delete meanVector;
				for (int i1 = 0; i1 < 2; ++i1) {
					delete covMatrix[i1];
				}
				delete[] covMatrix;

				return moment;
			}
			else {
				if (order2 == 2) {
					moment = covMatrix[0][0] + pow(meanVector[0], 2);

					delete meanVector;
					for (int i1 = 0; i1 < 2; ++i1) {
						delete covMatrix[i1];
					}
					delete[] covMatrix;

					return moment;
				}
				else {
					return 0.0;
				}
			}
		}
		else {
			if (order1 == 1) {
				if (order2 == 1) {

					moment = covMatrix[0][1] + meanVector[0] * meanVector[1];

					delete meanVector;
					for (int i1 = 0; i1 < 2; ++i1) {
						delete covMatrix[i1];
					}
					delete[] covMatrix;

					return moment;
				}
				else {
					if (order2 == 2) {

						moment = 2 * covMatrix[0][1] * meanVector[1] + covMatrix[1][1] * meanVector[0] + meanVector[0] * pow(meanVector[1], 2);

						delete meanVector;
						for (int i1 = 0; i1 < 2; ++i1) {
							delete covMatrix[i1];
						}
						delete[] covMatrix;

						return moment;
					}
					else {
						return 0.0;
					}
				}
			}
			else {
				if (order1 == 2) {
					if (order2 == 1) {

						moment = 2 * covMatrix[0][1] * meanVector[0] + covMatrix[0][0] * meanVector[1] + pow(meanVector[0], 2) * meanVector[1];

						delete meanVector;
						for (int i1 = 0; i1 < 2; ++i1) {
							delete covMatrix[i1];
						}
						delete[] covMatrix;

						return moment;
					}
					else {
						if (order2 == 2) {

							moment = 2 * pow(covMatrix[0][1], 2) + 2 * covMatrix[0][1] * meanVector[0] * meanVector[1] + covMatrix[0][0] * covMatrix[1][1] + covMatrix[0][0] * pow(meanVector[1], 2) + 2 * covMatrix[0][1] * meanVector[0] * meanVector[1] + pow(meanVector[0], 2) * covMatrix[1][1] + pow(meanVector[0] * meanVector[1], 2);

							delete meanVector;
							for (int i1 = 0; i1 < 2; ++i1) {
								delete covMatrix[i1];
							}
							delete[] covMatrix;

							return moment;
						}
						else {
							return 0.0;
						}
					}
				}
				else {
					return 0.0;
				}
			}
		}
	}
}

void test_matrixops()
{
	//Inverse
	int n = 3;
	double** A = new double*[3];
	int counter = 1;
	for (int i = 0; i < n; i++){
		A[i] = new double[3];
		for (int j = 0; j < n; j++) {
			if (i == j) {
				A[i][j] = 2;
			}
			else {
				A[i][j] = 0;
			}
			counter++;
		}
	}
	A[0][1] = 2;
	A[0][2] = 2;
	double** I =new double*[3];

	for (int i = 0; i < n; i++){
		I[i] = new double[3];
		for (int j = 0; j < n; j++) {
			if (i == j) {
				I[i][j] = 1;
			}
			else {
				I[i][j] = 0;
			}
			counter++;
		}
	}

	//Warning : there are issues with determinant and inverse
	calInverseMatrix(A, I,n);
	for(int i = 0; i < n; i++) {
		for (int j = 0;j < n; j++) {
			cout << I[i][j] << " ";
		}
		cout << endl;
	}
}