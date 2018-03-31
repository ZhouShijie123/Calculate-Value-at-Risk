#ifndef MATRIX_OPS
#define MATRIX_OPS
#include <string>
#include <map>
using namespace std;
void matrixMultiply(double **matrix1, double **matrix2, double **resultMatrix, int nRow1, int nRow2, int nCol);

void calInverseMatrix(double **matrix1, double **inverse, int n);

void getSquareRootofSymMat(double **matrix, double **squareRoot, int nRow);

double momentGenerate(string element1, string element2, map<string, double> mean, map<string, map<string, double> > cov, int order1, int order2);
#endif
