//#include <sys/time.h>
#include <stdio.h>
//#include <unistd.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <numeric>
#include <vector>
#include <map>
#include <math.h> 
#include<cstring>
#include <stdio.h>

using namespace std;

void matrixMultiply(double **matrix1, double **matrix2, double **resultMatrix, int nRow1, int nRow2, int nCol);

void calInverseMatrix(double **matrix1, double **inverse, int n);

double calDeterminant(double **matrix, int n);

void getSquareRootofSymMat(double **matrix, double **squareRoot, int nRow);

double invertingStandardNormalDistribution(double p);

void simulationOfMultinormal(double *mean, double **CovMat, string outputFileName, int nCol, int numOfPath, string *title);

double cumFuncOfNormalDistribution(double mean, double vol, double location);

double calPriceForForexOption(double s0, double strike, double rd, double rf, double mature, double vol, string optionType);

void EMfixData(double **myData, int **ifExist, int nRow, int nCol, double **fixedData, double **estiCovMat, double *estiMean);

double momentGenerate(string element1, string element2, map<string, double> mean, map<string, map<string, double> > cov, int order1, int order2);
