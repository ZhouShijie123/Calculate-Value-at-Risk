#ifndef NORMAL_DIST
#define NORMAL_DIST
#include <string>
using namespace std;

double invertingStandardNormalDistribution(double p);

void simulationOfMultinormal(double *mean, double **CovMat, string outputFileName, int nCol, int numOfPath, string *title);


#endif
