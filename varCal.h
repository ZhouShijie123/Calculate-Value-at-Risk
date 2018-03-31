#ifndef VARCAL
#define VARCAL
#include <string>
using namespace std;

class CovCal{
public:
	double **returns;
	int **returnMissing;

	void getnRowAndnCol(string inputFileName);

	int nRow;
	int nCol;

	double **OriginalData;
	int **Missing;
	string *date;
	string *title;
	double *todayData;
	double **fixedData;
	double **estiCovMat;
	double *estiMean;

	CovCal(string inputFileName);
	void getData();
	void getTodayData();
	void getMeanandVariance();
	~CovCal();

};
#endif