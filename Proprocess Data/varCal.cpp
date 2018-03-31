#include "stdafx.h"
#include "varCal.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

void CovCal::getnRowAndnCol(string inputFileName){
	//Initialize nRow and nCol
	nRow = -1; //not counting header
	nCol = -2; //not counting first and last column
	ifstream inputFile;
	inputFile.open(inputFileName.c_str());
	string line;
	while (getline(inputFile, line)){
		nRow++;
	}
	inputFile.clear();
	inputFile.seekg(0);
	string word;
	getline(inputFile, line);
	stringstream lineToProcess(line);
	while (getline(lineToProcess, word, ',')){
		nCol++;
	}
	inputFile.close();
}


CovCal::CovCal(string inputFileName){
	//N.B: rows in the file are sorted from oldest date to most recent(today)
	//Declarations
	getnRowAndnCol(inputFileName);
	date = new string[nRow];
	title = new string[nCol];
	OriginalData = new double*[nRow];
	Missing = new int*[nRow];
	returns = new double*[nRow-1]; //for n days we have n-1 returns
	fixedData = new double*[nRow-1]; //??
	returnMissing = new int *[nRow-1];
	todayData = new double[nCol];
	estiMean = new double[nCol];
	estiCovMat = new double *[nCol]; 
	//Initializations
	for (int i1 = 0; i1 < nCol; ++i1) {
		estiCovMat[i1] = new double[nCol];
	}
	for (int i1 = 0; i1 < nRow-1; ++i1){
		returns[i1]    = new double[nCol];
		fixedData[i1] = new double[nCol];
		returnMissing[i1]   = new int[nCol];
	}
	for (int i1 = 0; i1 < nRow; ++i1) {
		OriginalData[i1] = new double[nCol];
		Missing[i1] = new int[nCol];
	}
	//Looping through lines of the file
	ifstream inputFile;
	inputFile.open(inputFileName.c_str());
	string line;
	string word;
	getline(inputFile, line);
	stringstream lineToProcess(line);
	getline(lineToProcess, word, ',');
	//Filling Contract titles
	for (int col = 0; col< nCol; ++col){
		getline(lineToProcess, word, ',');
		title[col] = word;
	}
	//Fill OriginalData,date and Missing
	for (int row = 0; row < nRow; ++row){
		getline(inputFile, line);
		stringstream lineToProcess(line);
		getline(lineToProcess, word, ',');
		date[row] = word;
		for (int col = 0; col < nCol; ++col){
			getline(lineToProcess, word, ',');
			if(word.empty())
			{
				OriginalData[row][col] = 0.0;
				Missing[row][col] = 1;
			}
			else {
				OriginalData[row][col] = atof(word.c_str());
				Missing[row][col] = 0;
			}
		}
	}
	inputFile.close();
	//Last row  corresponds to today's data
	for(int col = 0; col < nCol; ++col){
		todayData[col] = OriginalData[nRow-1][col];
	}
	//Getting returns
	for (int row = 0; row< nRow-1; ++row){
		for (int col = 0; col < nCol; ++col) {
			if ((Missing[row][col] == 1) || (Missing[row + 1][col] == 1)) {
				returns[row][col] = 0;
				returnMissing[row][col] = 1;
			}
			else {
				returns[row][col] = (OriginalData[row + 1][col] - OriginalData[row][col]) / OriginalData[row][col];
				returnMissing[row][col] = 0;
			}
		}
	}
}
void CovCal::getData() {
	cout << "time" << " ";
	for (int i1 = 0; i1 < nCol; ++i1) {
		cout << title[i1] << " ";
	}
	cout << " " << endl;
	for (int i1 = 0; i1 < nRow - 1; ++i1) {
		cout << date[i1 + 1] << " ";
		for (int i2 = 0; i2 < nCol; ++i2) {
			if (returnMissing[i1][i2] == 0) {
				cout << returns[i1][i2] << " , ";
			}
			else {
				cout << "miss" << " , ";
			}
		}
		cout << " " << endl;
	}
}

void CovCal::getTodayData() {
	int position;
	for (int i1 = 0; i1 < nCol; ++i1) {
		position = 0;
		while (position < nRow) {
			if (Missing[nRow - 1 - position][i1] == 0) {
				todayData[i1] = OriginalData[nRow - 1 - position][i1];
				for (int i2 = 0; i2 < position; ++i2) {
					cout << "mydata: " << fixedData[nRow - 1 - i2 - 1][i1] << endl;
					todayData[i1] = todayData[i1] * (1.0 + fixedData[nRow - 1 - i2 - 1][i1]);
				}
				break;
			}
			else {
				position += 1;
			}
		}
	}
}

void CovCal::getMeanandVariance() {
	cout << "     mean     " << endl;
	for (int i1 = 0; i1 < nCol; ++i1) {
		cout << estiMean[i1] << " , ";
	}
	cout << " " << endl;
	cout << "     covariance matrix   " << endl;
	for (int i1 = 0; i1 < nCol; ++i1) {
		for (int i2 = 0; i2 < nCol; ++i2) {
			cout << estiCovMat[i1][i2] << " , ";
		}
		cout << " " << endl;
	}
}

CovCal::~CovCal() {
	for (int i1 = 0; i1 < nRow ; ++i1) {
		delete OriginalData[i1];
		delete Missing[i1];
	}
	for (int i1 = 0; i1 < nRow-1; ++i1) {
		delete returns[i1];
		delete fixedData[i1];
		delete returnMissing[i1];
	}
	for (int i1 = 0; i1 < nCol; ++i1) {
		delete estiCovMat[i1];
	}
	delete todayData;
	delete estiMean;
	delete[] OriginalData;
	delete[] Missing;
	delete[] estiCovMat;
	delete[] returns;
	delete[] fixedData;
	delete[] returnMissing;
	delete[] title;
	delete[] date;
}
















