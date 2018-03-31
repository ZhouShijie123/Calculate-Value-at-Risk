#include "stdafx.h"
#include "getPortfolio.h"
#include <fstream>

vector<Future> getFuturesPortfolio(){
	//Step1 : Processing data
	string futuresPositionsFile = "Sample CL Futures Position.csv";
	//Futures Portfolio
	io::CSVReader<4> futuresCSV(futuresPositionsFile);
	futuresCSV.read_header(io::ignore_extra_column, "Trader", "Position", "FutPx", "Sett");
	vector <Future> futuresPortfolio = vector <Future>();
	string instrument;
	double position;
	double price;
	double sett;
	string sym;
	while (futuresCSV.read_row(instrument, position, price, sett)) {
		Future futures;
		if (instrument == "Total:") continue;
		if (string::npos == instrument.find_first_of("0123456789")) { sym = instrument; }
		else {
			futures.title = sym + instrument;
			futures.position = position;
			futures.price = price / 100.0;
			futures.gamma = 0;
			futuresPortfolio.push_back(futures);
		}
	}
	return futuresPortfolio;
}

vector<vector<option>> getOptionsPortfolio(map<string, map<double, double> >& impliedVol_grid){
	string optionsPositionsFile = "Sample LO Options Position.csv";
	string volSmileFile = "Vol Smile.csv";
	//Options Portfolio
	io::CSVReader<21> optionsCSV(optionsPositionsFile);
	optionsCSV.read_header(io::ignore_extra_column, "Trader", "F17C", "F17P", "G17C", "G17P", "H17C", "H17P", "J17C", "J17P", "K17C", "K17P", "M17C", "M17P", "N17C", "N17P", "Q17C", "Q17P", "U17C", "U17P", "Z17C", "Z17P");
	string optionTitle[] = { "CLF17", "CLF17", "CLG17", "CLG17", "CLH17", "CLH17", "CLJ17", "CLJ17", "CLK17", "CLK17", "CLM17", "CLM17", "CLN17", "CLN17", "CLQ17", "CLQ17", "CLU17", "CLU17", "CLZ17", "CLZ17" };
	vector<vector<option>> optionsPortfolio;
	int optionRow = 0;
	int line = 0;
	string firstColumn;
	int const nContracts = 20;
	string column[nContracts];
	map<string, int> rowContract;
	map<string, int> colContract;
	
	while (optionsCSV.read_row(firstColumn, column[0], column[1], column[2], column[3], column[4], column[5], column[6], column[7], column[8], column[9], column[10], column[11], column[12], column[13], column[14], column[15], column[16], column[17], column[18], column[19])) {
		//Skip second row(contains total)
		if (firstColumn == "CL") {
			continue;
		}
		else {
			vector <option> row = vector<option>(nContracts);
			
			for (int i1 = 0; i1 < nContracts; ++i1) {
				row[i1].underlying = optionTitle[i1];
				if (i1 % 2 == 0) {
					row[i1].optionType = "call";
				}
				else {
					row[i1].optionType = "put";
				}
				row[i1].strike = atof(firstColumn.c_str()) / 100.0;
				row[i1].position = atof((column[i1]).c_str());
				string contract = row[i1].underlying + row[i1].optionType + firstColumn;
				rowContract[contract] = optionRow;
				colContract[contract] = i1;
			}
			optionRow++;
			optionsPortfolio.push_back(row);
		}
		line++;
	}
	//Volatility Smile
	io::CSVReader<10> volSmileCSV(volSmileFile);
	volSmileCSV.read_header(io::ignore_extra_column, "Opt", "Cntr", "Stk", "Call", "Put", "Vol", "FutPx", "DTE", "BSCDlta", "BSPDlta");
	string Opt;
	string Cntr;
	int Strike;
	double Call;
	double Put;
	double Vol;
	double FutPx;
	double DTE;
	double BSCDlta;
	double BSPDlta;
	//Optimal O(N) search
	while (volSmileCSV.read_row(Opt, Cntr, Strike, Call, Put, Vol, FutPx, DTE, BSCDlta, BSPDlta)){
				string key = Opt + Cntr + "call" + to_string(Strike);
				if (rowContract.count(key) == 1) {
					impliedVol_grid[Opt + Cntr][(Strike - FutPx) / 100.0] = Vol / 100.0;
					string keyCall = Opt + Cntr + "call" + to_string(Strike);
					string keyPut = Opt + Cntr + "put" + to_string(Strike);
					//Positions of Call and Put
					int i1Call = rowContract[keyCall];
					int i2Call = colContract[keyCall];
					int i1Put = rowContract[keyPut];
					int i2Put = colContract[keyPut];
					option callOption = optionsPortfolio[i1Call][i2Call];
					option putOption = optionsPortfolio[i1Put][i2Put];
					//Both Call and Put
					callOption.underlyingPrice = FutPx / 100.0;
					putOption.underlyingPrice = FutPx / 100.0;
					callOption.maturity = DTE / 252.0;
					putOption.maturity = DTE / 252.0;
					callOption.rd = 0.0;
					putOption.rd = 0.0;
					callOption.rf = 0.0;
					putOption.rf = 0.0;
					callOption.impliedVol = Vol / 100.0;
					putOption.impliedVol = Vol / 100.0;
					double d1 = (log(callOption.underlyingPrice) - log(callOption.strike) + (callOption.rd - callOption.rf + 1.0 / 2.0 * pow(callOption.impliedVol, 2)) * callOption.maturity) / (callOption.impliedVol * sqrt(callOption.maturity));
					callOption.d1 = d1;
					putOption.d1 = d1;
					double gamma = 1.0 / sqrt(2.0 * M_PI) * exp(-pow(d1, 2) / 2.0) * exp(-callOption.rf * callOption.maturity) / (callOption.underlyingPrice * callOption.impliedVol * sqrt(callOption.maturity));
					callOption.gamma = gamma;
					putOption.gamma = gamma;
					//Call Specific Information
					optionsPortfolio[i1Call][i2Call].price = Call / 100.0;
					optionsPortfolio[i1Call][i2Call].delta = BSCDlta / 100.0;
					optionsPortfolio[i1Call][i2Call].underlyingPrice = callOption.underlyingPrice;
					optionsPortfolio[i1Call][i2Call].maturity = callOption.maturity;
					optionsPortfolio[i1Call][i2Call].rd = callOption.rd;
					optionsPortfolio[i1Call][i2Call].rf = callOption.rf;
					optionsPortfolio[i1Call][i2Call].impliedVol = callOption.impliedVol;
					optionsPortfolio[i1Call][i2Call].gamma = callOption.gamma;
					//Put Specific Information
					optionsPortfolio[i1Put][i2Put].price = Put / 100.0;
					optionsPortfolio[i1Put][i2Put].delta = BSPDlta / 100.0;
					optionsPortfolio[i1Put][i2Put].underlyingPrice = putOption.underlyingPrice;
					optionsPortfolio[i1Put][i2Put].maturity = putOption.maturity;
					optionsPortfolio[i1Put][i2Put].rd = putOption.rd;
					optionsPortfolio[i1Put][i2Put].rf = putOption.rf;
					optionsPortfolio[i1Put][i2Put].impliedVol = putOption.impliedVol;
					optionsPortfolio[i1Put][i2Put].gamma = putOption.gamma;
				}

	}
	return optionsPortfolio;
	}

CovCal getHistoricalData() {
	string inputFileName = "TestPrices.csv";
	CovCal historicalData(inputFileName);
	//Output contract title
	int nCol = historicalData.nCol;
	ofstream contractfile("title.csv");
	for (int i = 0; i < nCol; i++) {
		contractfile << historicalData.title[i] << ",";
	}
	contractfile.close();
	//EM fixes data and outputs cov and mean estimates in csv files
	EMfixData(historicalData.returns, historicalData.returnMissing, historicalData.nRow - 1, historicalData.nCol, historicalData.fixedData, historicalData.estiCovMat, historicalData.estiMean,historicalData.title,"mean.csv","cov.csv");
	return historicalData;
}

double** getCovMat(int& col){
	ifstream covfile("cov.csv");
	//Get covariance matrix from covfile
	int nCol = 0;
	string line;
	while (getline(covfile, line)) {
		nCol++;
	}
	//Looping through the file
	covfile.clear();
	covfile.seekg(0);
	double** covMat = new double*[nCol];
	for (int i = 0; i < nCol; i++) {
		covMat[i] = new double[nCol];
	}
	int i = 0;
	int j = 0;
	string value;
	while (getline(covfile, line)) {
		stringstream lineToProcess(line);
		while (getline(lineToProcess, value, ',')) {
			covMat[i][j] = atof(value.c_str());
			j++;
		}
		i++;
		j = 0;
	}
	covfile.close();
	col = nCol;
	return covMat;
}

double* getMean() {
	ifstream meanfile("mean.csv");
	//Get mean vector from meanfile
	int nCol = 0;
	string line;
	string value;
	getline(meanfile, line);
	stringstream lineToProcess(line);
	while (getline(lineToProcess, value, ',')) {
		nCol++;
	}
	meanfile.clear();
	meanfile.seekg(0);
	double* meanPtr = new double[nCol];
	int j = 0;
	getline(meanfile, line);
	stringstream row(line);
	while (getline(row,value, ',')){
			meanPtr[j] = atof(value.c_str());
			j++;
	}
	meanfile.close();
	return meanPtr;
}

string* getTitle() {
	ifstream contractfile("title.csv");
	//Get mean vector from meanfile
	int nCol = 0;
	string line;
	string value;
	getline(contractfile, line);
	stringstream lineToProcess(line);
	while (getline(lineToProcess, value, ',')) {
		nCol++;
	}
	contractfile.clear();
	contractfile.seekg(0);
	string* title= new string[nCol];
	getline(contractfile, line);
	stringstream row(line);
	int j = 0;
	while (getline(row, value, ',')) {
		title[j] = value;
		j++;
	}
	contractfile.close();
	return title;
}