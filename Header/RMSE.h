#ifndef RMSE_H
#define RMSE_H

#pragma once
#include <vector>
#include "Structures.h"
#include "VGFiles.h"
using namespace std;

class CRMSE {
public:
	CRMSE(void);
	~CRMSE(void);

	// RMSE loss function for the Nelder Mead algorithm
	double CRMSE::RMSE(vector<double> B,RMSEinputs rmsein);

};

#endif // RMSE_H
