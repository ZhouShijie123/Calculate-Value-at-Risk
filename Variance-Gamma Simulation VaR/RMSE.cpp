#include "StdAfx.h"
#include "RMSE.h"



CRMSE::CRMSE(void) { }
CRMSE::~CRMSE(void) { }

//double CRMSE::RMSE(vector<double> B,RMSEinputs rmsein) {
//	CVGFiles VG;
//	int M = rmsein.K.size();
//	vector<double> VGPrice(M,0.0);
//	double Sum = 0.0;
//	// Penalty for negative parameter values
//	if ( (B[1]<0.0) || (B[2]<0.0) )
//		return 1.0e100;
//    else {
//		// Create the loss function
//		for (int i=0; i<=M-1; i++)	{
//			VGPrice[i] =  VG.VGCall(B[0],B[1],B[2],rmsein.S,rmsein.K[i],rmsein.r,rmsein.T);
////			Sum +=  pow(log(rmsein.MPrice[i]) - log(VGPrice[i]), 2.0);
//			Sum +=  pow( (rmsein.MPrice[i] - VGPrice[i])/rmsein.MPrice[i], 2.0);
//		}
//	return sqrt(Sum/M);
//	}
//}

