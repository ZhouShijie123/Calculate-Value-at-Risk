#include "stdafx.h"
#include "varCal.h"
#include "csv.h"
#include "matrixops.h"
#include "normalDist.h"
#include "EMtry.h"
#include "varMethods.h"
#include "getPortfolio.h"
#include "BlackScholes.h"
#include "normalDist.h"
#include <iostream>
#include <fstream>


void MonteCarloVaR(vector<Future> futuresPortfolio, vector<vector<option>> optionsPortfolio, map<string, map<double, double> > impliedVol_grid, 
	double** covMat, string* title, int nCol) {

    //Monte Carlo simulation
    double *riskFreeRate;
    riskFreeRate = new double [nCol];
    for(int i1=0; i1 < nCol; ++i1){
        //riskFreeRate[i1] = 0.005 / 252.0;
        riskFreeRate[i1] = 0.0;
    }
    
    int numOfPath;
    numOfPath = 1000;
    string simulationFileName;
    simulationFileName = "simulationData.csv";
    simulationOfMultinormal(riskFreeRate, covMat, simulationFileName, nCol, numOfPath, title);
    
    double *simulatedPortPriceDif;
	double *simfuturediff = new double[numOfPath];
    simulatedPortPriceDif = new double [numOfPath];
    for(int i1=0; i1 < numOfPath; ++i1){
        simulatedPortPriceDif[i1] = 0.0;
		simfuturediff[i1] = 0.0;
    }
    
    CovCal simulationData(simulationFileName);
    double lower_Vol_loc; // linear interpolation
    double lower_Vol;
    double upper_Vol_loc;
    double upper_Vol = 0.0;
    for(int i1=0; i1 < numOfPath; ++i1){
        for(int i2=0; i2 < futuresPortfolio.size(); ++i2){
            for(int i3=0; i3 < simulationData.nCol; ++i3){
                if(simulationData.title[i3] == futuresPortfolio[i2].title){
                    /*
&& (futuresPortfolio[i2].title != "CLF17") &&(futuresPortfolio[i2].title != "CLG17") && (futuresPortfolio[i2].title != "CLH17") && (futuresPortfolio[i2].title != "CLJ17") && (futuresPortfolio[i2].title != "CLK17") && (futuresPortfolio[i2].title != "CLM17") &&(futuresPortfolio[i2].title != "CLN17") && (futuresPortfolio[i2].title != "CLQ17") && (futuresPortfolio[i2].title != "CLU17") && (futuresPortfolio[i2].title != "CLZ17")
                                                                               */
                    futuresPortfolio[i2].priceChange = futuresPortfolio[i2].price * simulationData.OriginalData[i1][i3];
                }
            }
            simulatedPortPriceDif[i1] += futuresPortfolio[i2].priceChange * futuresPortfolio[i2].position;
			simfuturediff[i1] += futuresPortfolio[i2].priceChange * futuresPortfolio[i2].position;
        }

        for(int i2=0; i2 < optionsPortfolio.size(); ++i2){
            for(int i3=0; i3 < optionsPortfolio[i2].size(); ++i3){
                for(int i4=0; i4 < simulationData.nCol; ++i4){
                    if((optionsPortfolio[i2][i3].underlying == simulationData.title[i4]) && (optionsPortfolio[i2][i3].position != 0.0)){
                        optionsPortfolio[i2][i3].MCunderlyingPrice = optionsPortfolio[i2][i3].underlyingPrice * (1.0 + simulationData.OriginalData[i1][i4]);
                        lower_Vol = 0.0;
                        upper_Vol = 0.0;
                        for(map<string, map<double, double> >::iterator iter1 = impliedVol_grid.begin(); iter1 != impliedVol_grid.end(); ++iter1){
                            if(iter1->first == optionsPortfolio[i2][i3].underlying){
                                lower_Vol = iter1 -> second.begin() -> second;
                                for(map<double, double>::iterator iter2 = iter1->second.begin(); iter2 != iter1->second.end(); ++iter2){
    
                                    if(iter2->first <= (optionsPortfolio[i2][i3].strike - optionsPortfolio[i2][i3].MCunderlyingPrice)){
                                        lower_Vol_loc = iter2->first;
                                        lower_Vol = iter2->second;
                                    }
                                    else{
                                        upper_Vol_loc = iter2->first;
                                        upper_Vol = iter2->second;
                                        break;
                                    }
                                }
                            }
                        }
                        if(lower_Vol == upper_Vol){
                            optionsPortfolio[i2][i3].impliedVol_LP = lower_Vol;
                        }
                        else{
                            if(upper_Vol == 0.0){
                                optionsPortfolio[i2][i3].impliedVol_LP = lower_Vol;
                            }
                            else{
                                optionsPortfolio[i2][i3].impliedVol_LP = lower_Vol * (upper_Vol_loc - (optionsPortfolio[i2][i3].strike - optionsPortfolio[i2][i3].MCunderlyingPrice)) / (upper_Vol_loc - lower_Vol_loc) + upper_Vol * (optionsPortfolio[i2][i3].strike - optionsPortfolio[i2][i3].MCunderlyingPrice - lower_Vol_loc) / (upper_Vol_loc - lower_Vol_loc);
								
							}
                        }
                        optionsPortfolio[i2][i3].priceEvaluatedByBSformula = optionPriceBS(optionsPortfolio[i2][i3].MCunderlyingPrice, optionsPortfolio[i2][i3].strike, optionsPortfolio[i2][i3].rd, optionsPortfolio[i2][i3].rf, optionsPortfolio[i2][i3].maturity - (1.0 / 252.0), optionsPortfolio[i2][i3].impliedVol_LP, optionsPortfolio[i2][i3].optionType);
                        if(isnan(optionsPortfolio[i2][i3].priceEvaluatedByBSformula - optionsPortfolio[i2][i3].price)){
                            simulatedPortPriceDif[i1] += 0;
                        }
                        else{
                            simulatedPortPriceDif[i1] += optionsPortfolio[i2][i3].position * (optionsPortfolio[i2][i3].priceEvaluatedByBSformula - optionsPortfolio[i2][i3].price);
                        }
                    }
                }
            }
        }
		for (int i = 0; i < optionsPortfolio[0].size(); i += 2) {

			for (int j = 0; j < optionsPortfolio.size(); j++) {
				if (i == 0) {
					optionsPortfolio[j][i].impliedVol_LP *= 5 / 6;
					optionsPortfolio[j][i + i].impliedVol_LP *= 5 / 6;

				}
				else {
					optionsPortfolio[j][i].impliedVol_LP -= (optionsPortfolio[j][i].impliedVol - optionsPortfolio[j][i - 2].impliedVol) / (optionsPortfolio[j][i].maturity - optionsPortfolio[j][i - 2].maturity);
					optionsPortfolio[j][i + 1].impliedVol_LP = optionsPortfolio[j][i].impliedVol_LP;
				}

			}
		}
    }
	
    // save the option information
    string outputFileName_infoOnOptions = "optionInfo_MonteCarloVaR.csv";
    ofstream output_infoOnOptions;
    output_infoOnOptions.open(outputFileName_infoOnOptions);
    
    output_infoOnOptions << "underlying,";
    
    for(int i1=0; i1 < optionsPortfolio.size(); ++i1){
        for(int i2=0; i2 < optionsPortfolio[i1].size(); ++i2){
            output_infoOnOptions << optionsPortfolio[i1][i2].underlying << ",";
        }
    }
    
    output_infoOnOptions << 0 << endl;
    
    output_infoOnOptions << "underlying price,";
    
    for(int i1=0; i1 < optionsPortfolio.size(); ++i1){
        for(int i2=0; i2 < optionsPortfolio[i1].size(); ++i2){
            output_infoOnOptions << optionsPortfolio[i1][i2].underlyingPrice << ",";
        }
    }
    
    output_infoOnOptions << 0 << endl;
    /*
    output_infoOnOptions << "MC underlying price,";
    
    for(int i1=0; i1 < optionsPortfolio.size(); ++i1){
        for(int i2=0; i2 < optionsPortfolio[i1].size(); ++i2){
            output_infoOnOptions << optionsPortfolio[i1][i2].MCunderlyingPrice << ",";
        }
    }
    
    output_infoOnOptions << 0 << endl;
    */
    output_infoOnOptions << "position,";
    
    for(int i1=0; i1 <optionsPortfolio.size(); ++i1){
        for(int i2=0; i2 < optionsPortfolio[i1].size(); ++i2){
            output_infoOnOptions << to_string(optionsPortfolio[i1][i2].position) << ",";
        }
    }
    
    output_infoOnOptions << 0 << endl;
    
    output_infoOnOptions << "maturity,";
    
    for(int i1=0; i1 < optionsPortfolio.size(); ++i1){
        for(int i2=0; i2 < optionsPortfolio[i1].size(); ++i2){
            output_infoOnOptions << to_string(optionsPortfolio[i1][i2].maturity) << ",";
        }
    }
    
    output_infoOnOptions << 0 << endl;
    
    output_infoOnOptions << "divident yield,";
    
    for(int i1=0; i1 < optionsPortfolio.size(); ++i1){
        for(int i2=0; i2 < optionsPortfolio[i1].size(); ++i2){
            output_infoOnOptions << to_string(optionsPortfolio[i1][i2].rf) << ",";
        }
    }
    
    output_infoOnOptions << 0 << endl;
    
    output_infoOnOptions << "strike,";
    
    for(int i1=0; i1 < optionsPortfolio.size(); ++i1){
        for(int i2=0; i2 < optionsPortfolio[i1].size(); ++i2){
            output_infoOnOptions << to_string(optionsPortfolio[i1][i2].strike) << ",";
        }
    }
    
    output_infoOnOptions << 0 << endl;
    
    output_infoOnOptions << "implied Vol,";
    
    for(int i1=0; i1 < optionsPortfolio.size(); ++i1){
        for(int i2=0; i2 < optionsPortfolio[i1].size(); ++i2){
            output_infoOnOptions << to_string(optionsPortfolio[i1][i2].impliedVol) << ",";
        }
    }
    
    output_infoOnOptions << 0 << endl;
    /*
    output_infoOnOptions << "implied Vol PL,";
    
    for(int i1=0; i1 < optionsPortfolio.size(); ++i1){
        for(int i2=0; i2 < optionsPortfolio[i1].size(); ++i2){
            output_infoOnOptions << to_string(optionsPortfolio[i1][i2].impliedVol_LP) << ",";
        }
    }
    
    output_infoOnOptions << 0 << endl;
    */
    output_infoOnOptions << "option Type,";
    
    for(int i1=0; i1 < optionsPortfolio.size(); ++i1){
        for(int i2=0; i2 < optionsPortfolio[i1].size(); ++i2){
            output_infoOnOptions << optionsPortfolio[i1][i2].optionType << ",";
        }
    }
    
    output_infoOnOptions << 0 << endl;
    
    output_infoOnOptions << "option Price,";
    
    for(int i1=0; i1 < optionsPortfolio.size(); ++i1){
        for(int i2=0; i2 < optionsPortfolio[i1].size(); ++i2){
            output_infoOnOptions << to_string(optionsPortfolio[i1][i2].price) << ",";
        }
    }
    output_infoOnOptions << 0 << endl;
    
    output_infoOnOptions.close();
    
    double media;
    for(int i1=0; i1 < numOfPath; ++i1){
        for(int i2=i1+1; i2 < numOfPath; ++i2){
            if(simulatedPortPriceDif[i1] > simulatedPortPriceDif[i2]){
                media = simulatedPortPriceDif[i1];
                simulatedPortPriceDif[i1] = simulatedPortPriceDif[i2];
                simulatedPortPriceDif[i2] = media;
            }
        }
    }
    
    string outputFileName_PL = "unrealistic_P&L_MonteCarloVaR.csv";
    ofstream output_PL(outputFileName_PL);
    output_PL << "MC simulation," << "P&L of the portfolio," << 0 << endl;
    for(int i1=0; i1 < numOfPath; ++i1){
        output_PL << "path" + to_string(i1 + 1) << "," << simulatedPortPriceDif[i1] << "," << simfuturediff[i1] << endl;
    }
    output_PL.close();
    
    double portfolioPrice = 0.0;
    
    for(int i1=0; i1 < futuresPortfolio.size(); ++i1){
        portfolioPrice += futuresPortfolio[i1].price * futuresPortfolio[i1].position;
    }
    
    for(int i1=0; i1 < optionsPortfolio.size(); ++i1){
        for(int i2=0; i2 < optionsPortfolio[i1].size(); ++i2){
            portfolioPrice += optionsPortfolio[i1][i2].price * optionsPortfolio[i1][i2].position;
        }
    }
	for (int i = 0; i < numOfPath; i++)
		cout << simulatedPortPriceDif[i] << endl;
    cout << "The value of the portfolio is: " << portfolioPrice << endl;
    cout << "The VaR of the portfolio is: " << -simulatedPortPriceDif[int(0.01 * numOfPath) -1] << endl;
    
    delete riskFreeRate;
    delete simulatedPortPriceDif;
}

