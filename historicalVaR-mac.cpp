
/*
 Author: Yiming Li, yl3573
 Start Date: Oct, 2017
 Discription: This .cpp conducts the historical VaR calculation for the given portfolio
 */

#include "varCal.h"
#include "functions.h"
#include "csv.h"

struct OptionsForCommod{
public:
    string underlying;
    double underlyingPrice;
    double MCunderlyingPrice;
    double position;
    double maturity;
    double rd;  //domestic risk free interest rate
    double rf;  //foreign risk free interest rate or divident
    double strike;
    double impliedVol;
    double impliedVol_LP;
    string optionType;
    double optionPrice;
    double priceEvaluatedByBSformula;
    
};

struct Future{
public:
    string title;
    double position;
    double price;
    double priceChange;
    
};

int main1(){
    
    double quantile = 0.05;
    //The quantile we use for historical VaR calculation
    
    cout << "=============================================" << endl;
    cout << "               data processing " << endl;
    cout << "=============================================" << endl;
    
    string path = "/Users/FZhou/Desktop/Yiming/VAR/";
    
    //io::CSVReader <4> futureInformation("Sample CL Futures Position.csv");
    io::CSVReader <4> futureInformation(path +"Sample CL Futures Position.csv");
    futureInformation.read_header(io::ignore_extra_column, "Trader", "Position", "FutPx", "Sett");
    
    string instrument;
    double position;
    double price;
    double Sett;
    
    vector <Future> futureInPortfolio = vector <Future> ();
    // futureInPortfolio: is a vector of Future struct, see the definition above
    // Future has member title, position, price
    
    string sym;
    while(futureInformation.read_row(instrument, position, price, Sett)){
        
        Future futures;
        if(instrument == "Total:") continue;
        if(string::npos == instrument.find_first_of("0123456789")){sym = instrument;}
        else{
            futures.title = sym + instrument;
            futures.position = position;
            futures.price = price / 100.0;
            futureInPortfolio.push_back(futures);
        }
    }
    /*
     cout << futureInPortfolio.size() << endl;
     
     double portfolioPrice1=0.0;
     for(int i1=0; i1 < futureInPortfolio.size(); ++i1){
     portfolioPrice1 += futureInPortfolio[i1].position * futureInPortfolio[i1].price;
     
     }
     cout << portfolioPrice1 << endl;
    */
    /*
    for (unsigned i=0; i<futureInPortfolio.size(); i++){
        cout<<i<<endl;
        cout << ' ' << futureInPortfolio.at(i).title<<' ' << futureInPortfolio.at(i).position<<
            ' ' << futureInPortfolio.at(i).price;
        cout << '\n';
    }
    */
    
    
    io::CSVReader<21> optionPosition(path +"Sample LO Options Position.csv");
    
    optionPosition.read_header(io::ignore_extra_column, "Trader", "F17C", "F17P", "G17C", "G17P", "H17C", "H17P", "J17C", "J17P", "K17C", "K17P", "M17C", "M17P", "N17C", "N17P", "Q17C", "Q17P", "U17C", "U17P", "Z17C", "Z17P");
    
    string optionTitle[] = {"CLF17", "CLF17", "CLG17", "CLG17", "CLH17", "CLH17", "CLJ17", "CLJ17", "CLK17", "CLK17", "CLM17", "CLM17", "CLN17", "CLN17", "CLQ17", "CLQ17", "CLU17", "CLU17", "CLZ17", "CLZ17"};
    
    vector < vector <OptionsForCommod> > optionInPortfolio;
    //optionInPortfolio is a vector of vector of OptionsForCommod
    //The inner vector, vector<OptionsForCommod> has a same Strike
    //outer vector is of size equal to the # of strikes, inner vector is of size 20 or # of kinds of options
    
    int optionRow = 0;
    string firstColumn;
    vector<string> underlyingAsset = vector <string> (20);
    while (optionPosition.read_row(firstColumn, underlyingAsset[0], underlyingAsset[1], underlyingAsset[2], underlyingAsset[3], underlyingAsset[4], underlyingAsset[5], underlyingAsset[6], underlyingAsset[7], underlyingAsset[8], underlyingAsset[9], underlyingAsset[10], underlyingAsset[11], underlyingAsset[12], underlyingAsset[13], underlyingAsset[14], underlyingAsset[15], underlyingAsset[16], underlyingAsset[17], underlyingAsset[18], underlyingAsset[19])){
        if(firstColumn == "CL"){
            continue;
        }
        else{
            vector <OptionsForCommod> options = vector <OptionsForCommod> (20);
            for(int i1=0; i1 < 20; ++i1){
                options[i1].underlying = optionTitle[i1];
                if(i1 % 2 == 0){
                    options[i1].optionType = "call";
                }
                else{
                    options[i1].optionType = "put";
                }
                options[i1].strike = atof(firstColumn.c_str()) / 100.0;
                options[i1].position = atof((underlyingAsset[i1]).c_str());
            }
            optionInPortfolio.push_back(options);
            optionRow += 1;
        }
    }
    /*
     for(int i1=0; i1 < optionRow; ++i1){
         for(int i2=0; i2 < 20; ++i2){
                cout << optionInPortfolio[i1][i2].strike<<" , ";
                cout << optionInPortfolio[i1][i2].position << " , ";
         }
         cout << " " << endl<<endl;
     }
     */
    /*
    cout<<optionInPortfolio[0][3].optionType<<endl;
    cout<<optionInPortfolio[0][3].strike<<endl;
    cout<<optionInPortfolio[0][3].underlying<<endl;
    cout<<optionInPortfolio[0][3].position<<endl;
    */
    io::CSVReader<10> optionInformation(path +"Vol Smile.csv");
    
    optionInformation.read_header(io::ignore_extra_column, "Opt", "Cntr", "Stk", "Call", "Put", "Vol", "FutPx", "IntRate", "DTE", "DividendYield");
    
    string Opt;
    string Cntr;
    double Stk;
    double Call;
    double Put;
    double Vol;
    double FutPx;
    string IntRate;
    double DTE;
    string DividendYield;
    
    map<string, map<double, double> > impliedVol_grid;
    //impliedVol_grid is a map, key is a string and value is another map
    // for map in value, the key is a double and the value is also a double
    
    while(optionInformation.read_row(Opt, Cntr, Stk, Call, Put, Vol, FutPx, IntRate, DTE, DividendYield)){
        impliedVol_grid[Opt + Cntr][(Stk - FutPx) / 100.0] = Vol / 100.0;
        for(int i1=0; i1 < optionRow; ++i1){
            for(int i2=0; i2 < 20; ++i2){
                if((Opt + Cntr == optionInPortfolio[i1][i2].underlying) && (Stk / 100.0 == optionInPortfolio[i1][i2].strike)){
                    optionInPortfolio[i1][i2].underlyingPrice = FutPx / 100.0;
                    optionInPortfolio[i1][i2].maturity = DTE / 252.0;
                    optionInPortfolio[i1][i2].rd = 0;// This should be risk free interest rate which I use three month US treasury interest rate on 12/8/2016 instead
                    //optionInPortfolio[i1][i2].rf = 0.005;
                    optionInPortfolio[i1][i2].rf = atof(DividendYield.c_str()) / 100.0;
                    optionInPortfolio[i1][i2].impliedVol = Vol / 100.0;
                    if(optionInPortfolio[i1][i2].optionType == "call"){
                        optionInPortfolio[i1][i2].optionPrice = Call / 100.0;
                    }
                    else{
                        optionInPortfolio[i1][i2].optionPrice = Put / 100.0;
                    }
                }
            }
        }
    }
    //cout<<impliedVol_grid["CLF17"][-40.84]<<endl;
    // Will show 1.1229, the first element in the impliedVol_grid
    
    /*
     for(map<string, map<double, double> >::iterator iter1 = impliedVol_grid.begin(); iter1 != impliedVol_grid.end(); ++iter1){
         for(map<double, double>::iterator iter2 = iter1->second.begin(); iter2 != iter1->second.end(); ++iter2){
                cout << iter1->first << " " << iter2->first << " " << iter2->second << endl;
         }
     }
     */
    /*
     double portfolioPrice2=0.0;
     
     for(int i1=0; i1 < optionInPortfolio.size(); ++i1){
         for(int i2=0; i2 < optionInPortfolio[i1].size(); ++i2){
                portfolioPrice2 += optionInPortfolio[i1][i2].position * optionInPortfolio[i1][i2].optionPrice;
                cout << optionInPortfolio[i1][i2].optionPrice << endl;
         }
     }
     cout << portfolioPrice2 << endl;
     */
    
    cout << "=============================================" << endl;
    cout << "                EM algorithm" << endl;
    cout << "=============================================" << endl;
    
    string inputFileName = path +"TestPrices.csv";
    CovCal historicalData(inputFileName);
    //historicalData.getData();
    EMfixData(historicalData.myData, historicalData.ifExist, historicalData.nRow, historicalData.nCol, historicalData.fixedData, historicalData.estiCovMat, historicalData.estiMean);
    //historicalData.getMeanandVariance();
    //historicalData.getFixedData();
        
        
        
    //save the historical return matrix in to csv file, no need for running EMalgo again
        
    string outputFileName_HistoricalMeanBefore = "HistoricalMeanBefore.csv";
    ofstream output_HistoricalMeanBefore;
    output_HistoricalMeanBefore.open(path + outputFileName_HistoricalMeanBefore);
    output_HistoricalMeanBefore << "HistorialMeanBefore,";
    for(int i1=0; i1 < historicalData.nCol; ++i1){
            output_HistoricalMeanBefore << historicalData.title[i1] << ",";
    }
    output_HistoricalMeanBefore << 0 << endl;
    for(int i1=0 ;i1 < historicalData.nRow; ++i1){
        output_HistoricalMeanBefore << historicalData.date[i1] << ",";
        for(int i2=0; i2 < historicalData.nCol; ++i2){
            output_HistoricalMeanBefore << to_string(historicalData.myData[i1][i2]) << ",";
        }
        output_HistoricalMeanBefore << 0 << endl;
    }
    output_HistoricalMeanBefore.close();
    
    
    vector<vector<double>> historicalReturn;
    
    string outputFileName_HistoricalMeanAfter = "HistoricalMeanAfter.csv";
    ofstream output_HistoricalMeanAfter;
    output_HistoricalMeanAfter.open(path + outputFileName_HistoricalMeanAfter);
    output_HistoricalMeanAfter << "HistorialMeanAfter,";
    for(int i1=0; i1 < historicalData.nCol; ++i1){
        output_HistoricalMeanAfter << historicalData.title[i1] << ",";
    }
    output_HistoricalMeanAfter << 0 << endl;
    for(int i1=0 ;i1 < historicalData.nRow; ++i1){
        output_HistoricalMeanAfter << historicalData.date[i1] << ",";
        vector <double> dailyReturn = vector <double> (historicalData.nCol);
        for(int i2=0; i2 < historicalData.nCol; ++i2){
            output_HistoricalMeanAfter << to_string(historicalData.fixedData[i1][i2]) << ",";
            dailyReturn[i2] =historicalData.fixedData[i1][i2];
        }
        output_HistoricalMeanAfter << 0 << endl;
        historicalReturn.push_back(dailyReturn);
    }
    output_HistoricalMeanAfter.close();

    /*
    for(int i1=0 ;i1 < historicalData.nRow; ++i1){
        for(int i2=0; i2 < historicalData.nCol; ++i2){
            cout<<historicalReturn[i1][i2]<<",";
        }
        cout<<endl;
    }
    */
        
        

        
    //Save the best estimated covariance matrix
    string outputFileName_Covariance = "covarianceMatrix.csv";
    ofstream output_Covariance;
    output_Covariance.open(path + outputFileName_Covariance);
    output_Covariance << "Cov,";
    for(int i1=0; i1 < historicalData.nCol; ++i1){
        output_Covariance << historicalData.title[i1] << ",";
    }
    output_Covariance << 0 << endl;
    for(int i1=0 ;i1 < historicalData.nCol; ++i1){
        output_Covariance << historicalData.title[i1] << ",";
        for(int i2=0; i2 < historicalData.nCol; ++i2){
            output_Covariance << to_string(historicalData.estiCovMat[i1][i2]) << ",";
        }
        output_Covariance << 0 << endl;
    }
    output_Covariance.close();

    
    
    cout << "=============================================" << endl;
    cout << "               VaR Calculation" << endl;
    cout << "=============================================" << endl;
    /*
    vector<vector<double>> historicalReturn1(historicalReturn[0].size(),
                                             vector<double>(historicalReturn.size()));
    for (int i=0;i<historicalReturn.size();++i){
        for (int j=0; j<historicalReturn[0].size();++j){
            historicalReturn1[j][i] = historicalReturn[i][j];
        }
    }
    
    cout<<historicalReturn1.size()<<" and "<<historicalReturn1[0].size()<<endl;
    
    for(int i=0;i<historicalReturn1.size();++i){
        std::sort(historicalReturn1[i].begin(), historicalReturn1[i].end());
    }
    
    
    
    for(int i1=0 ;i1 <historicalReturn1.size(); ++i1){
        for(int i2=0; i2 < historicalReturn1[0].size(); ++i2){
            cout<<historicalReturn[i1][i2]<<",";
        }
        cout << "=============================================" << endl;
        for(int i2=0; i2 < historicalReturn1[0].size(); ++i2){
            cout<<historicalReturn1[i1][i2]<<",";
        }
        cout<<endl;
    }
    
    int positionOfQuantile = int(quantile * historicalReturn1[0].size());
    cout<<endl<<positionOfQuantile<<endl;
    */
    //cout<<historicalData.nRow<<endl;  264
    //cout<<historicalData.nCol<<endl;  40
    //cout<<futureInPortfolio.size()<<endl;  40
    
    for(int i1=0 ;i1 < historicalData.nRow; ++i1){
        for(int i2=0; i2 < historicalData.nCol; ++i2){
            //historicalData.fixedData[i1][i2];
        }
    }
    
    int numOfVar = historicalData.nRow;
    // 264, equal to the number of set of return we have
    
    
    double *simulatedPortPriceDif;
    simulatedPortPriceDif = new double [numOfVar];
    for(int i1=0; i1 < numOfVar; ++i1){
        simulatedPortPriceDif[i1] = 0.0;
    }
    // save the final Var for each date in history
    
    string historicalVarFileName;
    historicalVarFileName = path +"HistoricalMeanAfter.csv";

    
    CovCal historicalReturnData(historicalVarFileName);
    double lower_Vol_loc; // linear interpolation
    double lower_Vol;
    double upper_Vol_loc;
    double upper_Vol = 0.0;
    for(int i1=0; i1 < numOfVar; ++i1){
        for(int i2=0; i2 < futureInPortfolio.size(); ++i2){
            for(int i3=0; i3 < historicalReturnData.nCol; ++i3){
                if(historicalReturnData.title[i3] == futureInPortfolio[i2].title){
                    /*
                     && (futureInPortfolio[i2].title != "CLF17") &&(futureInPortfolio[i2].title != "CLG17") && (futureInPortfolio[i2].title != "CLH17") && (futureInPortfolio[i2].title != "CLJ17") && (futureInPortfolio[i2].title != "CLK17") && (futureInPortfolio[i2].title != "CLM17") &&(futureInPortfolio[i2].title != "CLN17") && (futureInPortfolio[i2].title != "CLQ17") && (futureInPortfolio[i2].title != "CLU17") && (futureInPortfolio[i2].title != "CLZ17")
                     */
                    futureInPortfolio[i2].priceChange = futureInPortfolio[i2].price * historicalReturnData.myOriginalData[i1][i3];
                    //cout<<"the futureInPortfolio "<<i2<<" 's"<<futureInPortfolio[i2].title<<"'s price change is "<<futureInPortfolio[i2].priceChange<<endl;
                }
            }
            simulatedPortPriceDif[i1] += futureInPortfolio[i2].priceChange * futureInPortfolio[i2].position;
        }
        
        for(int i2=0; i2 < optionInPortfolio.size(); ++i2){
            for(int i3=0; i3 < optionInPortfolio[i2].size(); ++i3){
                for(int i4=0; i4 < historicalReturnData.nCol; ++i4){
                    if((optionInPortfolio[i2][i3].underlying == historicalReturnData.title[i4]) && (optionInPortfolio[i2][i3].position != 0.0)){
                        optionInPortfolio[i2][i3].MCunderlyingPrice = optionInPortfolio[i2][i3].underlyingPrice * (1.0 + historicalReturnData.myOriginalData[i1][i4]);
                        lower_Vol = 0.0;
                        upper_Vol = 0.0;
                        for(map<string, map<double, double> >::iterator iter1 = impliedVol_grid.begin(); iter1 != impliedVol_grid.end(); ++iter1){
                            if(iter1->first == optionInPortfolio[i2][i3].underlying){
                                lower_Vol = iter1 -> second.begin() -> second;
                                for(map<double, double>::iterator iter2 = iter1->second.begin(); iter2 != iter1->second.end(); ++iter2){
                                    
                                    if(iter2->first <= (optionInPortfolio[i2][i3].strike - optionInPortfolio[i2][i3].MCunderlyingPrice)){
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
                            optionInPortfolio[i2][i3].impliedVol_LP = lower_Vol;
                        }
                        else{
                            if(upper_Vol == 0.0){
                                optionInPortfolio[i2][i3].impliedVol_LP = lower_Vol;
                            }
                            else{
                                optionInPortfolio[i2][i3].impliedVol_LP = lower_Vol * (upper_Vol_loc - (optionInPortfolio[i2][i3].strike - optionInPortfolio[i2][i3].MCunderlyingPrice)) / (upper_Vol_loc -lower_Vol_loc) + upper_Vol * (optionInPortfolio[i2][i3].strike - optionInPortfolio[i2][i3].MCunderlyingPrice - lower_Vol_loc) / (upper_Vol_loc - lower_Vol_loc);
                            }
                        }
                        optionInPortfolio[i2][i3].priceEvaluatedByBSformula = calPriceForForexOption(optionInPortfolio[i2][i3].MCunderlyingPrice, optionInPortfolio[i2][i3].strike, optionInPortfolio[i2][i3].rd, optionInPortfolio[i2][i3].rf, optionInPortfolio[i2][i3].maturity - (1.0 / 252.0), optionInPortfolio[i2][i3].impliedVol_LP, optionInPortfolio[i2][i3].optionType);
                        if(isnan(optionInPortfolio[i2][i3].priceEvaluatedByBSformula - optionInPortfolio[i2][i3].optionPrice)){
                            simulatedPortPriceDif[i1] += 0;
                        }
                        else{
                            simulatedPortPriceDif[i1] += optionInPortfolio[i2][i3].position * (optionInPortfolio[i2][i3].priceEvaluatedByBSformula - optionInPortfolio[i2][i3].optionPrice);
                        }
                    }
                }
            }
        }
    }
    
    // save the option information
    string outputFileName_infoOnOptions = "optionInfo.csv";
    ofstream output_infoOnOptions;
    output_infoOnOptions.open(path+outputFileName_infoOnOptions);
    
    output_infoOnOptions << "underlying,";
    
    for(int i1=0; i1 < optionInPortfolio.size(); ++i1){
        for(int i2=0; i2 < optionInPortfolio[i1].size(); ++i2){
            output_infoOnOptions << optionInPortfolio[i1][i2].underlying << ",";
        }
    }
    
    output_infoOnOptions << 0 << endl;
    
    output_infoOnOptions << "underlying price,";
    
    for(int i1=0; i1 < optionInPortfolio.size(); ++i1){
        for(int i2=0; i2 < optionInPortfolio[i1].size(); ++i2){
            output_infoOnOptions << optionInPortfolio[i1][i2].underlyingPrice << ",";
        }
    }
    
    output_infoOnOptions << 0 << endl;
    /*
     output_infoOnOptions << "MC underlying price,";
     
     for(int i1=0; i1 < optionInPortfolio.size(); ++i1){
     for(int i2=0; i2 < optionInPortfolio[i1].size(); ++i2){
     output_infoOnOptions << optionInPortfolio[i1][i2].MCunderlyingPrice << ",";
     }
     }
     
     output_infoOnOptions << 0 << endl;
     */
    output_infoOnOptions << "position,";
    
    for(int i1=0; i1 <optionInPortfolio.size(); ++i1){
        for(int i2=0; i2 < optionInPortfolio[i1].size(); ++i2){
            output_infoOnOptions << to_string(optionInPortfolio[i1][i2].position) << ",";
        }
    }
    
    output_infoOnOptions << 0 << endl;
    
    output_infoOnOptions << "maturity,";
    
    for(int i1=0; i1 < optionInPortfolio.size(); ++i1){
        for(int i2=0; i2 < optionInPortfolio[i1].size(); ++i2){
            output_infoOnOptions << to_string(optionInPortfolio[i1][i2].maturity) << ",";
        }
    }
    
    output_infoOnOptions << 0 << endl;
    
    output_infoOnOptions << "divident yield,";
    
    for(int i1=0; i1 < optionInPortfolio.size(); ++i1){
        for(int i2=0; i2 < optionInPortfolio[i1].size(); ++i2){
            output_infoOnOptions << to_string(optionInPortfolio[i1][i2].rf) << ",";
        }
    }
    
    output_infoOnOptions << 0 << endl;
    
    output_infoOnOptions << "strike,";
    
    for(int i1=0; i1 < optionInPortfolio.size(); ++i1){
        for(int i2=0; i2 < optionInPortfolio[i1].size(); ++i2){
            output_infoOnOptions << to_string(optionInPortfolio[i1][i2].strike) << ",";
        }
    }
    
    output_infoOnOptions << 0 << endl;
    
    output_infoOnOptions << "implied Vol,";
    
    for(int i1=0; i1 < optionInPortfolio.size(); ++i1){
        for(int i2=0; i2 < optionInPortfolio[i1].size(); ++i2){
            output_infoOnOptions << to_string(optionInPortfolio[i1][i2].impliedVol) << ",";
        }
    }
    
    output_infoOnOptions << 0 << endl;
    /*
     output_infoOnOptions << "implied Vol PL,";
     
     for(int i1=0; i1 < optionInPortfolio.size(); ++i1){
     for(int i2=0; i2 < optionInPortfolio[i1].size(); ++i2){
     output_infoOnOptions << to_string(optionInPortfolio[i1][i2].impliedVol_LP) << ",";
     }
     }
     
     output_infoOnOptions << 0 << endl;
     */
    output_infoOnOptions << "option Type,";
    
    for(int i1=0; i1 < optionInPortfolio.size(); ++i1){
        for(int i2=0; i2 < optionInPortfolio[i1].size(); ++i2){
            output_infoOnOptions << optionInPortfolio[i1][i2].optionType << ",";
        }
    }
    
    output_infoOnOptions << 0 << endl;
    
    output_infoOnOptions << "option Price,";
    
    for(int i1=0; i1 < optionInPortfolio.size(); ++i1){
        for(int i2=0; i2 < optionInPortfolio[i1].size(); ++i2){
            output_infoOnOptions << to_string(optionInPortfolio[i1][i2].optionPrice) << ",";
        }
    }
    output_infoOnOptions << 0 << endl;
    
    output_infoOnOptions.close();
    
    double media;
    for(int i1=0; i1 < numOfVar; ++i1){
        for(int i2=i1+1; i2 < numOfVar; ++i2){
            if(simulatedPortPriceDif[i1] > simulatedPortPriceDif[i2]){
                media = simulatedPortPriceDif[i1];
                simulatedPortPriceDif[i1] = simulatedPortPriceDif[i2];
                simulatedPortPriceDif[i2] = media;
            }
        }
    }
    
    string outputFileName_PL = path+"unrealistic P&L.csv";
    ofstream output_PL(outputFileName_PL);
    output_PL << "MC simulation," << "P&L of the portfolio," << 0 << endl;
    for(int i1=0; i1 < numOfVar; ++i1){
        output_PL << "path" + to_string(i1 + 1) << "," << simulatedPortPriceDif[i1] << "," << 0 << endl;
    }
    output_PL.close();
    
    double portfolioPrice = 0.0;
    
    for(int i1=0; i1 < futureInPortfolio.size(); ++i1){
        portfolioPrice += futureInPortfolio[i1].price * futureInPortfolio[i1].position;
    }
    
    for(int i1=0; i1 < optionInPortfolio.size(); ++i1){
        for(int i2=0; i2 < optionInPortfolio[i1].size(); ++i2){
            portfolioPrice += optionInPortfolio[i1][i2].optionPrice * optionInPortfolio[i1][i2].position;
        }
    }
    
    cout << "The value of the portfolio is: " << portfolioPrice << endl;
    cout << "The VaR of the portfolio is: " << -simulatedPortPriceDif[int(0.01 * numOfVar) -1] << endl;
    //0.01 is the VaR, significant level
    delete simulatedPortPriceDif;

    

    
    
    
    cout << "*********************************************" << endl;
    cout<<"             We are at the end now!            "<<endl;
    cout << "*********************************************" << endl;
    
    return 0;

}