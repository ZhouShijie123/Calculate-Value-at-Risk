#ifndef BlackScholes
#define BlackScholes
#include <string>

using namespace std;
double Phi(double mean, double vol, double location);
double optionPriceBS(double s0, double strike, double rd, double rf, double mature, double vol, string optionType);

#endif
