#include "stdafx.h"
#include "BlackScholes.h"
#include <cmath>
#include <iostream>

using namespace std;

double Phi(double mean, double vol, double location) {
	//Cumulative distribution of a random normal variable
	double scaledLocation = (location - mean) / vol;

	double integral;
	integral = 0.5 + 1.0 / 2.0 * erf(scaledLocation / (sqrt(2.0)));
	return integral;
}

double optionPriceBS(double s0, double strike, double rd, double rf, double maturity, double vol, string optionType) {
	//Input:
	//rd domestic risk free interest rate
	//rf foreign risk free interest rate or dividend yield
	double d1;
	double d2;
	double priceOfOption;
	d1 = ((rd - rf + 1.0 / 2.0 * pow(vol, 2)) * maturity + log(s0) - log(strike)) / (vol * sqrt(maturity));
	d2 = ((rd - rf - 1.0 / 2.0 * pow(vol, 2)) * maturity + log(s0) - log(strike)) / (vol * sqrt(maturity));
	if (!optionType.compare("call")) {
		priceOfOption = s0 * exp(-rf * maturity) * Phi(0.0, 1.0, d1) - exp(-rd * maturity) * strike * Phi(0.0, 1.0, d2);
		return priceOfOption;
	}
	else {
		if (!optionType.compare("put")) {
			priceOfOption = exp(-rd * maturity) * strike * Phi(0.0, 1.0, -d2) - s0 * exp(-rf * maturity) * Phi(0.0, 1.0, -d1);
			return priceOfOption;
		}
		else {
			cout << "The last paramater should be \"put\" or \"call\"" << endl;
			return 0.0;
		}
	}
}

int test(){
cout << optionPriceBS(15147.35, 15147.35, 0.05, 0, 0.25, 0.0001, "call") << endl;
return 0;
}

