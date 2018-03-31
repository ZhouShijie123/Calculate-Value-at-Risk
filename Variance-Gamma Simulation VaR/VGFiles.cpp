#include "StdAfx.h"
#include "VGFiles.h"

CVGFiles::CVGFiles(void) { }
CVGFiles::~CVGFiles(void) { }

// N(0,1) density
double CVGFiles::normpdf(double x) {
	double pi = 3.141592653589793;
	return exp(-x*x*0.5)/sqrt(2*pi);
}

// Waissi and Rossin normal cdf approximation
double CVGFiles::normcdf(double z) {
	double pi = 3.141592653589793;
	double b1 = -0.0004406;
	double b2 =  0.0418198;
	double b3 =  0.9;
	if (z < -8)
		return 0.0;
	else if (z > 8)
		return 1.0;
	else
		return 1.0 / (1.0 + exp(-sqrt(pi)*(b1*pow(z,5.0) + b2*pow(z,3.0) + b3*z)));
}

// Variance Gamma call price
double CVGFiles::VGCall(double alpha,double sigma,double v,double S,double K,double r,double T) {
	double a = pow(alpha + sigma,2.0);
	double num = 1.0 - v*a/2.0;
	double den = 1.0 - v*alpha*alpha/2.0;
	double d1 = log(S/K)/sigma/T + ((r + 1.0/v*log(num/den))/sigma + alpha+sigma)*sqrt(T);
	double d2 = d1 - sigma*sqrt(T);
	return S*exp(a*T/2.0) * pow((1.0 - v*a/2.0),T/v)*normcdf(d1)
         - K*exp(-r*T + alpha*alpha*T/2.0) * pow((1.0-v*alpha*alpha/2.0),T/v)*normcdf(d2);
}

// Black Scholes call price
double CVGFiles::BSCall(double S,double K,double r,double q,double v,double T) {
	double d1 = (log(S/K) + T*(r - q + v*v/2.0))/v/sqrt(T);
	double d2 = d1 - v*sqrt(T);
	return S*exp(-q*T)*normcdf(d1) - exp(-r*T)*K*normcdf(d2);
}

// Random uniform number in  (a,b)
double CVGFiles::U(double a, double b) {
	std::random_device rd;
	std::default_random_engine e1(rd());
	std::uniform_real_distribution<double> uniform_dist(a,b);
	double U = uniform_dist(e1);
	return U;
}





