#ifndef VGFiles_H
#define VGFiles_H

#pragma once
#include <math.h>
#include <ctime>
#include <random>

class CVGFiles
{
public:
	CVGFiles(void);
	~CVGFiles(void);

	// N(0,1) density
	double CVGFiles::normpdf(double x);

	// Waissi and Rossin normal cdf approximation
	double CVGFiles::normcdf(double z);

	// Variance Gamma call price
	double CVGFiles::VGCall(double alpha,double sigma,double v,double S,double K,double r,double T);

	// Black Scholes call price
	double CVGFiles::BSCall(double S,double K,double r,double q,double v,double T);

	// Random uniform number in (a,b)
	double U(double a, double b);
		
};

#endif // VGFiles_H