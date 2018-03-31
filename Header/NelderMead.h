#ifndef NelderMead_H
#define NelderMead_H

#pragma once
#include <vector>
#include <algorithm>
#include "Structures.h"
using namespace std;

class CNelderMead
{
public:
	CNelderMead(void);
	~CNelderMead(void);

	// Function to find the standard deviation of a vector
	double CNelderMead::stdev(vector<double> x);

	// Function to calculate the mean value of
	// a set of n vectors each of dimension n
	// namely a (n x n) matrix
	vector<double> CNelderMead::VMean(vector<vector<double> > X, int n);

	// Function to add two vectors together
	vector<double> CNelderMead::VAdd(vector<double> x, vector<double> y);

	// Function to subtract two vectors
	vector<double> CNelderMead::VSub(vector<double> x, vector<double> y);

	// Function to multiply a vector by a constant
	vector<double> CNelderMead::VMult(vector<double> x, double a);

	// Nelder Mead Algorithm
	vector<double> CNelderMead::NelderMeadFunction(double (*f)(vector<double>, RMSEinputs rmsein), NMsettings nmset, vector<vector<double> > x);
};

#endif // NelderMead_H