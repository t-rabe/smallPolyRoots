#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H
#define _USE_MATH_DEFINES

#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <cmath>
#include <thread>
#include <numeric>
#include <complex>

using namespace std;

/**
 *  Class that models a sequence of i.i.d. standard gaussians
 *  vector of distributions that can be sampled with random seeds, giving
 *  a model of uncorrelated Gaussian variables.
 */
class Polynomial{
	private:
        int arrSize;
        vector<double> realComp; // the real component of each element of the polynomial
        vector<double> imgComp; // the imaginary component
        complex<double>* aPoly;
        vector<complex<double>> vPoly;
	public:
		Polynomial(vector<double> realPart, vector<double> imgPart, int size, complex<double>* blankPoly, bool isArr);
		Polynomial() = default; // pass in default constructor
		~Polynomial() = default;
        void arrPoly();
        void vectPoly();
		complex<double>* getArrPoly();
        vector<complex<double>> getVectPoly();
};

#endif