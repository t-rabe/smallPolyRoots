#ifndef POLYEVAL_H
#define POLYEVAL_H
#define _USE_MATH_DEFINES

#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <cmath>
#include <thread>
#include <numeric>
#include <complex>

#include "./tools.hpp"

using namespace std;

class PolyEval{
	private:
		int numPoints; // total # of points to plot (length of the total random vect)
        int numSamples; // # of vectors for the total vect to be split into
		int genNum; // used to differentiate between real and img parts so they dont get same vals
		int imgSize;
        int realSpan;
        vector<complex<double>> vectPolyToUse;
        Tools kit;
        vector<double> realSpaced;
        vector<double> imgSpaced;
        int polySize;
        vector<vector<int>> binCountVect; // a large vect of rand ordered indices split into shorter vects
		void threadSafe_Sample();
		void evalPixel(int startReal, int endReal);
	public:
		PolyEval(vector<complex<double>> vectPolyToUse_, Tools kit_, vector<double> realSpaced_, vector<double> imgSpaced_, int polySize_, int imgSize_, int numSamples_);
		PolyEval() = default; // pass in default constructor
		~PolyEval() = default;
		vector<vector<int>> getBinCount();
};

#endif