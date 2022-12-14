#ifndef AXESEVAL_H
#define AXESEVAL_H
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

class AxesEval{
	private:
        int numSamples; // # of vectors for the total vect to be split into
		int imgSize;
        int halfImgSize;
        int realSpan;
		int polySpan;
		int polySize;
		int numOfPoly;
        double stepSize;
		double bigBoy;
        Tools kit;
        vector<complex<double>> vectPolyToUse;
        vector<double> realSpaced;
        vector<double> imgSpaced;
        vector<vector<unsigned short int>> totBinCountVect;
		vector<vector<vector<float>>> imgPixValVect;
		vector<vector<vector<float>>> realPixValVect;
		vector<vector<int>> imgMultiBinCountVect;
		vector<vector<int>> realMultiBinCountVect;
		void threadSafe_Sample();
		void threadSafe_Sample2();
		void threadSafe_Sample3();
		void threadSafe_Sample4();
		void imgEvalPixel(int startImg, int endImg);
		void realEvalPixel(int startReal, int endReal);
		void findImgMins(int startPoly, int endPoly);
		void findRealMins(int startPoly, int endPoly);
		void combineAxes();
	public:
		AxesEval(vector<complex<double>> vectPolyToUse_, Tools kit_, vector<double> realSpaced_, vector<double> imgSpaced_,
					int polySize_, int imgSize_, int numSamples_, int numOfPoly_, double bigBoy_);
		AxesEval() = default; // pass in default constructor
		~AxesEval() = default;
		vector<vector<unsigned short int>> getBinCount();
};

#endif