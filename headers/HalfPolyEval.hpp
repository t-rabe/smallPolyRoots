#ifndef HALFPOLYEVAL_H
#define HALFPOLYEVAL_H
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

class HalfPolyEval{
	private:
		int numPoints; // total # of points to plot (length of the total random vect)
        int numSamples; // # of vectors for the total vect to be split into
		int genNum; // used to differentiate between real and img parts so they dont get same vals
		int imgSize;
        int twoImgSize;
        int realSpan;
		int polySpan;
		int polySize;
		int numOfPoly;
		double bigBoy;
        vector<complex<double>> vectPolyToUse;
        Tools kit;
        vector<double> realSpaced;
        vector<double> imgSpaced;
        vector<vector<int>> binCountVect; // a large vect of rand ordered indices split into shorter vects
        vector<vector<int>> binCountVect2;
		vector<vector<int>> binCountVect3;
		vector<vector<int>> maxBinCountVect;
		vector<vector<double>> pixelValVect; // a vect of vects to store absolute val of poly evaluated at each pixel
		vector<vector<int>> colMinVect;
		vector<vector<int>> rowMinVect;
		vector<vector<int>> colMaxVect;
		vector<vector<int>> rowMaxVect;
		vector<vector<int>> minPeaksVect;
		vector<vector<vector<float>>> multPixValVect;
		vector<vector<vector<int>>> multBinValVect;
		void threadSafe_Sample();
		void threadSafe_Sample2();
		void threadSafe_Sample3();
		void threadSafe_Sample4();
		void threadSafe_Sample5();
		void threadSafe_Sample6();
		void threadSafe_Sample7();
		void evalPixel(int startReal, int endReal);
		void createMat(int startReal, int endReal);
		void createMat2(int startReal, int endReal);
		void findAllMins(int startPoly, int endPoly);
		void findColMin(int startReal, int endReal);
		void findRowMin(int startImg, int endImg);
		void findColMax(int startReal, int endReal);
		void findRowMax(int startImg, int endImg);
		void findMinPeaks(int startReal, int endReal);
		void combineColRow();
		void combineColRow2();
		void combineColRow3();
	public:
		HalfPolyEval(vector<complex<double>> vectPolyToUse_, Tools kit_, vector<double> realSpaced_, vector<double> imgSpaced_,
					int polySize_, int imgSize_, int numSamples_, int numOfPoly_, double bigBoy_);
		HalfPolyEval() = default; // pass in default constructor
		~HalfPolyEval() = default;
		vector<vector<int>> getBinCount();
		vector<vector<int>> getBinCount2();
		vector<vector<int>> getBinCount3();
		vector<vector<int>> getMaxBinCount();
		vector<vector<double>> getPixelVal();
		vector<vector<int>> getMinPeaks();
};

#endif