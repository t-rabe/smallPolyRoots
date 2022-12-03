#ifndef FLATHALFPOLY_H
#define FLATHALFPOLY_H
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

class FlatHalfPoly{
	private:
		int numPoints; // total # of points to plot (length of the total random vect)
        int numSamples; // # of vectors for the total vect to be split into
		int genNum; // used to differentiate between real and img parts so they dont get same vals
		int imgSize;
        int twoImgSize;
        int vectSize;
        int realSpan;
		double bigBoy;
        vector<complex<double>> vectPolyToUse;
        Tools kit;
        vector<double> realSpaced;
        vector<double> imgSpaced;
        int polySize;
        vector<int> binCountVect; // a large vect of rand ordered indices split into shorter vects
        vector<int> binCountVect2;
		vector<int> maxBinCountVect;
		vector<double> pixelValVect; // a vect of vects to store absolute val of poly evaluated at each pixel
		vector<int> colMinVect;
		vector<int> rowMinVect;
		vector<int> colMaxVect;
		vector<int> rowMaxVect;
		vector<int> minPeaksVect;
		void threadSafe_Sample();
		void threadSafe_Sample2();
		void threadSafe_Sample3();
		void threadSafe_Sample4();
		void threadSafe_Sample5();
		void evalPixel(int startReal, int endReal);
		void createMat(int startReal, int endReal);
		void findColMin(int startReal, int endReal);
		void findRowMin(int startImg, int endImg);
		void findColMax(int startReal, int endReal);
		void findRowMax(int startImg, int endImg);
		void findMinPeaks(int startReal, int endReal);
		void combineColRow();
		void combineColRow2();
	public:
		FlatHalfPoly(vector<complex<double>> vectPolyToUse_, Tools kit_, vector<double> realSpaced_, vector<double> imgSpaced_, int polySize_, int imgSize_, int numSamples_, double bigBoy_);
		FlatHalfPoly() = default; // pass in default constructor
		~FlatHalfPoly() = default;
		vector<int> getBinCount();
		vector<int> getBinCount2();
		vector<int> getHalfBinCount2();
		vector<int> getMaxBinCount();
		vector<double> getPixelVal();
		vector<int> getMinPeaks();
};

#endif