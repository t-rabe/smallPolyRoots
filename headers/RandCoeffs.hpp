#ifndef RANDCOEFFS_H
#define RANDCOEFFS_H
#define _USE_MATH_DEFINES

#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <cmath>
#include <thread>
#include <numeric>

using namespace std;

/**
 *  Class that models a sequence of i.i.d. standard gaussians
 *  vector of distributions that can be sampled with random seeds, giving
 *  a model of uncorrelated Gaussian variables.
 */
class RandCoeffs{
	private:
		int numPoints; // total # of points to plot (length of the total random vect)
        int numSamples; // # of vectors for the total vect to be split into
		int genNum; // used to differentiate between real and img parts so they dont get same vals
		vector<vector<double>> samples; // a large vect of rand ordered indices split into shorter vects
		vector<vector<double>> samplesOnes;
		void threadSafe_Sample();
		void threadSafe_SampleOnes();
		void getRandomValue(unsigned int seed, int index);
		void getRandomValueOnes(unsigned int seed, int index);
	public:
		RandCoeffs(int numPoints_, int numSamples_, int genNum_);
		RandCoeffs() = default; // pass in default constructor
		~RandCoeffs() = default;
		vector<vector<double>> getSample();
        vector<double> getTotSample();
		vector<vector<double>> getSampleOnes();
        vector<double> getTotSampleOnes();
};

#endif