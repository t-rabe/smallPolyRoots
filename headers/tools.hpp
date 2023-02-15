#ifndef TOOLS_H
#define TOOLS_H
#define _USE_MATH_DEFINES

#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <cmath>
#include <thread>
#include <fstream>
#include <numeric>
#include <complex>

using namespace std;

class Tools{
	public:
		int horner(int poly[], int n, int x);
        complex<double> horner2(complex<double> poly2[], int n2, complex<double> x2);
        float horner3(float poly3[], int n3, float x3);
        vector<float> horner4(vector<float> poly5, int n5, float x5);
        float* horner5(float poly6[], int n6, float x6);
        vector<complex<double>> horner6(vector<complex<double>> poly7, int n7, complex<double> x7);
        complex<double> horner7(vector<complex<double>> poly7, int n7, complex<double> x7, double bigBoy);
        float* horner8(vector<complex<double>> poly8, int numPoly8, int n8, complex<double> x8);
        vector<float> horner9(vector<complex<double>> poly9, int numPoly9, int n9, complex<double> x9);
        vector<double> linspace(double start_in, double end_in, int num_in);
        void writeBlankFile(string whichFile_, int sideLen_);
};

#endif