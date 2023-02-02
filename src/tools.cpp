#include "../headers/tools.hpp"

// returns value of poly[0]x(n-1) + poly[1]x(n-2) + .. + poly[n-1]
int Tools::horner(int poly[], int n, int x)
{
    int result = poly[0]; // Initialize result
 
    // Evaluate value of polynomial using Horner's method
    for (int i=1; i<n; i++)
        result = result*x + poly[i];
 
    return result;
}

complex<double> Tools::horner2(complex<double> poly2[], int n2, complex<double> x2)
{
    complex<double> result2 = poly2[0]; // Initialize result
 
    // Evaluate value of polynomial using Horner's method
    for (int i=1; i<n2; i++)
        result2 = result2*x2 + poly2[i];
 
    return result2;
}

float Tools::horner3(float poly3[], int n3, float x3)
{
    float result3 = poly3[0]; // Initialize result
 
    // Evaluate value of polynomial using Horner's method
    for (int i=1; i<n3; i++)
        result3 = result3*x3 + poly3[i];
 
    return result3;
}

// slower version to get new polynomial (use horner5 instead of this)
vector<float> Tools::horner4(vector<float> poly5, int n5, float x5) {
    vector<float> vals;
    float temp = poly5.back();
    poly5.pop_back();
    vals.push_back(temp);

    while (!poly5.empty()) {
        temp = (temp *x5) + poly5.back();
        poly5.pop_back();
        vals.push_back(temp);
    }
    return vals;
}

// This is the fastest function (so far) for getting
// the new polynomial after factoring out a root
float* Tools::horner5(float poly6[], int n6, float x6) {
    float *coeffs = new float[n6 -1];
    
    float temp = poly6[0];
    coeffs[0] = temp;

    for (int i=1; i<(n6 -1); i++) {
        temp = (temp *x6) + poly6[i];
        coeffs[i] = temp;
    }
    return coeffs;
}

// this does the same thing as the previous function, but it uses complex instead of float
// RETURNS THE REMAINDER AS PART OF THE FUNCTION
vector<complex<double>> Tools::horner6(vector<complex<double>> poly7, int n7, complex<double> x7) {
    vector<complex<double>> coeffs;
    complex<double> temp = poly7[0];
    coeffs.push_back(temp);
    for (int i=1; i<(n7); i++) {
    // for (int i=(n7-2); i>=0; i--) {
        temp = (temp *x7) + poly7[i];
        coeffs.push_back(temp);
    }
    return coeffs;
}

// just like horner2 BUT it takes in a complex vector instead of complex array
complex<double> Tools::horner7(vector<complex<double>> poly7, int n7, complex<double> x7, double bigBoy)
{
    complex<double> result7 = poly7[0]; // Initialize result
    // double bigBoy = pow(10,200);
    // Evaluate value of polynomial using Horner's method
    for (int i=1; i<n7; i++)
        result7 = result7*x7 + poly7[i];
    if (abs(result7) < bigBoy) {
        return result7;
    }
    else {
        return 0;
    }
}

// used to evaluate multiple polynomials at a time and return an array of all results
// at a given pixel. Splits poly8 into numPoly sections of length n8.
float* Tools::horner8(vector<complex<double>> poly8, int numPoly8, int n8, complex<double> x8, double bigBoy)
{
    float *valsAtPix = new float[numPoly8];
    int index = 0;
    
    for (int k=0; k<numPoly8; k++) {
        complex<double> tempResult = poly8[index]; // the first value of each sub polynomial
        // Evaluate value of polynomial using Horner's method
        for (int i=1; i<n8; i++) {
            index++;
            tempResult = tempResult*x8 + poly8[index];
        }
        valsAtPix[k] = abs(tempResult);
        index++;
    }
    
    return valsAtPix;
}

// same as above except RETURNS A VECTOR
vector<float> Tools::horner9(vector<complex<double>> poly9, int numPoly9, int n9, complex<double> x9, double bigBoy)
{
    vector<float> valsAtPix2;
    complex<double> tempResult2;
    // float result9 = 0.0;
    int index2 = 0;
    
    for (int k=0; k<numPoly9; k++) {
        // index2 = k *n9;
        tempResult2 = poly9[index2]; // the first value of each sub polynomial
        // cout << tempResult2.real() << " " << tempResult2.imag() << endl;
        // Evaluate value of polynomial using Horner's method
        for (int i=1; i<n9; i++) {
            index2++;
            tempResult2 = tempResult2*x9 + poly9[index2];
            
        }
        index2++;
        valsAtPix2.push_back(abs(tempResult2));
        // result9 = abs(tempResult2);
        // if (result9 < bigBoy) {
        //     valsAtPix2[k] = result9;
        // }
        // else {
        //     valsAtPix2[k] = 0.0;
        // }
    }
    
    return valsAtPix2;
}

// Used to get evenly spaced numbers over a certain range
vector<double> Tools::linspace(double start_in, double end_in, int num_in) {
	vector<double> linspaced;
	
	double start = start_in;
	double end = end_in;
	double num = num_in *1.0;
	
	if (num_in == 0) { 
		return linspaced; 
	}
	if (num_in == 1) {
		linspaced.push_back(start);
		return linspaced;
	}
	
	double delta = (end - start) / (num - 1);
	for(int i=0; i < num_in; ++i) {
		linspaced.push_back((start + delta * i));
	}
	
	return linspaced;
}

void Tools::writeBlankFile(string whichFile, int sideLen) {
    ofstream myFile;
    myFile.open(whichFile);
    for (int g=0; g<sideLen; g++) {
        myFile << "0";
        for (int h=1; h<sideLen; h++) {
            myFile << ",0";
        }
        myFile << '\n';
    }
    myFile.close();
}