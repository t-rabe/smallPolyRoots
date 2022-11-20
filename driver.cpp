#include <iomanip>
#include <iostream>
#include <complex>
// #include <string>
#include <vector>
// #include <random>
// #include <cmath>
#include <thread>
// #include <functional>
#include <fstream>

#include "./tools.hpp"
#include "./Polynomial.hpp"
#include "./RandCoeffs.hpp"
#include "./PolyEval.hpp"

using namespace std;

/**
 * @brief MAKE TOLERANCE SMALL, then factor out the found roots and repeat the process
 * adding to the data with each iteration so that it is more contoured
 * 
 * @return int 
 */

int main()
{
    auto start = std::chrono::high_resolution_clock::now();
    Tools kit;
    int sideLen = 1800;
    int polySize = 4000;
    int numSamples = 15;
    bool polyIsArr = false;
    int poly[] = {2, -6, 2, -1, 1, -1, 1};
    float poly3[] = {2., -6., 2., -1., 1., -1., 1.};
    int poly4[] = {2, -6, 2, -1, 1, -1, 1, 1,-1, 2, -6, 2, -1, 1, -1, 1, 1,-1};
    vector<float> poly5 = {1., -1., 1., -1., 2., -6., 2.};
    float poly6[] = {2., -6., 2., -1., 1., -1., 1.};
    // float poly6[] = {1., -1., 1., -1., 2., -6., 2.};
    // float poly6[] = {1., -6., 11., -6.};
    // float poly6[] = {-6., 11., -6., 1.};
    int x = 3;
    float x3 = 3.;
    float x6 = 3.;
    complex<double> poly2[] = {{2.,0.}, {-6.,0.}, {2.,0.}, {-1.,0.}, {1.,0.}, {-1.,0.}, {1.,0.}, {-6.,0.},
                                {2.,0.}, {-6.,0.}, {2.,0.}, {-1.,0.}, {1.,0.}, {-1.,0.}, {1.,0.}, {-6.,0.},
                                {2.,0.}, {-6.,0.}, {2.,0.}, {-1.,0.}, {1.,0.}, {-1.,0.}, {1.,0.}, {-6.,0.},
                                {2.,0.}, {-6.,0.}, {2.,0.}, {-1.,0.}, {1.,0.}, {-1.,0.}, {1.,0.}, {-6.,0.},
                                {2.,0.}, {-6.,0.}, {2.,0.}, {-1.,0.}, {1.,0.}, {-1.,0.}, {1.,0.}, {-6.,0.},
                                {2.,0.}, {-6.,0.}, {2.,0.}, {-1.,0.}, {1.,0.}, {-1.,0.}, {1.,0.}, {-6.,0.},
                                {2.,0.}, {-6.,0.}, {2.,0.}, {-1.,0.}, {1.,0.}, {-1.,0.}, {1.,0.}, {-6.,0.},
                                {2.,0.}, {-6.,0.}, {2.,0.}, {-1.,0.}, {1.,0.}, {-1.,0.}, {1.,0.}, {-6.,0.},
                                {2.,0.}, {-6.,0.}, {2.,0.}, {-1.,0.}, {1.,0.}, {-1.,0.}, {1.,0.}, {-6.,0.},
                                {2.,0.}, {-6.,0.}, {2.,0.}, {-1.,0.}, {1.,0.}, {-1.,0.}, {1.,0.}, {-6.,0.},
                                {2.,0.}, {-6.,0.}, {2.,0.}, {-1.,0.}, {1.,0.}, {-1.,0.}, {1.,0.}, {-6.,0.},
                                {2.,0.}, {-6.,0.}, {2.,0.}, {-1.,0.}, {1.,0.}, {-1.,0.}, {1.,0.}, {-6.,0.},
                                {2.,0.}, {-6.,0.}, {2.,0.}, {-1.,0.}, {1.,0.}, {-1.,0.}, {1.,0.}, {-6.,0.},
                                {2.,0.}, {-6.,0.}, {2.,0.}, {-1.,0.}, {1.,0.}, {-1.,0.}, {1.,0.}, {-6.,0.},
                                {2.,0.}, {-6.,0.}, {2.,0.}, {-1.,0.}, {1.,0.}, {-1.,0.}, {1.,0.}, {-6.,0.},
                                {2.,0.}, {-6.,0.}, {2.,0.}, {-1.,0.}, {1.,0.}, {-1.,0.}, {1.,0.}, {-6.,0.},
                                {2.,0.}, {-6.,0.}, {2.,0.}, {-1.,0.}, {1.,0.}, {-1.,0.}, {1.,0.}, {-6.,0.},
                                {2.,0.}, {-6.,0.}, {2.,0.}, {-1.,0.}, {1.,0.}, {-1.,0.}, {1.,0.}, {-6.,0.},
                                {2.,0.}, {-6.,0.}, {2.,0.}, {-1.,0.}, {1.,0.}, {-1.,0.}, {1.,0.}, {-6.,0.},
                                {2.,0.}, {-6.,0.}, {2.,0.}, {-1.,0.}, {1.,0.}, {-1.,0.}, {1.,0.}, {-6.,0.},
                                {2.,0.}, {-6.,0.}, {2.,0.}, {-1.,0.}, {1.,0.}, {-1.,0.}, {1.,0.}, {-6.,0.},
                                {2.,0.}, {-6.,0.}, {2.,0.}, {-1.,0.}, {1.,0.}, {-1.,0.}, {1.,0.}, {-6.,0.},
                                {2.,0.}, {-6.,0.}, {2.,0.}, {-1.,0.}, {1.,0.}, {-1.,0.}, {1.,0.}, {-6.,0.},
                                {2.,0.}, {-6.,0.}, {2.,0.}, {-1.,0.}, {1.,0.}, {-1.,0.}, {1.,0.}, {-6.,0.},
                                {2.,0.}, {-6.,0.}, {2.,0.}, {-1.,0.}, {1.,0.}, {-1.,0.}, {1.,0.}, {-6.,0.},
                                {2.,0.}, {-6.,0.}, {2.,0.}, {-1.,0.}, {1.,0.}, {-1.,0.}, {1.,0.}, {-6.,0.},
                                {2.,0.}, {-6.,0.}, {2.,0.}, {-1.,0.}, {1.,0.}, {-1.,0.}, {1.,0.}, {-6.,0.}};
    complex<double> x2 {3.,0};
    int n = sizeof(poly)/sizeof(poly[0]);
    int n2 = sizeof(poly2)/sizeof(poly2[0]);
    int n3 = sizeof(poly3)/sizeof(poly3[0]);
    int n4 = sizeof(poly4)/sizeof(poly4[0]);
    int n5 = poly5.size();
    int n6 = sizeof(poly6)/sizeof(poly6[0]);
    // for (int i=0; i<20000000; i++) {
    //     float xTest = i / 10000000.0;
    //     if (abs(horner3(poly3, n3, xTest)) < 0.0000001) {
    //         cout << horner3(poly3, n3, xTest) << "    " << xTest << endl;
    //     }
    // }
    complex<double>* polyToUse = new complex<double>[polySize];
    vector<complex<double>> vectPolyToUse;
    RandCoeffs realcoeffs(polySize,1,1);
    
    RandCoeffs imgcoeffs(polySize,1,2);


    // // THIS CREATES NEW COEFFS AND WRITES THEM TO A FILE TO BE REUSED
    // vector<double> realP = realcoeffs.getTotSample();
    // vector<double> imgP = imgcoeffs.getTotSample();
    // ofstream coeffFile;
    // coeffFile.open("testCoeffs.csv");
    // for (int m=0; m<realP.size(); m++) {
    //     coeffFile << realP[m] << '\n';
    // }
    // for (int n=0; n<imgP.size() -1; n++) {
    //     coeffFile << imgP[n] << '\n';
    // }
    // coeffFile << imgP[imgP.size() -1];
    // coeffFile.close();

    // // THIS IS THE ALTERNATIVE TO THE ABOVE
    // // IT JUST READS PREVIOUSLY SAVED COEFFS FROM A FILE CALLED "testCoeffs.csv"
    vector<double> realP;
    vector<double> imgP;
    ifstream coeffFile;
    coeffFile.open("testCoeffs.csv");
    string coeff;
    double coeffDoub;
    int lineNum = 0;
    if (coeffFile.is_open()) {
        while (coeffFile) {
            getline(coeffFile,coeff);
            coeffDoub = stod(coeff);
            if (lineNum < polySize) {
                realP.push_back(coeffDoub);
            }
            else {
                imgP.push_back(coeffDoub);
            }
            lineNum ++;
        }
    }
    // cout << realP.size() << " " << imgP.size() << endl;

    Polynomial polynomial(realP, imgP, polySize, polyToUse, polyIsArr);
    polyToUse = polynomial.getArrPoly();
    vectPolyToUse = polynomial.getVectPoly();
    // int nToUse = sizeof(polyToUse); //sizeof(polyToUse[0]);
    // cout << nToUse << " here\n";

    ofstream myFile;
    std::cout << "Done here." << endl;

    auto start2 = std::chrono::high_resolution_clock::now();
    // for (int i=0; i<1000000; i++) {
    //     kit.horner3(poly3, n3, x3);
    // }
    auto end2 = std::chrono::high_resolution_clock::now();
    auto start3 = std::chrono::high_resolution_clock::now();
    // for (int i=0; i<1000000; i++) {
    //     kit.horner2(polyToUse, polySize, x2);
    // }
    auto end3 = std::chrono::high_resolution_clock::now();
    // cout << "Value of polynomial is " << kit.horner(poly, n, x) << endl;
    // cout << "Value of polynomial is " << kit.horner3(poly3, n3, x3) << endl;
    // cout << "Value of polynomial is ";
    // for (int j=0; j<n6; j++) {
    //     cout << kit.horner5(poly6, n6, x6)[j] << " ";
    // }
    // cout << endl;
    
    // complex<double> tempEval = {0.,0.};
    // vector<complex<double>> tempCoeffs;
    vector<double> realSpaced = kit.linspace(-1.1,1.1,sideLen);
    vector<double> imgSpaced = kit.linspace(-1.1,1.1,sideLen);
    // vector<double> realSpaced = kit.linspace(.55,.75,sideLen);
    // vector<double> imgSpaced = kit.linspace(-.85,-.65,sideLen);
    // complex<double> currLoc = {0.,0.};
    // int rootNum = 0;
    // bool going = true;
    // int binNum = 0;
    // int numLoops = 0;
    // int imgSize = imgSpaced.size(); // the number of imaginary bins (so as not to keep calculating it)
    // vector<int> binCountVect((sideLen*sideLen),0);
    // int binCountArr[sideLen *sideLen] = {0};

    PolyEval polyeval(vectPolyToUse,kit,realSpaced,imgSpaced,polySize,sideLen,numSamples);
    vector<vector<int>> binCountVect = polyeval.getBinCount();

    myFile.open("test40.csv");
    // tempCoeffs = kit.horner6(vectPolyToUse, polySize, {1.,0.2});
    // for (int j=0; j<polySize; j++) {
    //     cout << tempCoeffs[j] << '\n';
    // }
    // cout << tempCoeffs.back() << endl;

    // while (going) {
    //     for (int k=0; k<realSpaced.size(); k++) {
    //         for (int m=0; m<imgSize; m++) {
    //             binNum = (k*imgSize) +m;
    //             currLoc = {realSpaced[k],imgSpaced[m]};
    //             tempCoeffs = kit.horner6(vectPolyToUse, polySize, currLoc);
    //             if (abs(tempCoeffs.back()) < 0.1) {
    //             // if (abs(tempCoeffs.back()) < (0.1 + (numLoops /10.0))) {
    //             //     cout << rootNum << " real part = " << realSpaced[k] << " img part = " << imgSpaced[m] <<
    //             // " value = " << tempCoeffs.back() << '\n';
    //                 // binCountArr[binNum] ++;
    //                 binCountVect[binNum] ++;
    //                 // tempCoeffs.pop_back();
    //                 // vectPolyToUse = tempCoeffs;
    //                 rootNum++;
    //             }
    //         }
            
    //         // tempEval = horner3(poly3, n3, spacedNums[k]);
    //         // if (abs(tempEval) < 0.00001) {
    //         //     cout << "x val = " << spacedNums[k] << " value = " << tempEval << endl;
    //         // }
    //     }
    //     numLoops ++;
    //     std::cout << "Loop # " << numLoops << "\nnumRoots " << rootNum << '\n';
    //     if (rootNum >= 900) {
    //         going = false;
    //     }
    //     else if (numLoops >= 1) {
    //         going = false;
    //     }
    // }

    // for (int k=0; k<realSpaced.size(); k++) {
    //     for (int m=0; m<imgSpaced.size(); m++) {
    //         currLoc = {realSpaced[k],imgSpaced[m]};
    //         tempEval = kit.horner2(polyToUse, polySize, currLoc);
    //         if (abs(tempEval) < 0.1) {
    //             rootNum++;
                // cout << rootNum << " real part = " << realSpaced[k] << " img part = " << imgSpaced[m] <<
                // " value = " << tempEval << '\n';
    //             myFile << 1 << '\n';
    //         }
    //         else {
    //             myFile << 0 << '\n';
    //         }
    //     }
        
    //     // tempEval = horner3(poly3, n3, spacedNums[k]);
    //     // if (abs(tempEval) < 0.00001) {
    //     //     cout << "x val = " << spacedNums[k] << " value = " << tempEval << endl;
    //     // }
    // }

    // cout << "Root num: " << rootNum << '\n';
    for (int g=0; g<sideLen; g++) {
        myFile << binCountVect[g][0];
        for (int h=1; h<sideLen; h++) {
            myFile << ',' << binCountVect[g][h];
        }
        // myFile << binCountArr[g] << '\n';
        myFile << '\n';
    }
    myFile.close();
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end-start);
    auto duration2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end2-start2);
    auto duration3 = std::chrono::duration_cast<std::chrono::nanoseconds>(end3-start3);
    std::cout << "Iteration time: " << duration2.count() / 1000000.0 << "ns\nIteration2 time: " << 
        duration3.count() / 1000000.0 << "ns\nTotal time: " << duration.count() << "s\n";
    return 0;
}