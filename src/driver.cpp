#include <iomanip>
#include <iostream>
#include <complex>
#include <string>
#include <vector>
#include <thread>
#include <fstream>

#include "../headers/tools.hpp"
#include "../headers/Polynomial.hpp"
#include "../headers/RandCoeffs.hpp"
#include "../headers/PolyEval.hpp"
#include "../headers/HalfPolyEval.hpp"
#include "../headers/QuartPolyEval.hpp"
#include "../headers/FlatHalfPoly.hpp"
#include "../headers/AxesEval.hpp"

using namespace std;

/**
 * @brief MAKE TOLERANCE SMALL, then factor out the found roots and repeat the process
 * adding to the data with each iteration so that it is more contoured
 *
 * @return int
 */

int main()
{
    auto comeca = std::chrono::high_resolution_clock::now();
    for (int i=0; i<24; i++) {
        auto start = std::chrono::high_resolution_clock::now();
        string fileNum = "2400_500_axes";

        int sideLen = 2400; // side length (in pixels) of the resulting image
        int polySize = 24; // degree of the polynomial to be used
        int numSamples = 15;
        int startPoly = 0;
        int numPolys = 500; // NEEDS TO BE LARGER THAN NUMSAMPLES !!!!!
        int offset = i; // offsets the coeff's indices. Must be smaller than 24 !!!!
        int coeffSize = (polySize *numPolys); // num coeffs to load for real/img
        bool polyIsArr = false;

        string fileName = "../output/pixEval_" + fileNum + to_string(offset+0) + ".csv";

        Tools kit;
        double largeNum = pow(10,300); // upper bound for polyEval vals (above ~ infinity)

        vector<complex<double>> vectPolyToUse;

        /**
         * THIS CREATES NEW COEFFS AND WRITES THEM TO A FILE TO BE REUSED
        */
        // RandCoeffs realcoeffs(polySize,1,1);
        // RandCoeffs imgcoeffs(polySize,1,2);
        // vector<double> realP = realcoeffs.getTotSample();
        // vector<double> imgP = imgcoeffs.getTotSample();
        // ofstream coeffFile0;
        // coeffFile0.open("../src/testCoeffs3.csv");
        // for (int m=0; m<realP.size(); m++) {
        //     coeffFile0 << realP[m] << '\n';
        // }
        // for (int n=0; n<imgP.size() -1; n++) {
        //     coeffFile0 << imgP[n] << '\n';
        // }
        // coeffFile0 << imgP[imgP.size() -1];
        // coeffFile0.close();

        /**
         * THIS CREATES NEW COEFFS AND WRITES THEM TO A FILE TO BE REUSED
         * EXCLUSIVELY 1 OR -1, ONE MILLION OF THEM
        */
        // RandCoeffs onescoeffs(1000000,1,1);
        // vector<double> realP = onescoeffs.getTotSampleOnes();
        // ofstream coeffFile0;
        // coeffFile0.open("../src/testCoeffsOnes.csv");
        // for (int m=0; m<realP.size() -1; m++) {
        //     coeffFile0 << realP[m] << '\n';
        // }
        // coeffFile0 << realP[realP.size() -1];
        // coeffFile0.close();


        /**
         * THIS IS THE ALTERNATIVE TO THE ABOVE
         * IT JUST READS PREVIOUSLY SAVED COEFFS FROM A FILE CALLED "testCoeffs.csv"
        */
        vector<double> realP; // real part of each coeff
        vector<double> imgP; // complex part of each coeff
        ifstream coeffFile;
        coeffFile.open("../src/testCoeffsOnesBigRandom.csv"); // file holding old coeffs
        string coeff;
        double coeffDoub;
        int lineNum = 0;
        if (coeffFile.is_open()) {
            while (coeffFile) {
                getline(coeffFile,coeff);
                coeffDoub = stod(coeff);
                if ((lineNum>(offset+startPoly)) && (lineNum < (coeffSize+offset+startPoly))) {
                    realP.push_back(coeffDoub);
                    imgP.push_back(0.0); // uncomment to keep real coeffs
                }
                // else if (lineNum < (coeffSize*2)) {
                //     // imgP.push_back(coeffDoub); // uncomment to add complex coeffs
                //     imgP.push_back(0.0); // uncomment to keep real coeffs
                // }
                lineNum ++;
            }
        }
        Polynomial polynomial(realP, imgP, coeffSize, polyIsArr);
        vectPolyToUse = polynomial.getVectPoly();

        std::cout << "Done here. Number of Samples: " << numSamples << "\nFile number: " << offset << endl;

        vector<double> realSpaced = kit.linspace(-1.6,1.6,sideLen);
        vector<double> imgSpaced = kit.linspace(-1.6,1.6,sideLen);
        // vector<double> realSpaced = kit.linspace(-1.55,1.55,sideLen);
        // vector<double> imgSpaced = kit.linspace(-1.55,1.55,sideLen);
        // vector<double> realSpaced = kit.linspace(0.55,0.65,sideLen);
        // vector<double> imgSpaced = kit.linspace(0.25,0.35,sideLen);
        
        // PolyEval polyeval(vectPolyToUse,kit,realSpaced,imgSpaced,polySize,sideLen,numSamples,numPolys,largeNum);
        // HalfPolyEval halfpolyeval(vectPolyToUse,kit,realSpaced,imgSpaced,polySize,sideLen,numSamples,numPolys,largeNum);
        // QuartPolyEval quartpolyeval(vectPolyToUse,kit,realSpaced,imgSpaced,polySize,sideLen,numSamples,numPolys,largeNum);
        // FlatHalfPoly flathalfpoly(vectPolyToUse,kit,realSpaced,imgSpaced,polySize,sideLen,numSamples,largeNum);
        AxesEval axeseval(vectPolyToUse,kit,realSpaced,imgSpaced,polySize,sideLen,numSamples,numPolys,largeNum);
        ofstream myFile;
        myFile.open(fileName);
        auto start2 = std::chrono::high_resolution_clock::now();
        
        /**
         * @brief creates the minCol/minRow vector
         *
         * @return FLAT VECTOR and ONLY HALF of the image
         * @note INNEFICIENT AND CAN ONLY BE USED WITH COMPLEX CONJ. ROOT THRM. (CCRT)
         */
        // vector<int> flatBinVect = flathalfpoly.getBinCount2();
        // int sizeOfVect = sideLen *sideLen /2;
        // auto start3 = std::chrono::high_resolution_clock::now();
        // for (int g=0; g<sizeOfVect; g++) {
        //     myFile << flatBinVect[g] << '\n';
        // }

        /**
         * @brief creates the minCol/minRow vector
         *
         * @return 2D VECTOR and ONLY HALF of the image
         * @note VERY EFFICIENT AND CAN ONLY BE USED WITH COMPLEX CONJ. ROOT THRM. (CCRT)
         */
        // vector<vector<int>> halfBinVect = halfpolyeval.getBinCount3();
        // auto start3 = std::chrono::high_resolution_clock::now();
        // for (int g=0; g<sideLen; g++) {
        //     myFile << halfBinVect[g][0];
        //     for (int h=1; h<(sideLen/2); h++) {
        //         myFile << ',' << halfBinVect[g][h];
        //     }
        //     for (int i=((sideLen/2)-1); i>=0; i--) {
        //         myFile << ',' << halfBinVect[g][i];
        //     }
        //     myFile << '\n';
        // }

        /**
         * @brief creates the minCol/minRow ARRAY
         *
         * @return 2D VECTOR and ONLY HALF of the image same as above but returns an ARRAY
         * @note VERY EFFICIENT AND CAN ONLY BE USED WITH COMPLEX CONJ. ROOT THRM. (CCRT)
         */
        // int** halfBinVect = halfpolyeval.getBinCount4();
        // auto start3 = std::chrono::high_resolution_clock::now();
        // for (int g=0; g<sideLen; g++) {
        //     myFile << halfBinVect[g][0];
        //     for (int h=1; h<(sideLen/2); h++) {
        //         myFile << ',' << halfBinVect[g][h];
        //     }
        //     for (int i=((sideLen/2)-2); i>=0; i--) {
        //         myFile << ',' << halfBinVect[g][i];
        //     }
        //     myFile << ",0\n";
        // }

        /**
         * @brief creates the local min vector
         *
         * @return 2D VECTOR and ONLY QUARTER of the image
         * @note MOST EFFICIENT AND CAN ONLY BE USED WITH COMPLEX CONJ. ROOT THRM. (CCRT)
         * 
         * @note MUST BE USED WITH LARGE SAMPLES BC INDIVIDUAL POLYS ARE NOT SYMMETRIC
         *       OVER THE IMAGINARY AXIS !!!!!!!!!!!!!!
         */
        // vector<vector<unsigned short int>> quartBinVect = quartpolyeval.getBinCount3();
        // auto start3 = std::chrono::high_resolution_clock::now();
        // for (int g=0; g<(sideLen/2); g++) {
        //    myFile << quartBinVect[g][0];
        //    for (int h=1; h<(sideLen/2); h++) {
        //        myFile << ',' << quartBinVect[g][h];
        //    }
        //    for (int i=((sideLen/2)-2); i>=0; i--) {
        //        myFile << ',' << quartBinVect[g][i];
        //    }
        //    myFile << "\n";
        // }
        // for (int g=((sideLen/2)-2); g>=0; g--) {
        //    myFile << quartBinVect[g][0];
        //    for (int h=1; h<(sideLen/2); h++) {
        //        myFile << ',' << quartBinVect[g][h];
        //    }
        //    for (int i=((sideLen/2)-2); i>=0; i--) {
        //        myFile << ',' << quartBinVect[g][i];
        //    }
        //    myFile << "\n";
        // }

        /**
         * @brief creates the local min vector
         *
         * @return 2D VECTOR and ONLY AXES of the image
         * @note USED WITH THE FUNCTION ABOVE 
         *       AND CAN ONLY BE USED WITH COMPLEX CONJ. ROOT THRM. (CCRT)
         * 
         * @note MUST BE USED WITH LARGE SAMPLES BC INDIVIDUAL POLYS ARE NOT SYMMETRIC
         *       OVER THE IMAGINARY AXIS !!!!!!!!!!!!!!
         */
        vector<vector<unsigned short int>> axesBinVect = axeseval.getBinCount();
        auto start3 = std::chrono::high_resolution_clock::now();
        for (int g=0; g<(sideLen-1); g++) {
           myFile << axesBinVect[g][0];
           for (int h=1; h<(sideLen-1); h++) {
               myFile << ',' << axesBinVect[g][h];
           }
           myFile << "\n";
        }

        /**
         * @brief creates a vector of min peaks (the isolated points and values of 2)
         *        from the minCol/minRow combined vect
         *
         * @return 2D VECTOR and FULL image
         * @note OVERCOUNTS BECAUSE OF FADING CONTOUR LINES (SOME ISOLATED POINTS ARENT ROOTS)
         */
        // vector<vector<int>> minPeakValVect = polyeval.getMinPeaks();
        // auto start3 = std::chrono::high_resolution_clock::now();
        // for (int g=0; g<sideLen; g++) {
        //     myFile << minPeakValVect[g][0];
        //     for (int h=1; h<sideLen; h++) {
        //         myFile << ',' << minPeakValVect[g][h];
        //     }
        //     myFile << '\n';
        // }

        /**
         * @brief creates the minCol/minRow vector for non-symmetric images
         *
         * @return 2D VECTOR and FULL image
         * @note LESS EFICIENT BUT CAN BE USED WITHOUT COMPLEX CONJ. ROOT THRM. (CCRT)
         */
        // vector<vector<int>> binCountVect = polyeval.getBinCount3();
        // auto start3 = std::chrono::high_resolution_clock::now();
        // for (int g=0; g<sideLen; g++) {
        //     myFile << binCountVect[g][0];
        //     for (int h=1; h<sideLen; h++) {
        //         myFile << ',' << binCountVect[g][h];
        //     }
        //     myFile << '\n';
        // }

        /**
         * @brief creates the minCol/minRow vector for non-symmetric images
         *
         * @return 2D VECTOR and FULL image AS ARRAY
         * @note LESS EFICIENT BUT CAN BE USED WITHOUT COMPLEX CONJ. ROOT THRM. (CCRT)
         */
        // int** binCountArr = polyeval.getBinCount4();
        // auto start3 = std::chrono::high_resolution_clock::now();
        // for (int g=0; g<sideLen; g++) {
        //     myFile << binCountArr[g][0];
        //     for (int h=1; h<sideLen; h++) {
        //         myFile << ',' << binCountArr[g][h];
        //     }
        //     myFile << '\n';
        // }

        /**
         * @brief creates a vect of pixel values (essentially a better topo map)
         *
         * @return 2D VECTOR and FULL image
         * @note NEEDS TO BE TRIMMED AROUND THE EDGES IN PYTHON BECAUSE OF OVER-SATURATION
         */
        // vector<vector<double>> pixValVect = polyeval.getPixelVal();
        // auto start3 = std::chrono::high_resolution_clock::now();
        // for (int g=0; g<sideLen; g++) {
        //     myFile << pixValVect[g][0];
        //     for (int h=1; h<sideLen; h++) {
        //         myFile << ',' << pixValVect[g][h];
        //     }
        //     myFile << '\n';
        // }

        myFile.close();

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(start3-start2);
        auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(end-start3);
        auto duration3 = std::chrono::duration_cast<std::chrono::seconds>(end-start);

        std::cout << "Get bin counts: " << duration.count() << "s\nWrite file: " << duration2.count() <<
                    "ms\nTotal time: " << duration3.count() << "s\n";
    }
    auto fim = std::chrono::high_resolution_clock::now();
    auto duration4 = std::chrono::duration_cast<std::chrono::minutes>(fim-comeca);
    std::cout << "FINAL RUNTIME: " << duration4.count() << " mins\n";

    return 0;
}
