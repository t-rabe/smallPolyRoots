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
#include "../headers/FlatHalfPoly.hpp"

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
    string fileNum = "1800_4000";
    string fileName = "../output/pixEval_" + fileNum + ".csv";
    
    int sideLen = 1800; // side length (in pixels) of the resulting image
    int polySize = 4000; // degree of the polynomial to be used
    int numSamples = 15;
    bool polyIsArr = false;

    Tools kit;
    double largeNum = pow(10,300); // upper bound for polyEval vals (above ~ infinity)
        
    complex<double>* polyToUse = new complex<double>[polySize];
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
     * THIS IS THE ALTERNATIVE TO THE ABOVE
     * IT JUST READS PREVIOUSLY SAVED COEFFS FROM A FILE CALLED "testCoeffs.csv"
    */ 
    vector<double> realP; // real part of each coeff
    vector<double> imgP; // complex part of each coeff
    ifstream coeffFile;
    coeffFile.open("../src/testCoeffs.csv"); // file holding old coeffs
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
            else if (lineNum < (polySize*2)) {
                // imgP.push_back(coeffDoub); // uncomment to add complex coeffs
                imgP.push_back(0.1); // uncomment to keep real coeffs
            }
            lineNum ++;
        }
    }
    Polynomial polynomial(realP, imgP, polySize, polyToUse, polyIsArr);
    polyToUse = polynomial.getArrPoly();
    vectPolyToUse = polynomial.getVectPoly();

    std::cout << "Done here. Number of Samples: " << numSamples << endl;
    
    vector<double> realSpaced = kit.linspace(-1.1,1.1,sideLen);
    vector<double> imgSpaced = kit.linspace(-1.1,1.1,sideLen);
    // vector<double> realSpaced = kit.linspace(.55,.75,sideLen);
    // vector<double> imgSpaced = kit.linspace(-.85,-.65,sideLen);
    
    PolyEval polyeval(vectPolyToUse,kit,realSpaced,imgSpaced,polySize,sideLen,numSamples,largeNum);
    HalfPolyEval halfpolyeval(vectPolyToUse,kit,realSpaced,imgSpaced,polySize,sideLen,numSamples,largeNum);
    FlatHalfPoly flathalfpoly(vectPolyToUse,kit,realSpaced,imgSpaced,polySize,sideLen,numSamples,largeNum);
    auto start2 = std::chrono::high_resolution_clock::now();
    
    ofstream myFile;
    myFile.open(fileName);

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
    // vector<vector<int>> halfBinVect = halfpolyeval.getBinCount2();
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
    vector<vector<int>> binCountVect = polyeval.getBinCount2();
    auto start3 = std::chrono::high_resolution_clock::now();
    for (int g=0; g<sideLen; g++) {
        myFile << binCountVect[g][0];
        for (int h=1; h<sideLen; h++) {
            myFile << ',' << binCountVect[g][h];
        }
        myFile << '\n';
    }

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

    // ofstream myFile2;
    // myFile2.open("../output/SECONDFILENAME.csv"); // change file name

    /**
     * copy one of the above functions here to write two files per batch!
     */

    // myFile2.close();
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(start3-start2);
    auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(end-start3);
    auto duration3 = std::chrono::duration_cast<std::chrono::seconds>(end-start);

    std::cout << "Get bin counts: " << duration.count() << "s\nWrite file: " << duration2.count() <<
                "ms\nTotal time: " << duration3.count() << "s\n";
    return 0;
}
