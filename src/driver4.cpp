#include <iomanip>
#include <iostream>
#include <complex>
#include <string>
#include <vector>
#include <thread>
#include <fstream>
#include <sstream>
#include <ctime>

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
 * @brief This is another driver for the program. It solves polynomials from an input csv file. Currently
 * designed to analyze all polynomials of degree 12 with coeffs as 1 or -1 (includes coeff
 * for x^0). Very similar to the original design of the program with no bells or whistles.
 * 
 * NOTE: This driver is made to read from a single .csv file and write
 * to a single .csv file. It updates an array saved in memory and writes the file at the end.
 * The file which is read is a 1D .csv with each row representing a new coeff.
 *
 * @return int
 */

int main(int argc, char *argv[])
{
    auto comeca0 = std::chrono::high_resolution_clock::now();
    string inFileName = argv[1];
    string outFileName = argv[2];
    string scratchFile = argv[3];
    // int bNumber = 1; // the file number to read data in from

    // ideally sideLen is an odd number so that the middle val lands on a single pixel
    int sideLen = 693; // side length (in pixels) of the resulting image. Should be a multiple of numSamples
    int polySize = 19; // degree of the polynomial to be used
    int numSamples = 63; // how many threads to use
    int numPolys = 32768; // NEEDS TO BE LARGER THAN NUMSAMPLES !!!!!
    int numIters = 16; // how many times to evaluate the polynomials (numPolys per iteration)
    int coeffSize = (polySize *numPolys); // num coeffs to load for real/img
    int lineNum = 0; // line being read from in file
    int offset = stoi(argv[4]) * polySize; // param to start further into the in file (if a part of it has already been done)
    bool polyIsArr = false;
    string coeff; // stores a single coeff at a time from in file
    // string outFileName = "243_500_matrix_2Dfile_";

    Tools kit;

    vector<vector<int>> BFV; // the final vector with all data. written to out file at the end
    for (int i=0; i<sideLen; i++) { // populate the BFV with zeros to start with
        BFV.push_back(vector<int> (sideLen,0));
    }

    vector<complex<double>> vectPolyToUse;
    vector<vector<unsigned short int>> partBinVect;
    vector<double> realSpaced = kit.linspace(-2.0,2.0,sideLen);
    vector<double> imgSpaced = kit.linspace(-2.0,2.0,sideLen);
    // vector<double> realSpaced = kit.linspace(2.5,3.3,sideLen);
    // vector<double> imgSpaced = kit.linspace(4.1,5.9,sideLen);
    vector<double> realP; // real part of each coeff
    vector<double> imgP; // complex part of each coeff
    // vector<int> allRealP; // all of the real coeffs saved in one vector
    // vector<double> allImgP; // all of the complex coeffs saved in one vector
    double largeNum = pow(10,300); // upper bound for polyEval vals (above ~ infinity)

    auto dateTime = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    cout << ctime(&dateTime) << '\n'; 
    cout << "Starting with " << numSamples << " samples...\n";
    auto comeca = std::chrono::high_resolution_clock::now();
    auto fim = std::chrono::high_resolution_clock::now();
    auto duration4 = std::chrono::duration_cast<std::chrono::minutes>(fim-comeca);
    auto comecaBatch = std::chrono::high_resolution_clock::now();
    auto fimBatch = std::chrono::high_resolution_clock::now();
    auto durationBatch = std::chrono::duration_cast<std::chrono::minutes>(fimBatch-comecaBatch);
    
    ifstream coeffFile;
    string name = "../coeffFiles/" + inFileName; // for use with normal desktop
    // string name = "../../../.." + scratchFile + "/" + inFileName; // for use on hpc

// ######################### START READ BIN FILE ###############################
    // FILE * pFile;
    // long lSize;
    // int * buffer;
    // size_t result;
    // int numCount;

    // pFile = fopen (name.c_str(), "rb" );
    // if (pFile==NULL) {fputs ("File error",stderr); exit (1);}

    // // obtain file size:
    // fseek (pFile , 0 , SEEK_END);
    // lSize = ftell (pFile);
    // rewind (pFile);

    // // allocate memory to contain the whole file:
    // buffer = (int*) malloc (sizeof(int)*lSize);
    // if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

    // // copy the file into the buffer:
    // result = fread (buffer,1,lSize,pFile);
    // if (result != lSize) {fputs ("Reading error",stderr); exit (3);}

    // numCount = (lSize / sizeof(buffer[0]));
    // vector<int> allRealP(buffer, buffer+numCount);
    // // cout << allRealP.size() << endl;
    // fclose (pFile);
    // free (buffer);
// ########################## END READ BIN FILE ################################

// ######################### START READ CSV FILE ###############################
    vector<int> allRealP; // real part of each coeff
    coeffFile.clear();
    coeffFile.open(name); // file holding old coeffs
    lineNum = 0;
    if (coeffFile.is_open()) {
        while (coeffFile && (lineNum<((coeffSize*numIters)+offset))) {
            getline(coeffFile,coeff);
            if (coeff != "") {
                allRealP.push_back(stoi(coeff));
            }
            lineNum ++;
        }
    } 
    coeffFile.close();
// ############################ END READ CSV FILE ####################################
    
    for (int b=0; b<numIters; b++) {
        comeca = std::chrono::high_resolution_clock::now();
        cout << "Starting iteration # " << to_string(b+1) << "/" << to_string(numIters) <<
            " at T = " << std::chrono::duration_cast<std::chrono::minutes>(comeca-comeca0).count()
            << " mins..." << endl;

        int startPoly = b*coeffSize;
        
        vectPolyToUse.clear();
        partBinVect.clear();

        /**
         * THIS READS PREVIOUSLY SAVED COEFFS FROM !!!RANDOM MATRICES/CHARACTERISTIC POLYNOMIALS!!!
        */
        realP.clear();
        imgP.clear();
        lineNum = offset+startPoly; // index for the first coeffs of the batch      
        
        while (lineNum < (coeffSize+offset+startPoly)) {
            realP.push_back(allRealP[lineNum] *1.0);
            imgP.push_back(0.0); // uncomment to keep real coeffs
            lineNum ++;
        }
        
        Polynomial polynomial(realP, imgP, coeffSize, polyIsArr);
        vectPolyToUse = polynomial.getVectPoly();
        PolyEval polyeval(vectPolyToUse,kit,realSpaced,imgSpaced,polySize,sideLen,numSamples,numPolys,largeNum);
        // HalfPolyEval halfpolyeval(vectPolyToUse,kit,realSpaced,imgSpaced,polySize,sideLen,numSamples,numPolys,largeNum);
        // FlatHalfPoly flathalfpoly(vectPolyToUse,kit,realSpaced,imgSpaced,polySize,sideLen,numSamples,largeNum);
        
        // auto start2 = std::chrono::high_resolution_clock::now();
        
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
         * @brief creates the minCol/minRow vector for non-symmetric images
         *
         * @return 2D VECTOR and FULL image
         * @note LESS EFICIENT BUT CAN BE USED WITHOUT COMPLEX CONJ. ROOT THRM. (CCRT)
         */
        partBinVect = polyeval.getBinCount3();
        for (int g=0; g<sideLen; g++) {
            for (int h=0; h<sideLen; h++) {
                BFV[g][h] += partBinVect[g][h];
            }
        }
        
        /**
         * @brief creates the minCol/minRow vector for non-symmetric images
         *
         * @return 2D and FULL image AS ARRAY
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
    }
    
    string fileName = "../output/" + outFileName; // for general use
    // string fileName = "../../../.." + scratchFile + "/" + outFileName; // for use with scratch on hpc
    ofstream myFile;
    myFile.open(fileName);
    for (int g=0; g<sideLen; g++) {
        myFile << BFV[g][0];
        for (int h=1; h<sideLen; h++) {
            myFile << ',' << BFV[g][h];
        }
        myFile << '\n';
    }
    myFile.close();

    auto fim0 = std::chrono::high_resolution_clock::now();
    auto duration5 = std::chrono::duration_cast<std::chrono::minutes>(fim0-comeca0);
    std::cout << "TOTAL RUNTIME: " << duration5.count() << " mins\n";
    return 0;
}