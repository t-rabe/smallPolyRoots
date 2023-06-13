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
 * @brief This driver is used to solve the roots which fall on the GRIDLINES 
 * of the final image. It should be used in addition to driver2.cpp in order
 * to remove the empty space between each square which is calculated in that driver
 *
 * argv[0] = ./runFile 
 * argv[1] = inputFile_#.bin (name of file with input coeffs)
 * argv[2] = outFile name (portion of the output file name)
 * argv[3] = scratchFile (part of outFile path) ONLY USED ON HPC! use "blankScratch" on desktop
 * argv[4] = offset (int for # of batches already run using current inputFile)
 * argv[5] = number of squares to a side of grid (e.g. 3 => 9 total squares in final img)
 * argv[6] = is left side ('T' or 'F' to say whether we are analyzing lines left(T) or right(F) of i-axis)
 * 
 * 
 * @return int
 */
int main(int argc, char *argv[])
{
    auto comeca0 = std::chrono::high_resolution_clock::now();
    if (argc!=7) {
        std::cout << "ERROR: WRONG INPUT PARAMETERS\n";
        return -1;
    }
    // ideally sideLen is an odd number so that the middle val lands on a single pixel
    Tools kit;
    int numSqrs = stoi(argv[5]);
    int sqrSideLen = 189;
    int sideLen = sqrSideLen * numSqrs; // side length (in pixels) of the resulting image. Should be a multiple of numSamples
    int polySize = 8; // degree of the polynomial to be used + 1 (for the coeff of x^0)
    int numSamples = 21; // how many threads to use
    int numPolys = 40000; // NEEDS TO BE LARGER THAN NUMSAMPLES !!!!!
    int numIters = 3; // how many times to evaluate the polynomials (numPolys per iteration)
    // int numPolys = 1200000; // NEEDS TO BE LARGER THAN NUMSAMPLES !!!!!
    // int numIters = 35; // how many times to evaluate the polynomials (numPolys per iteration)
    int coeffSize = (polySize *numPolys); // num coeffs to load for real/img
    int lineNum = 0; // line being read from in file
    int offset = stoi(argv[4]) * polySize * numPolys * numIters; // param to start further into the in file (if a part of it has already been done)
    bool polyIsArr = false;
    bool isVert = false;
    bool isLeft = false;
    if ((char)tolower(argv[6][0]) == 't') {
        isLeft = true;
    }
    string multiCoeffs; // stores a row at a time from in file
    string coeff; // stores a single coeff at a time from in file
    string inFileName = argv[1];
    string outFileName;
    
    double lineCenter;
    double num = sqrSideLen *1.0;
    double delta = (0.8) / (num -1.0);
    vector<double> realSpaced;
    vector<double> imgSpaced;
    vector<vector<int>> BFV; // the final vector with all data. written to out file at the end
    string scratchFile = argv[3];

    vector<complex<double>> vectPolyToUse;
    vector<vector<int>> axesBinVect;

    vector<double> realP; // real part of each coeff
    vector<double> imgP; // complex part of each coeff
    double largeNum = pow(10,300); // upper bound for polyEval vals (above ~ infinity)

    auto dateTime = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::cout << ctime(&dateTime) << '\n'; 
    std::cout << "Starting with " << numSamples << " samples...\n";
    auto comeca = std::chrono::high_resolution_clock::now();
    auto comecaBatch = std::chrono::high_resolution_clock::now();
    auto fimBatch = std::chrono::high_resolution_clock::now();
    
    ifstream coeffFile;
    string name = "../coeffFiles/" + inFileName; // for use with normal desktop
    // string name = "../../../.." + scratchFile + "/" + inFileName; // for use on hpc

// ######################### START READ BIN FILE ###############################
    FILE * pFile;
    long unsigned lSize;
    int * buffer;
    size_t result;
    int numCount;

    pFile = fopen (name.c_str(), "rb" );
    if (pFile==NULL) {fputs ("File error",stderr); exit (1);}

    // obtain file size:
    fseek (pFile , 0 , SEEK_END);
    lSize = ftell (pFile);
    rewind (pFile);

    // allocate memory to contain the whole file:
    buffer = (int*) malloc (sizeof(int)*lSize);
    if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

    // copy the file into the buffer:
    result = fread (buffer,1,lSize,pFile);
    if (result != lSize) {fputs ("Reading error",stderr); exit (3);}

    numCount = (lSize / sizeof(buffer[0]));
    vector<int> allRealP(buffer, buffer+numCount);
    // cout << allRealP.size() << endl;
    fclose (pFile);
    free (buffer);
// ########################## END READ BIN FILE ################################

    for (int r=5; r<7; r++) {
        comecaBatch = std::chrono::high_resolution_clock::now();
        realSpaced.clear();
        imgSpaced.clear();
        BFV.clear(); // the final vector with all data. written to out file at the end
        if (r<5) {
            isVert = true;
            sideLen = sqrSideLen *5;
            if (isLeft) {
                outFileName = "leftVertLine" + to_string(r+1) + "_" + string(argv[2]) + ".csv";
                lineCenter = (-3.6) + (0.4*(r+1));
            }
            else {
                outFileName = "rightVertLine" + to_string(r+1) + "_" + string(argv[2]) + ".csv";
                lineCenter = (3.6) - (0.4*(r+1));
            }
            realSpaced.push_back((lineCenter)-(1.0*delta));
            realSpaced.push_back((lineCenter)+(0.0*delta));
            realSpaced.push_back((lineCenter)+(1.0*delta));

            imgSpaced.push_back(0.0);
            imgSpaced.push_back(0.0);
            vector<double> imgSpaced1 = kit.linspace((-1.0),(-0.6),sqrSideLen);
            vector<double> imgSpaced2 = kit.linspace((-0.6),(-0.2),sqrSideLen);
            vector<double> imgSpaced3 = kit.linspace((-0.2),(0.2),sqrSideLen);
            vector<double> imgSpaced4 = kit.linspace((0.2),(0.6),sqrSideLen);
            vector<double> imgSpaced5 = kit.linspace((0.6),(1.0),sqrSideLen);
            imgSpaced.insert(imgSpaced.end(),imgSpaced1.begin(),imgSpaced1.end());
            imgSpaced.insert(imgSpaced.end(),imgSpaced2.begin()+1,imgSpaced2.end());
            imgSpaced.insert(imgSpaced.end(),imgSpaced3.begin()+1,imgSpaced3.end());
            imgSpaced.insert(imgSpaced.end(),imgSpaced4.begin()+1,imgSpaced4.end());
            imgSpaced.insert(imgSpaced.end(),imgSpaced5.begin()+1,imgSpaced5.end());
            imgSpaced.push_back(0.0);
            imgSpaced.push_back(0.0);

            for (int i=0; i<sideLen; i++) { // populate the BFV with zeros to start with
                BFV.push_back(vector<int> (3,0));
            }
        }
        else {
            isVert = false;
            sideLen = sqrSideLen * numSqrs;
            if (isLeft) {
                outFileName = "leftHorizLine" + to_string(r-4) + "_" + string(argv[2]) + ".csv";
                
                realSpaced.push_back(0.0);
                realSpaced.push_back(0.0);
                vector<double> realSpaced1 = kit.linspace((-3.6),(-3.2),sqrSideLen);
                vector<double> realSpaced2 = kit.linspace((-3.2),(-2.8),sqrSideLen);
                vector<double> realSpaced3 = kit.linspace((-2.8),(-2.4),sqrSideLen);
                vector<double> realSpaced4 = kit.linspace((-2.4),(-2.0),sqrSideLen);
                vector<double> realSpaced5 = kit.linspace((-2.0),(-1.6),sqrSideLen);
                vector<double> realSpaced6 = kit.linspace((-1.6),(-1.2),sqrSideLen);
                realSpaced.insert(realSpaced.end(),realSpaced1.begin(),realSpaced1.end());
                realSpaced.insert(realSpaced.end(),realSpaced2.begin()+1,realSpaced2.end());
                realSpaced.insert(realSpaced.end(),realSpaced3.begin()+1,realSpaced3.end());
                realSpaced.insert(realSpaced.end(),realSpaced4.begin()+1,realSpaced4.end());
                realSpaced.insert(realSpaced.end(),realSpaced5.begin()+1,realSpaced5.end());
                realSpaced.insert(realSpaced.end(),realSpaced6.begin()+1,realSpaced6.end());
                realSpaced.push_back(0.0);
                realSpaced.push_back(0.0);
                realSpaced.push_back(0.0);
            }
            else {
                outFileName = "rightHorizLine" + to_string(r-4) + "_" + string(argv[2]) + ".csv";
                
                realSpaced.push_back(0.0);
                realSpaced.push_back(0.0);
                vector<double> realSpaced1 = kit.linspace((1.2),(1.6),sqrSideLen);
                vector<double> realSpaced2 = kit.linspace((1.6),(2.0),sqrSideLen);
                vector<double> realSpaced3 = kit.linspace((2.0),(2.4),sqrSideLen);
                vector<double> realSpaced4 = kit.linspace((2.4),(2.8),sqrSideLen);
                vector<double> realSpaced5 = kit.linspace((2.8),(3.2),sqrSideLen);
                vector<double> realSpaced6 = kit.linspace((3.2),(3.6),sqrSideLen);
                realSpaced.insert(realSpaced.end(),realSpaced1.begin(),realSpaced1.end());
                realSpaced.insert(realSpaced.end(),realSpaced2.begin()+1,realSpaced2.end());
                realSpaced.insert(realSpaced.end(),realSpaced3.begin()+1,realSpaced3.end());
                realSpaced.insert(realSpaced.end(),realSpaced4.begin()+1,realSpaced4.end());
                realSpaced.insert(realSpaced.end(),realSpaced5.begin()+1,realSpaced5.end());
                realSpaced.insert(realSpaced.end(),realSpaced6.begin()+1,realSpaced6.end());
                realSpaced.push_back(0.0);
                realSpaced.push_back(0.0);
                realSpaced.push_back(0.0);
            }
            lineCenter = (1.0) - (0.4*(r-4));

            imgSpaced.push_back((lineCenter)-(1.0*delta));
            imgSpaced.push_back((lineCenter)+(0.0*delta));
            imgSpaced.push_back((lineCenter)+(1.0*delta));
            for (int i=0; i<3; i++) { // populate the BFV with zeros to start with
                BFV.push_back(vector<int> (sideLen,0));
            }
        }
        /*
        if (r==0) {
            isVert = true;
            if (isLeft) {
                outFileName = "leftVertLine1_" + string(argv[2]) + "_" + string(argv[5]) + "_" + string(argv[4]) + ".csv";
                
                realSpaced.push_back((-2.8)-(1.0*delta));
                realSpaced.push_back((-2.8)+(0.0*delta));
                realSpaced.push_back((-2.8)+(1.0*delta));
            }
            else {
                outFileName = "rightVertLine1_" + string(argv[2]) + "_" + string(argv[5]) + "_" + string(argv[4]) + ".csv";
                
                realSpaced.push_back((2.8)-(1.0*delta));
                realSpaced.push_back((2.8)+(0.0*delta));
                realSpaced.push_back((2.8)+(1.0*delta));
            }
            imgSpaced.push_back(0.0);
            vector<double> imgSpaced1 = kit.linspace((-1.2),(-0.4),sqrSideLen);
            vector<double> imgSpaced2 = kit.linspace((-0.4),(0.4),sqrSideLen);
            vector<double> imgSpaced3 = kit.linspace((0.4),(1.2),sqrSideLen);
            imgSpaced.insert(imgSpaced.end(),imgSpaced1.begin(),imgSpaced1.end());
            imgSpaced.insert(imgSpaced.end(),imgSpaced2.begin()+1,imgSpaced2.end());
            imgSpaced.insert(imgSpaced.end(),imgSpaced3.begin()+1,imgSpaced3.end());
            imgSpaced.push_back(0.0);

            for (int i=0; i<sideLen; i++) { // populate the BFV with zeros to start with
                BFV.push_back(vector<int> (3,0));
            }
        }
        else if (r==1) {
            isVert = true;
            if (isLeft) {
                outFileName = "leftVertLine2_" + string(argv[2]) + "_" + string(argv[5]) + "_" + string(argv[4]) + ".csv";
                
                realSpaced.push_back((-2.0)-(1.0*delta));
                realSpaced.push_back((-2.0)+(0.0*delta));
                realSpaced.push_back((-2.0)+(1.0*delta));
            }
            else {
                outFileName = "rightVertLine2_" + string(argv[2]) + "_" + string(argv[5]) + "_" + string(argv[4]) + ".csv";
                
                realSpaced.push_back((2.0)-(1.0*delta));
                realSpaced.push_back((2.0)+(0.0*delta));
                realSpaced.push_back((2.0)+(1.0*delta));
            }
            
            imgSpaced.push_back(0.0);
            vector<double> imgSpaced1 = kit.linspace((-1.2),(-0.4),sqrSideLen);
            vector<double> imgSpaced2 = kit.linspace((-0.4),(0.4),sqrSideLen);
            vector<double> imgSpaced3 = kit.linspace((0.4),(1.2),sqrSideLen);
            imgSpaced.insert(imgSpaced.end(),imgSpaced1.begin(),imgSpaced1.end());
            imgSpaced.insert(imgSpaced.end(),imgSpaced2.begin()+1,imgSpaced2.end());
            imgSpaced.insert(imgSpaced.end(),imgSpaced3.begin()+1,imgSpaced3.end());
            imgSpaced.push_back(0.0);

            
            for (int i=0; i<sideLen; i++) { // populate the BFV with zeros to start with
                BFV.push_back(vector<int> (3,0));
            }
        }
        else if (r==2) {
            isVert = false;
            if (isLeft) {
                outFileName = "leftHorizLine1_" + string(argv[2]) + "_" + string(argv[5]) + "_" + string(argv[4]) + ".csv";
                
                realSpaced.push_back(0.0);
                vector<double> realSpaced1 = kit.linspace((-3.6),(-2.8),sqrSideLen);
                vector<double> realSpaced2 = kit.linspace((-2.8),(-2.0),sqrSideLen);
                vector<double> realSpaced3 = kit.linspace((-2.0),(-1.2),sqrSideLen);
                realSpaced.insert(realSpaced.end(),realSpaced1.begin(),realSpaced1.end());
                realSpaced.insert(realSpaced.end(),realSpaced2.begin()+1,realSpaced2.end());
                realSpaced.insert(realSpaced.end(),realSpaced3.begin()+1,realSpaced3.end());
                realSpaced.push_back(0.0);
            }
            else {
                outFileName = "rightHorizLine1_" + string(argv[2]) + "_" + string(argv[5]) + "_" + string(argv[4]) + ".csv";
                
                realSpaced.push_back(0.0);
                vector<double> realSpaced1 = kit.linspace((1.2),(2.0),sqrSideLen);
                vector<double> realSpaced2 = kit.linspace((2.0),(2.8),sqrSideLen);
                vector<double> realSpaced3 = kit.linspace((2.8),(3.6),sqrSideLen);
                realSpaced.insert(realSpaced.end(),realSpaced1.begin(),realSpaced1.end());
                realSpaced.insert(realSpaced.end(),realSpaced2.begin()+1,realSpaced2.end());
                realSpaced.insert(realSpaced.end(),realSpaced3.begin()+1,realSpaced3.end());
                realSpaced.push_back(0.0);
            }
            

            imgSpaced.push_back((0.4)-(1.0*delta));
            imgSpaced.push_back((0.4)+(0.0*delta));
            imgSpaced.push_back((0.4)+(1.0*delta));
            for (int i=0; i<3; i++) { // populate the BFV with zeros to start with
                BFV.push_back(vector<int> (sideLen,0));
            }
        }
        else {
            cout << "TOO MANY LINES !!\n";
            return -1;
        }
        */
        for (int b=0; b<numIters; b++) {
            comeca = std::chrono::high_resolution_clock::now();
            if (b==0 || (b+1)%5==0 || (b+1)==numIters) {
                std::cout << "Starting iteration # " << to_string(b+1) << "/" << to_string(numIters) <<
                    " at T = " << std::chrono::duration_cast<std::chrono::minutes>(comeca-comeca0).count()
                    << " mins..." << endl;
            }

            int startPoly = b*coeffSize;
            
            vectPolyToUse.clear();
            axesBinVect.clear();

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
            AxesEval axeseval(vectPolyToUse,kit,realSpaced,imgSpaced,isVert,polySize,sideLen,numSamples,numPolys,largeNum);
            

            /**
             * @brief creates the local min vector
             *
             * @return 2D VECTOR and ONLY GRIDLINES of the image
             */
            
            if (isVert) {
                axesBinVect = axeseval.vertGetBinCount();
            }
            else {
                axesBinVect = axeseval.horiGetBinCount();
            }
            for (long unsigned int g=0; g<BFV.size(); g++) {
                for (long unsigned int h=0; h<BFV[0].size(); h++) {
                    BFV[g][h] += axesBinVect[g][h];
                }
            }
            
        }
        
        string fileName = "../output/" + outFileName; // for general use
        // string fileName = "../../../.." + scratchFile + "/" + outFileName; // for use with scratch on hpc
        ofstream myFile;
        myFile.open(fileName);
        for (long unsigned int g=0; g<BFV.size(); g++) {
            myFile << BFV[g][0];
            for (long unsigned int h=1; h<BFV[0].size(); h++) {
                myFile << ',' << BFV[g][h];
            }
            myFile << '\n';
        }
        myFile.close();
        fimBatch = std::chrono::high_resolution_clock::now();
        std::cout << "End of line # " << to_string(r+1) << "/" << to_string(7) <<
                " --- Duration = " << std::chrono::duration_cast<std::chrono::minutes>(fimBatch-comecaBatch).count()
                << " mins...\n" << endl;
    }    
    auto fim0 = std::chrono::high_resolution_clock::now();
    auto duration5 = std::chrono::duration_cast<std::chrono::minutes>(fim0-comeca0);
    std::cout << "TOTAL RUNTIME: " << duration5.count() << " mins\n";
    return 0;
}