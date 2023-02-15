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
 * @brief This is the PREFERRED DRIVER for the program. It splits an image into a grid
 * and analyzes each of the squares one at a time. This allows for high res pics with a 
 * reasonably large number of polynomials analyzed per batch.
 * 
 * 
 *
 * @return int
 */

int main(int argc, char *argv[])
{
    string bNumberStr;
    if (argc != 2) {
        cout << "Specify which coeff file to use. Enter a number 1-4.\n";
        getline(cin,bNumberStr);
        cout << '\n';
    }
    else {
        bNumberStr = argv[1];
    }
    auto comeca0 = std::chrono::high_resolution_clock::now();
    int bNumber = stoi(bNumberStr);

    int sideLen = 1521; // side length (in pixels) of the resulting image
    int polySize = 13; // degree of the polynomial to be used
    int numSamples = 1; // how many threads to use
    int numPolys = 1; // NEEDS TO BE LARGER THAN NUMSAMPLES !!!!!
    int coeffSize = (polySize *numPolys); // num coeffs to load for real/img
    bool polyIsArr = false;
    // string fileNum = "507_8000_matrix_zoomOut";
    string fileNum = "1521_144_matrix_zoomInMidBasin4";

    Tools kit;

    vector<double> realSpaced;
    vector<double> imgSpaced;
    double largeNum = pow(10,300); // upper bound for polyEval vals (above ~ infinity)

    auto dateTime = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    cout << ctime(&dateTime) << '\n'; 
    cout << "Using coeff file # " << bNumberStr << ". Starting with " << numSamples << " samples...\n";
    auto comeca = std::chrono::high_resolution_clock::now();
    auto fim = std::chrono::high_resolution_clock::now();
    auto duration4 = std::chrono::duration_cast<std::chrono::minutes>(fim-comeca);
    auto comecaBatch = std::chrono::high_resolution_clock::now();
    auto fimBatch = std::chrono::high_resolution_clock::now();
    auto durationBatch = std::chrono::duration_cast<std::chrono::minutes>(fimBatch-comecaBatch);
    
    for (int b=1; b<2; b++) {
        bNumber = b;
        for (int q=0; q<1; q++) { // each row of batches in the image. Change back to q=0, q<6
            // if ((q==1) || (q==2)) { // remove to go back to whole rows
            //     cout << "ROW " << to_string(q) << " ALREADY WRITTEN\n";
            //     continue;
            // }
            comecaBatch = std::chrono::high_resolution_clock::now();
            imgSpaced.clear();
            // imgSpaced = kit.linspace((-0.4+(0.8*q)),(-0.4+(0.8*(q+1))),sideLen);
            // imgSpaced = kit.linspace(1.98,2.02,sideLen);
            imgSpaced = kit.linspace(-6.0,6.0,sideLen);
            cout << "Starting row # " << b << " at T = " << 
                std::chrono::duration_cast<std::chrono::minutes>(comecaBatch-comeca0).count()
                << " mins...\n";
            
            for (int r=144; r<288; r++) { // change back to r=0, r<15
                // if ((q==0) && (r==3)) { // remove to go back to whole rows
                //     cout << "FILE " << to_string(r+(15*q)) << " ALREADY WRITTEN\n";
                //     continue;
                // }
                // else if ((q==0) && (r==4)) { // remove to go back to whole rows
                //     cout << "FILE " << to_string(r+(15*q)) << " ALREADY WRITTEN\n";
                //     continue;
                // }
                // else if ((q==0) && (r==7)) { // remove to go back to whole rows
                //     cout << "FILE " << to_string(r+(15*q)) << " ALREADY WRITTEN\n";
                //     continue;
                // }
                // else if ((q==3) && (r==7)) { // remove to go back to whole rows
                //     cout << "FILE " << to_string(r+(15*q)) << " ALREADY WRITTEN\n";
                //     continue;
                // }
                // if ((((r-7)*(r-7))+(q*q)) >= 50) {
                //     string fileName = "../output/mat12x12_" + fileNum + to_string(r+(15*q)) + "_Tot" + ".csv";
                //     kit.writeBlankFile(fileName, sideLen);
                //     cout << "FILE " << to_string(r+(15*q)) << " CONTAINS NO ROOTS\n";
                //     continue;
                // }
                // else if ((q==3&&r==1) || (q==3&&r==13)){
                //     string fileName = "../output/mat12x12_" + fileNum + to_string(r+(15*q)) + "_Tot" + ".csv";
                //     kit.writeBlankFile(fileName, sideLen);
                //     cout << "FILE " << to_string(r+(15*q)) << " CONTAINS NO ROOTS\n";
                //     continue;
                // }
                comeca = std::chrono::high_resolution_clock::now();

                int startPoly = r;

                realSpaced.clear();
                // realSpaced = kit.linspace((-6.0 +(0.8*r)),(-6.0 +(0.8*(r+1))),sideLen);
                // realSpaced = kit.linspace(-0.02,0.02,sideLen);
                realSpaced = kit.linspace(-6.0,6.0,sideLen);

                // if (r%4==0) {
                //     realSpaced = kit.linspace(1.905,1.945,sideLen);
                // }
                // else if (r%4==1) {
                //     realSpaced = kit.linspace(1.945,1.985,sideLen);
                // }
                // else if (r%4==2) {
                //     realSpaced = kit.linspace(1.985,2.025,sideLen);
                // }
                // else if (r%4==3) {
                //     realSpaced = kit.linspace(2.025,2.065,sideLen);
                // }
                
                // if (r<4) {
                //     imgSpaced = kit.linspace(0.02,0.06,sideLen);
                // }
                // else if (r<8) {
                //     imgSpaced = kit.linspace(-0.02,0.02,sideLen);
                // }
                // else if (r<12) {
                //     imgSpaced = kit.linspace(-0.06,-0.02,sideLen);
                // }

                // int startPoly = r *numPolys *polySize; //uncomment for polynomials instead of matrices
                // string fileName = "../output/deg24_" + fileNum + to_string(r) + ".csv";
                // vector<vector<unsigned short int>> quartBinVectTot ((sideLen/2), vector<unsigned short int> ((sideLen/2),0));
                for (int i=0; i<1; i++) {
                    // auto start = std::chrono::high_resolution_clock::now();
                    string fileName = "../output/mat12x12_" + fileNum + to_string(r+(15*q)) + "_" + to_string(i+(2*(bNumber-1))) + ".csv";
                    int offset = i *coeffSize; // offsets the coeff's indices. Must be smaller than 24 !!!!
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
                    // vector<double> realP; // real part of each coeff
                    // vector<double> imgP; // complex part of each coeff
                    // ifstream coeffFile;
                    // coeffFile.open("../src/testCoeffsOnesBigRandom.csv"); // file holding old coeffs
                    // string coeff;
                    // double coeffDoub;
                    // int lineNum = 0;
                    // if (coeffFile.is_open()) {
                    //     while (coeffFile) {
                    //         getline(coeffFile,coeff);
                    //         coeffDoub = stod(coeff);
                    //         if ((lineNum>(offset+startPoly)) && (lineNum < (coeffSize+offset+startPoly))) {
                    //             realP.push_back(coeffDoub);
                    //             imgP.push_back(0.0); // uncomment to keep real coeffs
                    //         }
                    //         // else if (lineNum < (coeffSize*2)) {
                    //         //     // imgP.push_back(coeffDoub); // uncomment to add complex coeffs
                    //         //     imgP.push_back(0.0); // uncomment to keep real coeffs
                    //         // }
                    //         lineNum ++;
                    //     }
                    // }

                    /**
                     * THIS IS ANOTHER ALTERNATIVE TO THE ABOVE
                     * IT READS PREVIOUSLY SAVED COEFFS FROM !!!RANDOM MATRICES/CHARACTERISTIC POLYNOMIALS!!!
                    */
                    vector<double> realP; // real part of each coeff
                    vector<double> imgP; // complex part of each coeff
                    ifstream coeffFile;
                    string name = "../src/charPolyCoeffsBig" + to_string(bNumber) + ".csv";
                    coeffFile.open(name); // file holding old coeffs
                    string coeff;
                    double coeffDoub;
                    int lineNum = 0;
                    if (coeffFile.is_open()) {
                        while (coeffFile) {
                            getline(coeffFile,coeff);
                            if (coeff != "") {
                                coeffDoub = stod(coeff);
                                if ((lineNum>=(offset+startPoly)) && (lineNum < (coeffSize+offset+startPoly))) {
                                    realP.push_back(coeffDoub);
                                    imgP.push_back(0.0); // uncomment to keep real coeffs
                                }
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

                    // std::cout << "Done here. Number of Samples: " << numSamples << "\nFile number: " << offset << endl;

                    // vector<double> realSpaced = kit.linspace(-1.6,1.6,sideLen);
                    // vector<double> imgSpaced = kit.linspace(-1.6,1.6,sideLen);
                    // vector<double> realSpaced = kit.linspace(-1.55,1.55,sideLen);
                    // vector<double> imgSpaced = kit.linspace(-1.55,1.55,sideLen);
                    // vector<double> realSpaced = kit.linspace(0.55,0.65,sideLen);
                    // vector<double> imgSpaced = kit.linspace(0.25,0.35,sideLen);
                    PolyEval polyeval(vectPolyToUse,kit,realSpaced,imgSpaced,polySize,sideLen,numSamples,numPolys,largeNum);
                    // HalfPolyEval halfpolyeval(vectPolyToUse,kit,realSpaced,imgSpaced,polySize,sideLen,numSamples,numPolys,largeNum);
                    // QuartPolyEval quartpolyeval(vectPolyToUse,kit,realSpaced,imgSpaced,polySize,sideLen,numSamples,numPolys,largeNum);
                    // FlatHalfPoly flathalfpoly(vectPolyToUse,kit,realSpaced,imgSpaced,polySize,sideLen,numSamples,largeNum);
                    // AxesEval axeseval(vectPolyToUse,kit,realSpaced,imgSpaced,polySize,sideLen,numSamples,numPolys,largeNum);
                    
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
                     * @brief creates the local min vector
                     *
                     * @return 2D VECTOR and ONLY QUARTER of the image
                     * @note MOST EFFICIENT AND CAN ONLY BE USED WITH COMPLEX CONJ. ROOT THRM. (CCRT)
                     * 
                     * @note MUST BE USED WITH LARGE SAMPLES BC INDIVIDUAL POLYS ARE NOT SYMMETRIC
                     *       OVER THE IMAGINARY AXIS !!!!!!!!!!!!!!
                     */
                    // vector<vector<unsigned short int>> quartBinVectTemp = quartpolyeval.getBinCount3();
                    // // auto start3 = std::chrono::high_resolution_clock::now();
                    // for (int g=0; g<(sideLen/2); g++) {
                    //    for (int h=0; h<(sideLen/2); h++) {
                    //        quartBinVectTot[g][h] += quartBinVectTemp[g][h];
                    //    }
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
                    // vector<vector<unsigned short int>> axesBinVect = axeseval.getBinCount();
                    // auto start3 = std::chrono::high_resolution_clock::now();
                    // for (int g=0; g<(sideLen-1); g++) {
                    // myFile << axesBinVect[g][0];
                    // for (int h=1; h<(sideLen-1); h++) {
                    //     myFile << ',' << axesBinVect[g][h];
                    // }
                    // myFile << "\n";
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
                    ofstream myFile;
                    myFile.open(fileName);
                    vector<vector<unsigned short int>> wholeBinVect = polyeval.getBinCount3();
                    // auto start3 = std::chrono::high_resolution_clock::now();
                    for (int g=0; g<sideLen; g++) {
                        myFile << wholeBinVect[g][0];
                        for (int h=1; h<sideLen; h++) {
                            myFile << ',' << wholeBinVect[g][h];
                        }
                        myFile << '\n';
                    }
                    myFile.close();
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

                    

                    // auto end = std::chrono::high_resolution_clock::now();
                    // auto duration = std::chrono::duration_cast<std::chrono::seconds>(start3-start2);
                    // auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(end-start3);
                    // auto duration3 = std::chrono::duration_cast<std::chrono::seconds>(end-start);

                    // std::cout << "Sample #" << offset << " runtime: " << duration3.count() << "s\n";
                }

                /**
                 * UNCOMMENT THE FILE WRITING SECTION BELOW TO UTILIZE QUARTPOLYEVAL
                 * 
                 */
                // ofstream myFile;
                // myFile.open(fileName);

                // for (int g=0; g<(sideLen/2); g++) {
                //    myFile << quartBinVectTot[g][0];
                //    for (int h=1; h<(sideLen/2); h++) {
                //        myFile << ',' << quartBinVectTot[g][h];
                //    }
                //    for (int i=((sideLen/2)-2); i>=0; i--) {
                //        myFile << ',' << quartBinVectTot[g][i];
                //    }
                //    myFile << "\n";
                // }
                // for (int g=((sideLen/2)-2); g>=0; g--) {
                //    myFile << quartBinVectTot[g][0];
                //    for (int h=1; h<(sideLen/2); h++) {
                //        myFile << ',' << quartBinVectTot[g][h];
                //    }
                //    for (int i=((sideLen/2)-2); i>=0; i--) {
                //        myFile << ',' << quartBinVectTot[g][i];
                //    }
                //    myFile << "\n";
                // }

                // myFile.close();

                fim = std::chrono::high_resolution_clock::now();
                duration4 = std::chrono::duration_cast<std::chrono::minutes>(fim-comeca);
                std::cout << "FILE " << to_string(r+(15*q)) << " RUNTIME: " << duration4.count() << " mins\n";
            }
            fimBatch = std::chrono::high_resolution_clock::now();
            durationBatch = std::chrono::duration_cast<std::chrono::minutes>(fimBatch-comecaBatch);
            std::cout << "ROW " << q << "  ---  RUNTIME: " << durationBatch.count() << " mins\n";
        }
    }
    auto fim0 = std::chrono::high_resolution_clock::now();
    auto duration5 = std::chrono::duration_cast<std::chrono::minutes>(fim0-comeca0);
    std::cout << "TOTAL RUNTIME: " << duration5.count() << " mins\n";
    return 0;
}
