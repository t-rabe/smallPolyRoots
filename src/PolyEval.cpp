#include "../headers/PolyEval.hpp"


PolyEval::PolyEval(vector<complex<double>> vectPolyToUse_, Tools kit_, vector<double> realSpaced_, vector<double> imgSpaced_,
                                                int polySize_, int imgSize_, int numSamples_, int numOfPoly_, double bigBoy_) {
    vectPolyToUse = vectPolyToUse_;
    kit = kit_;
    realSpaced = realSpaced_;
    imgSpaced = imgSpaced_;
    polySize = polySize_;
    imgSize = imgSize_;
    numSamples = numSamples_;
    numOfPoly = numOfPoly_;
    bigBoy = bigBoy_;
    realSpan = floor(realSpaced_.size() /numSamples_);
    polySpan = max(int(floor(numOfPoly /numSamples)), 1);

    multPixValArr = new float**[imgSize];
    multBinValArr = new int**[imgSize];
    binCountArr3 = new int*[imgSize];
    
    // creates a 2d vector with (real x img) indices

    for (int f=0; f<imgSize; f++) {
        vector<vector<float>> interMed00(imgSize, vector<float> (numOfPoly,0.0));
        multPixValVect.push_back(interMed00);
    }

    for (int f=0; f<imgSize; f++) {
        vector<vector<int>> interMed01(imgSize, vector<int> (numOfPoly,0));
        multBinValVect.push_back(interMed01);
    }

    for (int h=0; h<imgSize; h++) {
        vector<unsigned short int> interMed02(imgSize,0);
        binCountVect3.push_back(interMed02);
    }

    /*
    for (int g=0; g<imgSize; g++) {
        vector<int> interMed0(imgSize,0);
        binCountVect.push_back(interMed0);
    }

    for (int h=0; h<imgSize; h++) {
        vector<int> interMed1(imgSize,0);
        binCountVect2.push_back(interMed1);
    }

    for (int h2=0; h2<imgSize; h2++) {
        vector<int> interMed12(imgSize,0);
        maxBinCountVect.push_back(interMed12);
    }

    for (int i=0; i<imgSize; i++) {
        vector<double> interMed2(imgSize,0.0);
        pixelValVect.push_back(interMed2);
    }

    for (int j=0; j<imgSize; j++) {
        vector<int> interMed3(imgSize,0);
        colMinVect.push_back(interMed3);
    }

    for (int k=0; k<imgSize; k++) {
        vector<int> interMed4(imgSize,0);
        rowMinVect.push_back(interMed4);
    }

    for (int j2=0; j2<imgSize; j2++) {
        vector<int> interMed32(imgSize,0);
        colMaxVect.push_back(interMed32);
    }

    for (int k2=0; k2<imgSize; k2++) {
        vector<int> interMed42(imgSize,0);
        rowMaxVect.push_back(interMed42);
    }

    for (int k3=0; k3<imgSize; k3++) {
        vector<int> interMed43(imgSize,0);
        minPeaksVect.push_back(interMed43);
    }
    */

}

// creates a matrix of the vals of the polynomial at each pixel
void PolyEval::createMat(int startReal, int endReal) {
    double realVal = 0.0;
    double imgVal = 0.0;
    // double rad = 0.0;
    for (int k=startReal; k<endReal; k++) {
        for (int m=0; m<imgSize; m++) {
            realVal = realSpaced[k];
            imgVal = imgSpaced[m];
            // rad = sqrt(((realVal*realVal) + (imgVal*imgVal)));
            // if (rad < 1.05) {
            //     pixelValVect[k][m] = abs(kit.horner7(vectPolyToUse, polySize, {realVal,imgVal}));
            // }
            pixelValVect[k][m] = abs(kit.horner7(vectPolyToUse, polySize, {realVal,imgVal}, bigBoy));
        }
    }
}

// creates a matrix of the vals of the polynomial at each pixel for MULTIPLE POLYS
void PolyEval::createMat2(int startReal, int endReal) {
    double realVal = 0.0;
    double imgVal = 0.0;
    vector<float> tempVect;
    // double rad = 0.0;
    for (int k=startReal; k<endReal; k++) {
        for (int m=0; m<imgSize; m++) {
            realVal = realSpaced[k];
            imgVal = imgSpaced[m];
            tempVect = kit.horner9(vectPolyToUse, numOfPoly, polySize, {realVal,imgVal}, bigBoy);
            multPixValVect[k][m].swap(tempVect);
        }
    }
}

// creates a matrix of vals of the polynomial at each pixel for MULTIPLE POLYS AS ARRAY
void PolyEval::createMat3(int startReal, int endReal) {
    double realVal = 0.0;
    double imgVal = 0.0;
    for (int k=startReal; k<endReal; k++) {
        multPixValArr[k] = new float*[imgSize];
        for (int m=0; m<imgSize; m++) {
            realVal = realSpaced[k];
            imgVal = imgSpaced[m];
            multPixValArr[k][m] = new float[numOfPoly];
            multPixValArr[k][m] = kit.horner8(vectPolyToUse, numOfPoly, polySize, {realVal,imgVal}, bigBoy);
        }
    }
    // delete tempArr;
    // double realVal = 0.0;
    // double imgVal = 0.0;
    // float **tempArr = new float*[imgSize];
    // // double rad = 0.0;
    // for (int k=startReal; k<endReal; k++) {
    //     // float **tempArr = new float*[imgSize];
    //     for (int m=0; m<imgSize; m++) {
    //         realVal = realSpaced[k];
    //         imgVal = imgSpaced[m];
    //         tempArr[m] = kit.horner8(vectPolyToUse, numOfPoly, polySize, {realVal,imgVal}, bigBoy);
    //     }
    //     multPixValArr[k] = tempArr;
    //     // delete tempArr;
    // }
    // // delete tempArr;
}

// finds all local mins for a given poly within the vect of all polys
// NOTE: does not check edges of the image
void PolyEval::findAllMins(int startPoly, int endPoly) {
    float curr;
    float top;
    float left;
    float right;
    float bott;
    float tota;
    for (int k=1; k<(imgSize-1); k++) {
        for (int m=1; m<(imgSize-1); m++) {
            for (int n=startPoly; n<endPoly; n++) {
                curr = multPixValVect[k][m][n];
                top = multPixValVect[k-1][m][n];
                left = multPixValVect[k][m-1][n];
                right = multPixValVect[k][m+1][n];
                bott = multPixValVect[k+1][m][n];
                tota = curr+top+left+right+bott;
                // cout << curr << ' ' << top << ' ' << left << ' ' << right << ' ' << bott << '\n';
                if ((tota>0) && (curr<=top) && (curr<=left) && (curr<=right) && (curr<=bott)) {
                    multBinValVect[k][m][n] ++;
                }
            }
        }
    }
}

// finds all local mins for a given poly within the ARRAY of all polys as ARRAY
// NOTE: does not check edges of the image THIS LEAVES EDGE INDICES WITHOUT A VALUE
void PolyEval::findAllMins2(int startPoly, int endPoly) {
    float curr;
    float top;
    float left;
    float right;
    float bott;
    float tota;
    // int **tempArr2 = new int*[imgSize]; // second dimension
    // int *tempArr3 = new int[numOfPoly]; // third dimension, represents individual pixels

    multBinValArr[0] = new int*[imgSize];
    multBinValArr[imgSize -1] = new int*[imgSize];
    for (int y=0; y<imgSize; y++) {
        multBinValArr[0][y] = new int[numOfPoly];
        multBinValArr[imgSize -1][y] = new int[numOfPoly];
        for (int n=startPoly; n<endPoly; n++) {
            multBinValArr[0][y][n] = 0;
            multBinValArr[imgSize -1][y][n] = 0;
        }
    }

    for (int k=1; k<(imgSize-1); k++) {
        multBinValArr[k] = new int*[imgSize]; // second dimension
        // cout << "check" << endl;
        for (int m=1; m<(imgSize-1); m++) {
            multBinValArr[k][m] = new int[numOfPoly]; // third dimension, represents individual pixels
            // cout << "check2" << endl;
            for (int n=startPoly; n<endPoly; n++) {
                // cout << "k = " << k << " m = " << m << " n = " << n << '\n';
                curr = multPixValArr[k][m][n];
                top = multPixValArr[k-1][m][n];
                left = multPixValArr[k][m-1][n];
                right = multPixValArr[k][m+1][n];
                bott = multPixValArr[k+1][m][n];
                tota = curr+top+left+right+bott;
                // cout << "k = " << k << " m = " << m << " n = " << n << '\n';
                // cout << curr << ' ' << top << ' ' << left << ' ' << right << ' ' << bott << '\n';
                // cout << "check3" << endl;
                if ((tota>0) && (curr<=top) && (curr<=left) && (curr<=right) && (curr<=bott)) {
                    // cout << "1__";
                    multBinValArr[k][m][n] = 1;
                    // cout << "2\n";
                }
                else {
                    // cout << k << "_" << m << "_" << n << "__";
                    multBinValArr[k][m][n] = 0;
                    // cout << "2\n";
                }
            }
            // tempArr2[m] = tempArr3;
            // cout << "k = " << k << " m = " << m << '\n';
            // delete tempArr3;
        }
        // multBinValArr[k] = tempArr2;
        // cout << "check5\n";
        // delete tempArr2;
    }
    // delete tempArr2;
    // delete tempArr3;
    // float curr;
    // float top;
    // float left;
    // float right;
    // float bott;
    // float tota;
    // // int **tempArr2 = new int*[imgSize]; // second dimension
    // // int *tempArr3 = new int[numOfPoly]; // third dimension, represents individual pixels

    // for (int k=1; k<(imgSize-1); k++) {
    //     int **tempArr2 = new int*[imgSize]; // second dimension
    //     // cout << "check" << endl;
    //     for (int m=1; m<(imgSize-1); m++) {
    //         int *tempArr3 = new int[numOfPoly]; // third dimension, represents individual pixels
    //         // cout << "check2" << endl;
    //         for (int n=startPoly; n<endPoly; n++) {
    //             curr = multPixValArr[k][m][n];
    //             top = multPixValArr[k-1][m][n];
    //             left = multPixValArr[k][m-1][n];
    //             right = multPixValArr[k][m+1][n];
    //             bott = multPixValArr[k+1][m][n];
    //             tota = curr+top+left+right+bott;
    //             // cout << curr << ' ' << top << ' ' << left << ' ' << right << ' ' << bott << '\n';
    //             // cout << "check3" << endl;
    //             if ((tota>0) && (curr<=top) && (curr<=left) && (curr<=right) && (curr<=bott)) {
    //                 tempArr3[n] = 1;
    //             }
    //             else {
    //                 tempArr3[n] = 0;
    //             }
    //         }
    //         tempArr2[m] = tempArr3;
    //         // cout << "check4" << endl;
    //         // delete tempArr3;
    //     }
    //     multBinValArr[k] = tempArr2;
    //     // cout << "check5\n";
    //     // delete tempArr2;
    // }
    // // delete tempArr2;
    // // delete tempArr3;
}

// evals individual pixels and updates their bin count. cannot factor out the found roots
// counts a pixel as a root if the val of the poly is "close to zero" (within a certain tolerance)
void PolyEval::evalPixel(int startReal, int endReal) {
    // int binNum = 0;
    int rootNum = 0;
    // complex<double> currLoc = {0.,0.};
    // vector<complex<double>> tempCoeffs;
    double tolerance = 0.0;
    double realVal = 0.0;
    double imgVal = 0.0;
    double firstCoeff = abs(vectPolyToUse[0]);
    double secCoeff = abs(vectPolyToUse[1]);
    // complex<double> secCoeff = abs(vectPolyToUse[1]);    
    // complex<double> firstCoeff = abs(vectPolyToUse[0]-vectPolyToUse[1]);
    // complex<double> firstCoeff = vectPolyToUse[0]+vectPolyToUse[1];
    complex<double> loc = {0.0,0.0};

    for (int k=startReal; k<endReal; k++) {
        for (int m=0; m<imgSize; m++) {
            // currLoc = {realSpaced[k],imgSpaced[m]};
            // tempCoeffs = kit.horner7(vectPolyToUse, polySize, {realSpaced[k],imgSpaced[m]});
            realVal = realSpaced[k];
            imgVal = imgSpaced[m];
            loc = {realVal,imgVal};
            tolerance = (firstCoeff *pow((abs(loc)),4000)) - (secCoeff *pow((abs(loc)),3999));
            // tolerance = abs((firstCoeff)*((realVal*realVal)+(imgVal*imgVal)));
            // tolerance = k *1.0;
            if (abs(kit.horner7(vectPolyToUse, polySize, {realVal,imgVal}, bigBoy)) < tolerance) {
                binCountVect[k][m] ++;
                rootNum++;
            }
        }
    }

    // for (int k=startReal; k<endReal; k++) {
    //     for (int m=0; m<imgSize; m++) {
    //         // currLoc = {realSpaced[k],imgSpaced[m]};
    //         tempCoeffs = kit.horner6(vectPolyToUse, polySize, {realSpaced[k],imgSpaced[m]});
    //         if (abs(tempCoeffs.back()) < 0.1) {
    //             binCountVect[k][m] ++;
    //             rootNum++;
    //         }
    //     }
    // }
    cout << "Num roots: " << rootNum << endl;
    // for a 2d vector which is supposed to be flattened later. Never implemented
    // for (int k=0; k<realSpan; k++) {
    //     for (int m=0; m<imgSize; m++) {
    //         binNum = (k*imgSize) +m;
    //         currLoc = {realSpaced[k],imgSpaced[m]};
    //         tempCoeffs = kit.horner6(vectPolyToUse, polySize, currLoc);
    //
    //         if (abs(tempCoeffs.back()) < 0.1) {
    //             binCountVect[index][binNum] ++;
    //             // tempCoeffs.pop_back();
    //             // vectPolyToUse = tempCoeffs;
    //             rootNum++;
    //         }
    //     }
    // }
}

// finds the local mins in every column to be compared with data from findRowMin(...)
void PolyEval::findColMin(int startReal, int endReal) {
    double prev;
    double curr;
    double next;
    for (int k=startReal; k<endReal; k++) {
        prev = pixelValVect[k][0];
        curr = pixelValVect[k][0];
        for (int m=1; m<imgSize; m++) {
            next = pixelValVect[k][m];
            if ((curr <= prev) && (curr <= next)) {
                colMinVect[k][m-1] ++;
            }
            prev = curr;
            curr = next;
        }
        if (curr <= prev) {
            colMinVect[k][imgSize-1] ++;
        }
    }
}

// finds the local mins in every row to be compared with data from findColMin(...)
void PolyEval::findRowMin(int startImg, int endImg) {
    double prev2;
    double curr2;
    double next2;
    for (int k=startImg; k<endImg; k++) {
        prev2 = pixelValVect[0][k];
        curr2 = pixelValVect[0][k];
        for (int m=1; m<imgSize; m++) {
            next2 = pixelValVect[m][k];
            if ((curr2 <= prev2) && (curr2 <= next2)) {
                rowMinVect[m-1][k] ++;
            }
            prev2 = curr2;
            curr2 = next2;
        }
        if (curr2 <= prev2) {
            rowMinVect[imgSize-1][k] ++;
        }
    }
}

// finds the local max in every column to be compared with data from findRowMax(...)
void PolyEval::findColMax(int startReal, int endReal) {
    double prev;
    double curr;
    double next;
    for (int k=startReal; k<endReal; k++) {
        prev = pixelValVect[k][0];
        curr = pixelValVect[k][0];
        for (int m=1; m<imgSize; m++) {
            next = pixelValVect[k][m];
            if ((curr >= prev) && (curr >= next)) {
                colMaxVect[k][m-1] --;
            }
            prev = curr;
            curr = next;
        }
        if (curr >= prev) {
            colMaxVect[k][imgSize-1] --;
        }
    }
}

// finds the local max in every row to be compared with data from findColMax(...)
void PolyEval::findRowMax(int startImg, int endImg) {
    double prev2;
    double curr2;
    double next2;
    for (int k=startImg; k<endImg; k++) {
        prev2 = pixelValVect[0][k];
        curr2 = pixelValVect[0][k];
        for (int m=1; m<imgSize; m++) {
            next2 = pixelValVect[m][k];
            if ((curr2 >= prev2) && (curr2 >= next2)) {
                rowMaxVect[m-1][k] --;
            }
            prev2 = curr2;
            curr2 = next2;
        }
        if (curr2 >= prev2) {
            rowMaxVect[imgSize-1][k] --;
        }
    }
}

// finds the local max in every row to be compared with data from findColMax(...)
void PolyEval::findMinPeaks(int startReal, int endReal) {
    int curr;
    int totNear;
    for (int k=(startReal+1); k<(endReal-1); k++) {
        for (int m=1; m<(imgSize-1); m++) {
            curr = binCountVect2[k][m];
            if (curr == 2) {
                minPeaksVect[k][m] ++;
            }
            else if (curr == 1) {
                // the sum of all neighboring pixel values
                totNear = binCountVect2[k-1][m-1]+binCountVect2[k-1][m]+binCountVect2[k-1][m+1]
                        + binCountVect2[k][m-1]+binCountVect2[k][m+1]
                        + binCountVect2[k+1][m-1]+binCountVect2[k+1][m]+binCountVect2[k+1][m+1];
                if (totNear == 0) {
                    minPeaksVect[k][m] ++;
                }
            }
        }
    }
}

// safely multithreads each sample to find roots within a certain tolerance
void PolyEval::threadSafe_Sample() {
    vector<thread> tasks;
    int startReal_ = 0;
    int endReal_ = 0;
    for (int j=0; j<numSamples; j++) {
        startReal_ = endReal_;
        endReal_ += realSpan;
        tasks.push_back(thread(&PolyEval::evalPixel, this, startReal_, endReal_));
    }
    for (unsigned int f=0; f<tasks.size(); f++) {
        tasks[f].join();
    }
}

// safely multithreads to create a matrix of pixel vals
void PolyEval::threadSafe_Sample2() {
    vector<thread> tasks2;
    int startReal_ = 0;
    int endReal_ = 0;
    for (int j=0; j<numSamples; j++) {
        startReal_ = endReal_;
        endReal_ += realSpan;
        tasks2.push_back(thread(&PolyEval::createMat, this, startReal_, endReal_));
    }
    for (unsigned int f=0; f<tasks2.size(); f++) {
        tasks2[f].join();
    }
}

// safely multithreads to find local mins
void PolyEval::threadSafe_Sample3() {
    vector<thread> tasks3;
    int start_ = 0;
    int end_ = 0;
    for (int j=0; j<numSamples; j++) {
        start_ = end_;
        end_ += realSpan;
        tasks3.push_back(thread(&PolyEval::findColMin, this, start_, end_));
        tasks3.push_back(thread(&PolyEval::findRowMin, this, start_, end_));
    }
    for (unsigned int f=0; f<tasks3.size(); f++) {
        tasks3[f].join();
    }
}

// safely multithreads to find local max
void PolyEval::threadSafe_Sample4() {
    vector<thread> tasks4;
    int start_ = 0;
    int end_ = 0;
    for (int j=0; j<numSamples; j++) {
        start_ = end_;
        end_ += realSpan;
        tasks4.push_back(thread(&PolyEval::findColMax, this, start_, end_));
        tasks4.push_back(thread(&PolyEval::findRowMax, this, start_, end_));
    }
    for (unsigned int f=0; f<tasks4.size(); f++) {
        tasks4[f].join();
    }
}

// safely multithreads to find peaks of the min val vector
void PolyEval::threadSafe_Sample5() {
    vector<thread> tasks5;
    int start_ = 0;
    int end_ = realSpan;
    for (int j=0; j<numSamples; j++) {
        tasks5.push_back(thread(&PolyEval::findMinPeaks, this, start_, end_));
        start_ = end_ -2;
        end_ += realSpan;
    }
    for (unsigned int f=0; f<tasks5.size(); f++) {
        tasks5[f].join();
    }
}

// safely multithreads to create a matrix of pixel vals for MULTIPLE POLYNOMIALS
void PolyEval::threadSafe_Sample6() {
    vector<thread> tasks6;
    int startReal_ = 0;
    int endReal_ = 0;
    for (int j=0; j<numSamples; j++) {
        startReal_ = endReal_;
        endReal_ += realSpan;
        tasks6.push_back(thread(&PolyEval::createMat2, this, startReal_, endReal_));
    }
    for (unsigned int f=0; f<tasks6.size(); f++) {
        tasks6[f].join();
    }
}

// safely multithreads to find local mins for MULTIPLE POLYNOMIALS
void PolyEval::threadSafe_Sample7() {
    // cout << polySpan << '\n';
    vector<thread> tasks7;
    int start_ = 0;
    int end_ = 0;
    for (int j=0; j<numSamples; j++) {
        start_ = end_;
        end_ += polySpan;
        tasks7.push_back(thread(&PolyEval::findAllMins, this, start_, end_));
    }
    for (unsigned int f=0; f<tasks7.size(); f++) {
        tasks7[f].join();
    }
}

// safely multithreads to create a matrix of pixel vals for MULTIPLE POLYNOMIALS AS ARRAY
void PolyEval::threadSafe_Sample8() {
    cout << "Making matrix..." << endl;
    vector<thread> tasks8;
    int startReal_ = 0;
    int endReal_ = 0;
    for (int j=0; j<numSamples; j++) {
        startReal_ = endReal_;
        endReal_ += realSpan;
        tasks8.push_back(thread(&PolyEval::createMat3, this, startReal_, endReal_));
    }
    for (unsigned int f=0; f<tasks8.size(); f++) {
        tasks8[f].join();
    }
}

// safely multithreads to find local mins for MULTIPLE POLYNOMIALS AS ARRAY
void PolyEval::threadSafe_Sample9() {
    cout << "Finding mins..." << endl;
    findAllMins2(0,numOfPoly);
    // vector<thread> tasks9;
    // int start_ = 0;
    // int end_ = 0;
    // for (int j=0; j<1; j++) {
    //     start_ = end_;
    //     end_ += numOfPoly;
    //     tasks9.push_back(thread(&PolyEval::findAllMins2, this, start_, end_));
    // }
    // for (unsigned int f=0; f<tasks9.size(); f++) {
    //     tasks9[f].join();
    // }
}

// Used to create a single matrix of min vals
void PolyEval::combineColRow() {
    for (int p=0; p<imgSize; p++) {
        for (int f=0; f<imgSize; f++) {
            binCountVect2[p][f] += colMinVect[p][f] + rowMinVect[p][f];
        }
    }
}

// Used to create a single matrix of max vals
void PolyEval::combineColRow2() {
    for (int p=0; p<imgSize; p++) {
        for (int f=0; f<imgSize; f++) {
            maxBinCountVect[p][f] += colMaxVect[p][f] + rowMaxVect[p][f];
        }
    }
}

// Used to create a single matrix of max vals from MULTIPLE POLYNOMIALS
void PolyEval::combineColRow3() {
    for (int p=0; p<imgSize; p++) {
        for (int f=0; f<imgSize; f++) {
            for (int z=0; z<numOfPoly; z++) {
                binCountVect3[p][f] += multBinValVect[p][f][z];
            }
        }
    }
}

// Used to create a single matrix of max vals from MULTIPLE POLYNOMIALS AS ARRAY
void PolyEval::combineColRow4() {
    cout << "Combining array..." << endl;

    // pads top and bottom with arr of zeros
    binCountArr3[0] = new int[imgSize];
    binCountArr3[imgSize -1] = new int[imgSize];
    for (int y=0; y<imgSize; y++) {
        binCountArr3[0][y] = 0;
        binCountArr3[imgSize -1][y] = 0;
    }

    // NOTE: Indices are funky because edges aren't calculated in findAllMins2(...)
    for (int p=1; p<(imgSize-1); p++) {
        binCountArr3[p] = new int[imgSize];
        binCountArr3[p][0] = 0; // pads left with zero
        for (int f=1; f<(imgSize-1); f++) {
            binCountArr3[p][f] = multBinValArr[p][f][0];
            for (int z=1; z<numOfPoly; z++) {
                binCountArr3[p][f] += multBinValArr[p][f][z];
            }
        }
        binCountArr3[p][imgSize -1] = 0; // pads right with zero
        // binCountArr3[p] = tempArr4;
        // delete tempArr4;
    }
    // delete multBinValArr;
    // pads top and bottom with arr of zeros
    // int *tempArr5 = new int[imgSize];
    // for (int y=0; y<imgSize; y++) {
    //     tempArr5[y] = 0;
    // }
    // binCountArr3[0] = tempArr5;
    // binCountArr3[imgSize -1] = tempArr5;
    // delete tempArr5;
    // cout << "Combining array..." << endl;
    // // int *tempArr4 = new int[imgSize];
    // // NOTE: Indices are funky because edges aren't calculated in findAllMins2(...)
    // for (int p=1; p<(imgSize-1); p++) {
    //     int *tempArr4 = new int[imgSize];
    //     tempArr4[0] = 0; // pads left with zero
    //     for (int f=1; f<(imgSize-1); f++) {
    //         tempArr4[f] = multBinValArr[p][f][0];
    //         for (int z=1; z<numOfPoly; z++) {
    //             tempArr4[f] += multBinValArr[p][f][z];
    //         }
    //     }
    //     tempArr4[imgSize -1] = 0; // pads right with zero
    //     binCountArr3[p] = tempArr4;
    //     // delete tempArr4;
    // }

    // // pads top and bottom with arr of zeros
    // int *tempArr5 = new int[imgSize];
    // for (int y=0; y<imgSize; y++) {
    //     tempArr5[y] = 0;
    // }
    // binCountArr3[0] = tempArr5;
    // binCountArr3[imgSize -1] = tempArr5;
    // // delete tempArr5;
}

// returns bin count from original method (tolerances, no local mins)
vector<vector<int>> PolyEval::getBinCount() {
    threadSafe_Sample();
    return binCountVect;
}

// returns bin count from combined col/row local mins
vector<vector<int>> PolyEval::getBinCount2() {
    threadSafe_Sample2();
    threadSafe_Sample3();
    combineColRow();
    return binCountVect2;
}

// returns bin count from combined col/row local mins w/ MULTIPLE POLYNOMIALS
vector<vector<unsigned short int>> PolyEval::getBinCount3() {
    threadSafe_Sample6();
    threadSafe_Sample7();
    combineColRow3();
    return binCountVect3;
}

// returns bin count from combined col/row local mins w/ MULTIPLE POLYNOMIALS AS 2D ARRAY
int** PolyEval::getBinCount4() {
    threadSafe_Sample8();
    threadSafe_Sample9();
    combineColRow4();
    return binCountArr3;
}

// returns bin count from combined col/row local max
vector<vector<int>> PolyEval::getMaxBinCount() {
    threadSafe_Sample2();
    threadSafe_Sample4();
    combineColRow2();
    return maxBinCountVect;
}

// returns a vect of vects which each contain the val of the polynomial at a given pixel
vector<vector<double>> PolyEval::getPixelVal() {
    threadSafe_Sample2();
    return pixelValVect;
}

vector<vector<int>> PolyEval::getMinPeaks() {
    threadSafe_Sample2();
    threadSafe_Sample3();
    combineColRow();
    threadSafe_Sample5();
    return minPeaksVect;
}