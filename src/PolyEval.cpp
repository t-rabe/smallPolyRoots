#include "../headers/PolyEval.hpp"


PolyEval::PolyEval(vector<complex<double>> vectPolyToUse_, Tools kit_, vector<double> realSpaced_, vector<double> imgSpaced_, int polySize_, int imgSize_, int numSamples_, double bigBoy_) {
    vectPolyToUse = vectPolyToUse_;
    kit = kit_;
    realSpaced = realSpaced_;
    imgSpaced = imgSpaced_;
    polySize = polySize_;
    imgSize = imgSize_;
    numSamples = numSamples_;
    bigBoy = bigBoy_;
    realSpan = floor(realSpaced_.size() /numSamples_);
    
    // creates a 2d vector with (real x img) indices
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
    // threadSafe_Sample();
    // threadSafe_Sample2();
    // threadSafe_Sample3();
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

// safely multithreads each sample to find roots within a certain tolerance
void PolyEval::threadSafe_Sample() {
    vector<thread> tasks;
    int startReal_ = 0;
    int endReal_ = 0;
    for (int j=0; j<numSamples; j++) {
        startReal_ = endReal_;
        endReal_ += realSpan;
        // unsigned int seeder = (unsigned int) time(NULL);
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
        // unsigned int seeder = (unsigned int) time(NULL);
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
        // unsigned int seeder = (unsigned int) time(NULL);
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
        // unsigned int seeder = (unsigned int) time(NULL);
        tasks4.push_back(thread(&PolyEval::findColMax, this, start_, end_));
        tasks4.push_back(thread(&PolyEval::findRowMax, this, start_, end_));
    }
    for (unsigned int f=0; f<tasks4.size(); f++) {
        tasks4[f].join();
    }
}

// Used to create a single matrix of min vals
void PolyEval::combineColRow() {
    for (int p=0; p<imgSize; p++) {
        for (int f=0; f<imgSize; f++) {
            binCountVect2[p][f] += colMinVect[p][f] + rowMinVect[p][f];
        }
    }
}

// Used to create a single matrix of min vals
void PolyEval::combineColRow2() {
    for (int p=0; p<imgSize; p++) {
        for (int f=0; f<imgSize; f++) {
            maxBinCountVect[p][f] += colMaxVect[p][f] + rowMaxVect[p][f];
        }
    }
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