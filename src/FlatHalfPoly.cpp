#include "../headers/FlatHalfPoly.hpp"


FlatHalfPoly::FlatHalfPoly(vector<complex<double>> vectPolyToUse_, Tools kit_, vector<double> realSpaced_, vector<double> imgSpaced_, int polySize_, int imgSize_, int numSamples_, double bigBoy_) {
    vectPolyToUse = vectPolyToUse_;
    kit = kit_;
    realSpaced = realSpaced_;
    imgSpaced = imgSpaced_;
    polySize = polySize_;
    imgSize = imgSize_ /2;
    twoImgSize = imgSize_;
    vectSize = twoImgSize *imgSize;
    numSamples = numSamples_;
    bigBoy = bigBoy_;
    realSpan = floor(realSpaced_.size() /numSamples_);
    
    // creates a 2d vector with (real x img) indices
    for (int f=0; f<vectSize; f++) {
        binCountVect.push_back(0);
    }

    for (int g=0; g<vectSize; g++) {
        binCountVect2.push_back(0);
    }

    for (int h=0; h<vectSize; h++) {
        maxBinCountVect.push_back(0);
    }

    for (int i=0; i<vectSize; i++) {
        pixelValVect.push_back(0.0);
    }

    for (int j=0; j<vectSize; j++) {
        colMinVect.push_back(0);
    }

    for (int k=0; k<vectSize; k++) {
        rowMinVect.push_back(0);
    }

    for (int m=0; m<vectSize; m++) {
        colMaxVect.push_back(0);
    }

    for (int n=0; n<vectSize; n++) {
        rowMaxVect.push_back(0);
    }

    for (int p=0; p<vectSize; p++) {
        minPeaksVect.push_back(0);
    }
}

// creates a matrix of the vals of the polynomial at each pixel
void FlatHalfPoly::createMat(int startReal, int endReal) {
    double realVal = 0.0;
    double imgVal = 0.0;
    int index = 0;
    for (int k=startReal; k<endReal; k++) {
        for (int m=0; m<imgSize; m++) {
            realVal = realSpaced[k];
            imgVal = imgSpaced[m];
            index = (k*imgSize) +m;
            pixelValVect[index] = abs(kit.horner7(vectPolyToUse, polySize, {realVal,imgVal}, bigBoy));
        }
    }
}

// evals individual pixels and updates their bin count. cannot factor out the found roots
// counts a pixel as a root if the val of the poly is "close to zero" (within a certain tolerance)
void FlatHalfPoly::evalPixel(int startReal, int endReal) {
    // int binNum = 0;
    int rootNum = 0;
    int index = 0;
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
                index = (k*imgSize) +m;
                binCountVect[index] ++;
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
void FlatHalfPoly::findColMin(int startReal, int endReal) {
    double prev;
    double curr;
    double next;
    int index = 0;
    for (int k=startReal; k<endReal; k++) {
        prev = pixelValVect[(k*imgSize)];
        curr = pixelValVect[(k*imgSize)];
        for (int m=1; m<imgSize; m++) {
            index = (k*imgSize) +m;
            next = pixelValVect[index];
            if ((curr <= prev) && (curr <= next)) {
                colMinVect[index-1] ++;
            }
            prev = curr;
            curr = next;
        }
        if (curr <= prev) {
            index = (k*imgSize) +imgSize -1;
            colMinVect[index] ++;
        }
    }
}


// CHANGE STARTIMG AND ENDIMG SO THAT THEY MATCH IMGSIZE / 2



// finds the local mins in every row to be compared with data from findColMin(...)
void FlatHalfPoly::findRowMin(int startImg, int endImg) {
    double prev2;
    double curr2;
    double next2;
    int index = 0;
    int gap = imgSize;
    for (int k=startImg; k<endImg; k++) {
        prev2 = pixelValVect[k];
        curr2 = pixelValVect[k];
        for (int m=1; m<twoImgSize; m++) {
            index = (m*gap) +k;
            next2 = pixelValVect[index];
            if ((curr2 <= next2) && (curr2 <= prev2)) {
                rowMinVect[((m*gap) -gap+k)] ++;
                // cout << prev2 << " " << curr2 << " " << next2 << '\n';
            }
            prev2 = curr2;
            curr2 = next2;
        }
        if (curr2 <= prev2) {
            rowMinVect[((twoImgSize*gap) -imgSize+k)] ++;
        }
    }
}

// finds the local max in every column to be compared with data from findRowMax(...)
void FlatHalfPoly::findColMax(int startReal, int endReal) {
    double prev;
    double curr;
    double next;
    int index = 0;
    for (int k=startReal; k<endReal; k++) {
        prev = pixelValVect[(k*imgSize)];
        curr = pixelValVect[(k*imgSize)];
        for (int m=1; m<imgSize; m++) {
            index = (k*imgSize) +m;
            next = pixelValVect[index];
            if ((curr >= prev) && (curr >= next)) {
                colMaxVect[index-1] --;
            }
            prev = curr;
            curr = next;
        }
        if (curr >= prev) {
            index = (k*imgSize) +imgSize-1;
            colMaxVect[index] --;
        }
    }
}

// finds the local max in every row to be compared with data from findColMax(...)
void FlatHalfPoly::findRowMax(int startImg, int endImg) {
    double prev2;
    double curr2;
    double next2;
    int index = 0;
    int gap = endImg - startImg;
    for (int k=startImg; k<endImg; k++) {
        prev2 = pixelValVect[k];
        curr2 = pixelValVect[k];
        for (int m=1; m<twoImgSize; m++) {
            index = (m*gap) +k;
            next2 = pixelValVect[index];
            if ((curr2 >= prev2) && (curr2 >= next2)) {
                rowMaxVect[((m*gap) -gap+k)] --;
            }
            prev2 = curr2;
            curr2 = next2;
        }
        if (curr2 >= prev2) {
            rowMaxVect[((twoImgSize*gap) -gap+k)] --;
        }
    }
}

// finds the local max in every row to be compared with data from findColMax(...)
void FlatHalfPoly::findMinPeaks(int startReal, int endReal) {
    int curr;
    int totNear;
    int index = 0;
    for (int k=(startReal+1); k<(endReal-1); k++) {
        for (int m=1; m<(imgSize-1); m++) {
            index = (k*imgSize) +m;
            curr = binCountVect2[index];
            if (curr == 2) {
                minPeaksVect[index] ++;
            }
            else if (curr == 1) {
                // the sum of all neighboring pixel values
                totNear = binCountVect2[index-imgSize-1]+binCountVect2[index-imgSize]+binCountVect2[index-imgSize+1]
                        + binCountVect2[index-1]+binCountVect2[index+1]
                        + binCountVect2[index+imgSize-1]+binCountVect2[index+imgSize]+binCountVect2[index+imgSize+1];
                if (totNear == 0) {
                    minPeaksVect[index] ++;
                }
            }
        }
    }
}

// safely multithreads each sample to find roots within a certain tolerance
void FlatHalfPoly::threadSafe_Sample() {
    vector<thread> tasks;
    int startReal_ = 0;
    int endReal_ = 0;
    for (int j=0; j<numSamples; j++) {
        startReal_ = endReal_;
        endReal_ += realSpan;
        tasks.push_back(thread(&FlatHalfPoly::evalPixel, this, startReal_, endReal_));
    }
    for (unsigned int f=0; f<tasks.size(); f++) {
        tasks[f].join();
    }
}

// safely multithreads to create a matrix of pixel vals
void FlatHalfPoly::threadSafe_Sample2() {
    vector<thread> tasks2;
    int startReal_ = 0;
    int endReal_ = 0;
    for (int j=0; j<numSamples; j++) {
        startReal_ = endReal_;
        endReal_ += realSpan;
        tasks2.push_back(thread(&FlatHalfPoly::createMat, this, startReal_, endReal_));
    }
    for (unsigned int f=0; f<tasks2.size(); f++) {
        tasks2[f].join();
    }
}

// safely multithreads to find local mins
void FlatHalfPoly::threadSafe_Sample3() {
    vector<thread> tasks3;
    int start_ = 0;
    int end_ = 0;
    for (int j=0; j<numSamples; j++) {
        start_ = end_;
        end_ += realSpan;
        tasks3.push_back(thread(&FlatHalfPoly::findColMin, this, start_, end_));
        tasks3.push_back(thread(&FlatHalfPoly::findRowMin, this, start_ /2, end_ /2));
    }
    for (unsigned int f=0; f<tasks3.size(); f++) {
        tasks3[f].join();
    }
}

// safely multithreads to find local max
void FlatHalfPoly::threadSafe_Sample4() {
    vector<thread> tasks4;
    int start_ = 0;
    int end_ = 0;
    for (int j=0; j<numSamples; j++) {
        start_ = end_;
        end_ += realSpan;
        tasks4.push_back(thread(&FlatHalfPoly::findColMax, this, start_, end_));
        tasks4.push_back(thread(&FlatHalfPoly::findRowMax, this, start_, end_));
    }
    for (unsigned int f=0; f<tasks4.size(); f++) {
        tasks4[f].join();
    }
}

// safely multithreads to find peaks of the min val vector
void FlatHalfPoly::threadSafe_Sample5() {
    vector<thread> tasks5;
    int start_ = 0;
    int end_ = realSpan;
    for (int j=0; j<numSamples; j++) {
        tasks5.push_back(thread(&FlatHalfPoly::findMinPeaks, this, start_, end_));
        start_ = end_ -2;
        end_ += realSpan;
    }
    for (unsigned int f=0; f<tasks5.size(); f++) {
        tasks5[f].join();
    }
}

// Used to create a single matrix of min vals
void FlatHalfPoly::combineColRow() {
    for (int p=0; p<vectSize; p++) {
        binCountVect2[p] += colMinVect[p] + rowMinVect[p];
    }
}

// Used to create a single matrix of max vals
void FlatHalfPoly::combineColRow2() {
    for (int p=0; p<vectSize; p++) {
        maxBinCountVect[p] += colMaxVect[p] + rowMaxVect[p];
    }
}

// returns bin count from original method (tolerances, no local mins)
vector<int> FlatHalfPoly::getBinCount() {
    threadSafe_Sample();
    return binCountVect;
}

// returns bin count from combined col/row local mins
vector<int> FlatHalfPoly::getBinCount2() {
    threadSafe_Sample2();
    threadSafe_Sample3();
    combineColRow();
    return binCountVect2;
}

// returns bin count from combined col/row local max
vector<int> FlatHalfPoly::getMaxBinCount() {
    threadSafe_Sample2();
    threadSafe_Sample4();
    combineColRow2();
    return maxBinCountVect;
}

// returns a vect of vects which each contain the val of the polynomial at a given pixel
vector<double> FlatHalfPoly::getPixelVal() {
    threadSafe_Sample2();
    return pixelValVect;
}

vector<int> FlatHalfPoly::getMinPeaks() {
    threadSafe_Sample2();
    threadSafe_Sample3();
    combineColRow();
    threadSafe_Sample5();
    return minPeaksVect;
}