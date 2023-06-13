#include "../headers/AxesEval.hpp"


AxesEval::AxesEval(vector<complex<double>> vectPolyToUse_, Tools kit_, vector<double> realSpaced_, vector<double> imgSpaced_,
                            bool isVertical_, int polySize_, int imgSize_, int numSamples_, int numOfPoly_, double bigBoy_) {
    vectPolyToUse = vectPolyToUse_;
    kit = kit_;
    realSpaced = realSpaced_;
    imgSpaced = imgSpaced_;
    isVertical = isVertical_;
    polySize = polySize_;
    imgSize = imgSize_;
    numSamples = numSamples_;
    numOfPoly = numOfPoly_;
    bigBoy = bigBoy_;
    realSpan = floor(imgSize_ /numSamples_);
    polySpan = max(int(floor(numOfPoly /numSamples)), 1);
    // creates a 2d vector with (real x img) indices

    // for (int h=0; h<imgSize; h++) {
    //     totBinCountVect.push_back(vector<unsigned short int> (imgSize,0));
    // }

    if (isVertical) {
        for (int i=0; i<imgSize; i++) {
            pixValVect.push_back(vector<vector<float>>(3, vector<float>(numOfPoly, 0.0)));
            multiBinCountVect.push_back(vector<vector<int>>(3, vector<int>(numOfPoly, 0)));
            totBinCountVect.push_back(vector<int>(3,0));
        }
    }
    else {
        for (int i=0; i<3; i++) {
            pixValVect.push_back(vector<vector<float>>(imgSize, vector<float>(numOfPoly, 0.0)));
            multiBinCountVect.push_back(vector<vector<int>>(imgSize, vector<int>(numOfPoly, 0)));
            totBinCountVect.push_back(vector<int>(imgSize,0));
        }
    }
}

// creates a matrix of the vals of the polynomial at each pixel on img axis
void AxesEval::imgEvalPixel(int startImg, int endImg) {
    double realVal = 0.0;
    double imgVal = 0.0;

    for (int k=startImg; k<endImg; k++) {
        for (int m=0; m<3; m++) {
            realVal = realSpaced[k];
            imgVal = imgSpaced[m];
            // potentially use .swap here? look at QuartPolyEval.createMat2()
            pixValVect[m][k] = kit.horner9(vectPolyToUse, numOfPoly, polySize, {realVal,imgVal});
        }
    }
}

// creates a matrix of the vals of the polynomial at each pixel on real axis
void AxesEval::realEvalPixel(int startReal, int endReal) {
    double realVal = 0.0;
    double imgVal = 0.0;
    vector<float> tempVect;

    for (int k=startReal; k<endReal; k++) {
        for (int m=0; m<3; m++) {
            realVal = realSpaced[m];
            imgVal = imgSpaced[k];
            
            // potentially use .swap here? look at QuartPolyEval.createMat2()
            tempVect = kit.horner9(vectPolyToUse, numOfPoly, polySize, {realVal,imgVal});
            pixValVect[k][m].swap(tempVect);
        }
    }
}

// finds all local mins for a given poly along IMG axis
void AxesEval::findImgMins(int startPoly, int endPoly) {
    float curr;
    float top;
    float left;
    float right;
    float bott;
    float tota;
    for (int k=1; k<2; k++) {
        for (int m=1; m<(imgSize-1); m++) {
            for (int n=startPoly; n<endPoly; n++) {
                curr = pixValVect[k][m][n];
                top = pixValVect[k-1][m][n];
                left = pixValVect[k][m-1][n];
                right = pixValVect[k][m+1][n];
                bott = pixValVect[k+1][m][n];
                tota = curr+top+left+right+bott;
                if ((tota>0) && (curr<top) && (curr<left) && (curr<right) && (curr<bott)) {
                        multiBinCountVect[k][m][n] ++;
                }
            }
        }
    }
}

// finds all local mins for a given poly along REAL axis
void AxesEval::findRealMins(int startPoly, int endPoly) {
    float curr;
    float top;
    float left;
    float right;
    float bott;
    float tota;
    for (int m=1; m<(imgSize-1); m++) {
        for (int k=1; k<2; k++) {
            for (int n=startPoly; n<endPoly; n++) {
                curr = pixValVect[m][k][n];
                top = pixValVect[m][k-1][n];
                left = pixValVect[m-1][k][n];
                right = pixValVect[m+1][k][n];
                bott = pixValVect[m][k+1][n];
                tota = curr+top+left+right+bott;
                if ((tota>0) && (curr<top) && (curr<left) && (curr<right) && (curr<bott)) {
                    multiBinCountVect[m][k][n] ++;
                }
            }
        }
    }
}

// safely multithreads to find vals along img axis
void AxesEval::threadSafe_Sample() {
    vector<thread> tasks;
    int start_ = 0;
    int end_ = 0;
    for (int j=0; j<numSamples; j++) {
        start_ = end_;
        end_ += realSpan;
        tasks.push_back(thread(&AxesEval::imgEvalPixel, this, start_, end_));
    }
    for (unsigned int f=0; f<tasks.size(); f++) {
        tasks[f].join();
    }
}

// safely multithreads to find vals along real axis
void AxesEval::threadSafe_Sample2() {
    vector<thread> tasks;
    int start_ = 0;
    int end_ = 0;
    for (int j=0; j<numSamples; j++) {
        start_ = end_;
        end_ += realSpan;
        tasks.push_back(thread(&AxesEval::realEvalPixel, this, start_, end_));
    }
    for (unsigned int f=0; f<tasks.size(); f++) {
        tasks[f].join();
    }
}

// safely multithreads to find local mins
void AxesEval::threadSafe_Sample3() {
    vector<thread> tasks3;
    int start_ = 0;
    int end_ = 0;
    for (int j=0; j<numSamples; j++) {
        start_ = end_;
        end_ += polySpan;
        tasks3.push_back(thread(&AxesEval::findImgMins, this, start_, end_));
    }
    for (unsigned int f=0; f<tasks3.size(); f++) {
        tasks3[f].join();
    }
}

// safely multithreads to find local max
void AxesEval::threadSafe_Sample4() {
    vector<thread> tasks4;
    int start_ = 0;
    int end_ = 0;
    for (int j=0; j<numSamples; j++) {
        start_ = end_;
        end_ += polySpan;
        tasks4.push_back(thread(&AxesEval::findRealMins, this, start_, end_));
    }
    for (unsigned int f=0; f<tasks4.size(); f++) {
        tasks4[f].join();
    }
}

// Used to create a single matrix of max vals from MULTIPLE POLYNOMIALS AS ARRAY
void AxesEval::vertCombineAxes() {
    for (int p=0; p<imgSize; p++) {
        for (int f=0; f<3; f++) {
            for (int z=0; z<numOfPoly; z++) {
                totBinCountVect[p][f] += multiBinCountVect[p][f][z];
            }
        }
    }
}

void AxesEval::horiCombineAxes() {
    for (int p=0; p<3; p++) {
        for (int f=0; f<imgSize; f++) {
            for (int z=0; z<numOfPoly; z++) {
                totBinCountVect[p][f] += multiBinCountVect[p][f][z];
            }
        }
    }
}

// returns bin count from original method (tolerances, no local mins)
vector<vector<int>> AxesEval::vertGetBinCount() {
    threadSafe_Sample2();
    threadSafe_Sample4();
    vertCombineAxes();
    return totBinCountVect;
}

vector<vector<int>> AxesEval::horiGetBinCount() {
    threadSafe_Sample();
    threadSafe_Sample3();
    horiCombineAxes();
    return totBinCountVect;
}