#include "../headers/AxesEval.hpp"


AxesEval::AxesEval(vector<complex<double>> vectPolyToUse_, Tools kit_, vector<double> realSpaced_, vector<double> imgSpaced_,
                                                int polySize_, int imgSize_, int numSamples_, int numOfPoly_, double bigBoy_) {
    vectPolyToUse = vectPolyToUse_;
    kit = kit_;
    realSpaced = realSpaced_;
    imgSpaced = imgSpaced_;
    polySize = polySize_;
    imgSize = imgSize_ -1;
    halfImgSize = (imgSize_ /2) -1;
    numSamples = numSamples_;
    numOfPoly = numOfPoly_;
    bigBoy = bigBoy_;
    realSpan = floor(imgSize /numSamples_);
    polySpan = max(int(floor(numOfPoly /numSamples)), 1);

    // creates a 2d vector with (real x img) indices

    for (int h=0; h<imgSize; h++) {
        vector<int> interMed0(imgSize,0);
        totBinCountVect.push_back(interMed0);
    }

    for (int i=0; i<3; i++) {
        vector<vector<float>> interMed1(imgSize, vector<float>(numOfPoly, 0.0));
        imgPixValVect.push_back(interMed1);
    }
    
    for (int j=0; j<imgSize; j++) {
        vector<vector<float>> interMed2(3, vector<float>(numOfPoly, 0.0));
        realPixValVect.push_back(interMed2);
    }

    for (int k=0; k<imgSize; k++) {
        vector<int> interMed3(numOfPoly,0);
        imgMultiBinCountVect.push_back(interMed3);
    }

    for (int m=0; m<imgSize; m++) {
        vector<int> interMed4(numOfPoly, 0);
        realMultiBinCountVect.push_back(interMed4);
    }
}

// creates a matrix of the vals of the polynomial at each pixel on img axis
void AxesEval::imgEvalPixel(int startImg, int endImg) {
    double realVal = 0.0;
    double imgVal = 0.0;

    for (int k=startImg; k<endImg; k++) {
        for (int m=0; m<3; m++) {
            realVal = realSpaced[m+halfImgSize+1];
            imgVal = imgSpaced[k];
            // potentially use .swap here? look at QuartPolyEval.createMat2()
            imgPixValVect[m][k] = kit.horner9(vectPolyToUse, numOfPoly, polySize, {realVal,imgVal}, bigBoy);
        }
    }
}

// creates a matrix of the vals of the polynomial at each pixel on real axis
void AxesEval::realEvalPixel(int startReal, int endReal) {
    double realVal = 0.0;
    double imgVal = 0.0;

    for (int k=startReal; k<endReal; k++) {
        for (int m=0; m<3; m++) {
            realVal = realSpaced[k];
            imgVal = imgSpaced[m+halfImgSize+1];
            
            // potentially use .swap here? look at QuartPolyEval.createMat2()
            realPixValVect[k][m] = kit.horner9(vectPolyToUse, numOfPoly, polySize, {realVal,imgVal}, bigBoy);
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
    for (int m=1; m<(imgSize-1); m++) {
        for (int n=startPoly; n<endPoly; n++) {
            curr = imgPixValVect[1][m][n];
            top = imgPixValVect[0][m][n];
            left = imgPixValVect[1][m-1][n];
            right = imgPixValVect[1][m+1][n];
            bott = imgPixValVect[2][m][n];
            tota = curr+top+left+right+bott;
            if ((tota>0) && (curr<=top) && (curr<=left) && (curr<=right) && (curr<=bott)) {
                    imgMultiBinCountVect[m][n] ++;
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
        for (int n=startPoly; n<endPoly; n++) {
            curr = realPixValVect[m][1][n];
            top = realPixValVect[m][0][n];
            left = realPixValVect[m-1][1][n];
            right = realPixValVect[m+1][1][n];
            bott = realPixValVect[m][2][n];
            tota = curr+top+left+right+bott;
            if ((tota>0) && (curr<=top) && (curr<=left) && (curr<=right) && (curr<=bott)) {
                    realMultiBinCountVect[m][n] ++;
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
void AxesEval::combineAxes() {
    // cout << "Combining axes..." << endl;

    int binCount;
    vector<int> fillerZeros(halfImgSize, 0);
    
    // writes the top half of img axis
    for (int x=0; x<halfImgSize; x++) {
        binCount = 0;
        for (int i=0; i<numOfPoly; i++) {
            binCount += imgMultiBinCountVect[x][i];
        }
        totBinCountVect[x][halfImgSize] = binCount;
    }
    // this part writes the real axis
    for (int y=0; y<imgSize; y++) {
        binCount = 0;
        for (int j=0; j<numOfPoly; j++) {
            binCount += realMultiBinCountVect[y][j];
        }
        totBinCountVect[halfImgSize][y] = binCount;
    }
    // writes the bottom half of img axis
    for (int z=0; z<halfImgSize; z++) {
        binCount = 0;
        for (int k=0; k<numOfPoly; k++) {
            binCount += imgMultiBinCountVect[z+halfImgSize+1][k];
        }
        totBinCountVect[z+halfImgSize+1][halfImgSize] = binCount;
    }

}

// returns bin count from original method (tolerances, no local mins)
vector<vector<int>> AxesEval::getBinCount() {
    threadSafe_Sample();
    threadSafe_Sample2();
    threadSafe_Sample3();
    threadSafe_Sample4();
    combineAxes();
    return totBinCountVect;
}