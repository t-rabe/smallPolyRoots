#include "../headers/PolyEval.hpp"


PolyEval::PolyEval(vector<complex<double>> vectPolyToUse_, Tools kit_, vector<double> realSpaced_, vector<double> imgSpaced_, int polySize_, int imgSize_, int numSamples_) {
    vectPolyToUse = vectPolyToUse_;
    kit = kit_;
    realSpaced = realSpaced_;
    imgSpaced = imgSpaced_;
    polySize = polySize_;
    imgSize = imgSize_;
    numSamples = numSamples_;
    realSpan = floor(realSpaced_.size() /numSamples_);
    
    // creates a 2d vector with (real x img) indices
    for (int g=0; g<imgSize; g++) {
        vector<int> interMed(imgSize,0);
        binCountVect.push_back(interMed);
    }
    threadSafe_Sample();
}

// evals individual pixels and updates their bin count. cannot factor out the found roots
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
            if (abs(kit.horner7(vectPolyToUse, polySize, {realVal,imgVal})) < tolerance) {
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

// safely multithreads each sample and assigns unique seed to quickly create vects
void PolyEval::threadSafe_Sample() {
    vector<thread> tasks;
    int startReal_ = 0;
    int endReal_ = 0;
    for (int j=0; j<numSamples; j++) {
        startReal_ = endReal_;
        endReal_ += realSpan;
        unsigned int seeder = (unsigned int) time(NULL);
        tasks.push_back(thread(&PolyEval::evalPixel, this, startReal_, endReal_));
    }
    for (unsigned int f=0; f<tasks.size(); f++) {
        tasks[f].join();
    }
}

vector<vector<int>> PolyEval::getBinCount() {
    return binCountVect;
}