#include "../headers/Polynomial.hpp"

Polynomial::Polynomial(vector<double> realPart, vector<double> imgPart, int size, complex<double>* blankPoly, bool isArr) {
    realComp = realPart;
    imgComp = imgPart;
    aPoly = blankPoly;
    arrSize = size;
    if (isArr) {
        arrPoly();
    }
    else {
        vectPoly();
    }
}

void Polynomial::arrPoly() {
    for (int i=0; i<arrSize; i++) {
        aPoly[i] = {realComp[i],imgComp[i]};
    }
}

void Polynomial::vectPoly() {
    for (int i=0; i<arrSize; i++) {
        vPoly.push_back({realComp[i],imgComp[i]});
    }
}

complex<double>* Polynomial::getArrPoly() {
    return aPoly;
}

vector<complex<double>> Polynomial::getVectPoly() {
    return vPoly;
}