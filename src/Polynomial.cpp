#include "../headers/Polynomial.hpp"

Polynomial::Polynomial(vector<double> realPart, vector<double> imgPart, int size, bool isArr) {
    realComp = realPart;
    imgComp = imgPart;
    arrSize = size;
}

void Polynomial::vectPoly() {
    for (int i=0; i<arrSize; i++) {
        vPoly.push_back({realComp[i],imgComp[i]});
    }
}

vector<complex<double>> Polynomial::getVectPoly() {
    vectPoly();
    return vPoly;
}