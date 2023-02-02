#include <iomanip>
#include <iostream>
#include <complex>
#include <string>
#include <vector>
#include <thread>
#include <fstream>

using namespace std;

int main() {
    /**
     * THIS READS PREVIOUSLY SAVED COEFFS FROM A FILE
    */
    vector<double> vals; // values of each coefficient
    ifstream coeffFile;
    coeffFile.open("../src/testCoeffsOnesBigRandom2.csv"); // file holding old coeffs
    string coeff; // value to be read as string
    double coeffDoub; // value to be stored
    int lineNum = 0;
    int startIndex = 38265;
    int sideLen = 12;
    int matSize = pow(sideLen,2);
    if (coeffFile.is_open()) {
        while (coeffFile) {
            getline(coeffFile,coeff);
            coeffDoub = stod(coeff);
            if ((lineNum >= startIndex) && (lineNum < (startIndex + matSize))) {
                vals.push_back(coeffDoub);
            }
            lineNum ++;
        }
    }

    double currVal = 0;
    for (int i=0; i<sideLen; i++) {
        cout << "| ";
        for (int j=0; j<sideLen; j++) {
            currVal = vals[(i*sideLen)+j];
            if (currVal > 0) {
                cout << " ";
            }
            cout << currVal << " ";;
        }
        cout << "|\n";
    }
    return 0;   
}