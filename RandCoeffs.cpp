#include "./RandCoeffs.hpp"

// initializes a vect<vect<float>> of randomly ordered indices
RandCoeffs::RandCoeffs(int numPoints_, int numSamples_, int genNum_) {
    numPoints = numPoints_;
    numSamples = numSamples_;
    genNum = genNum_;
    for (int g=0; g<numSamples; g++) {
        vector<double> interMed(floor(numPoints /numSamples),0.0);
        samples.push_back(interMed);
    }
    threadSafe_Sample();
}

// fills individual samples with randomly ordered indices
void RandCoeffs::getRandomValue(unsigned int seeder, int index) {
    std::mt19937 gen(seeder);
    // std::mt19937 gen(10);
    // cout << "gen = " << gen << endl;
    std::uniform_real_distribution<double> dist(-1,1);
    
    // fill an entire column (one sample)
    for (int k=0; k<(numPoints /numSamples); k++) {
        samples[index][k] = (dist(gen));
        // samples[index][k] = (dist(seeder));
    }
}

// safely multithreads each sample and assigns unique seed to quickly create vects
void RandCoeffs::threadSafe_Sample() {
    vector<thread> tasks;
    int jStart = genNum; // a different num for each instantiation (real vs img)
    for (int j=jStart; j<(numSamples+jStart); j++) {
        unsigned int seeder = (unsigned int) time(NULL);
        tasks.push_back(thread(&RandCoeffs::getRandomValue, this, seeder, (j-jStart)));
    }
    for (unsigned int f=0; f<tasks.size(); f++) {
        tasks[f].join();
    }
}

// returns individual samples
vector<vector<double>> RandCoeffs::getSample() {
    return samples;
}

// returns the whole random vector in 1 dimension
vector<double> RandCoeffs::getTotSample() {
    vector<double> totSample;
    for (int i=0; i<samples.size(); i++) {
        for (int k=0; k<samples[0].size(); k++) {
            totSample.push_back(samples[i][k]);
        }
    }
    // cout << totSample[12] << endl;
    return totSample;
}