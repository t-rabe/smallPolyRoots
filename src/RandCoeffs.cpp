#include "../headers/RandCoeffs.hpp"

// initializes a vect<vect<float>> of randomly ordered indices
RandCoeffs::RandCoeffs(int numPoints_, int numSamples_, int genNum_) {
    numPoints = numPoints_;
    numSamples = numSamples_;
    genNum = genNum_;
    for (int g=0; g<numSamples; g++) {
        vector<double> interMed(floor(numPoints /numSamples),0.0);
        samples.push_back(interMed);
    }

    for (int g=0; g<numSamples; g++) {
        vector<double> interMed1(floor(numPoints /numSamples),0.0);
        samplesOnes.push_back(interMed1);
    }
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

// fills individual samples with randomly ordered indices OF 1 OR -1
void RandCoeffs::getRandomValueOnes(unsigned int seeder, int index) {
    std::mt19937 gen(seeder);
    // std::mt19937 gen(10);
    // cout << "gen = " << gen << endl;
    std::uniform_int_distribution<int> dist(0,10);
    
    // fill an entire column (one sample)
    for (int k=0; k<(numPoints /numSamples); k++) {
        samplesOnes[index][k] = (((dist(gen))%2) *2.0) -1.0;
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

// same as above but ONLY 1 OR -1
void RandCoeffs::threadSafe_SampleOnes() {
    vector<thread> tasks;
    int jStart = genNum; // a different num for each instantiation (real vs img)
    for (int j=jStart; j<(numSamples+jStart); j++) {
        unsigned int seeder = (unsigned int) time(NULL);
        tasks.push_back(thread(&RandCoeffs::getRandomValueOnes, this, seeder, (j-jStart)));
    }
    for (unsigned int f=0; f<tasks.size(); f++) {
        tasks[f].join();
    }
}

// returns individual samples
vector<vector<double>> RandCoeffs::getSample() {
    threadSafe_Sample();
    return samples;
}

// same as above but ONLY 1 OR -1
vector<vector<double>> RandCoeffs::getSampleOnes() {
    threadSafe_SampleOnes();
    return samplesOnes;
}

// returns the whole random vector in 1 dimension
vector<double> RandCoeffs::getTotSample() {
    threadSafe_Sample();
    vector<double> totSample;
    for (long unsigned int i=0; i<samples.size(); i++) {
        for (long unsigned int k=0; k<samples[0].size(); k++) {
            totSample.push_back(samples[i][k]);
        }
    }
    // cout << totSample[12] << endl;
    return totSample;
}

// same as above but ONLY 1 OR -1
vector<double> RandCoeffs::getTotSampleOnes() {
    threadSafe_SampleOnes();
    vector<double> totSample;
    for (long unsigned int i=0; i<samplesOnes.size(); i++) {
        for (long unsigned int k=0; k<samplesOnes[0].size(); k++) {
            totSample.push_back(samplesOnes[i][k]);
        }
    }
    return totSample;
}