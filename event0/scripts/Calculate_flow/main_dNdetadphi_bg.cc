#include<iostream>
#include <math.h>
#include<fstream>
#include<string>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include <vector>
#include <algorithm>

using namespace std;
// *** The code to calculate observables output from jat+parton cascade framework *** //
// Copyright @  Wenbin Zhao, 2023 //

// Function to calculate the difference in phi considering periodicity
double deltaPhi(double phi1, double phi2) {
    double dphi = abs(phi1 - phi2);
    if (dphi >= M_PI) dphi = 2 * M_PI - dphi;
    return dphi;
}


int main(int argc, char* argv[] )
{
    int pid, total_number_of_particles;

    // Read the binary output
       const char* filename = argv[1];

    // Open the binary file
    std::ifstream InStream;
    InStream.precision(15);
    InStream.open(filename, std::ios::in | std::ios::binary);

    // Check if the file was successfully opened
    if (!InStream.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return 2; // Return an error code
    }

    /*
    stringstream strparticle_name;
    strparticle_name << "final_state_hard_hadrons.bin";
    std::ifstream InStream;
    InStream.precision(15);
    string particle_name = strparticle_name.str();
    InStream.open(particle_name.c_str(), std::ios::in | std::ios::binary);
    */
    const int nBinsEta = 70; // Number of bins in the η direction
    const int nBinsPhi = 70; // Number of bins in the φ direction
    const double minEta = -3.50; // Minimum η value
    const double maxEta = 3.50; // Maximum η value
    const double minPhi = 0.; // Minimum φ value (-π)
    const double maxPhi = M_PI; // Maximum φ value (π)
    const int nBinsMultiplicity = 12;
    const int Nmix = 20;
    
    double mass, pT, phi, eta;
    float temp;
    double meanpT[11][4] = {0.}; int meanpT_event[11][4] = {0}; int Ntrg[12][3] = {0};
    //int MAMB_0_3[10] = {0}; int MAMB_0p3_3[10] = {0}; int MAMB_0p5_3[10] = {0}; 
    std::vector<int> NchVector; 
    double Nchbins[12] = {-1., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 500.,}; 
    
    
    int test[nBinsMultiplicity][nBinsEta][nBinsPhi] = {0};
    int Shistogram03[nBinsMultiplicity][nBinsEta][nBinsPhi] = {0};
    int Shistogram05[nBinsMultiplicity][nBinsEta][nBinsPhi] = {0};
    int Shistogram00[nBinsMultiplicity][nBinsEta][nBinsPhi] = {0};
    
    vector<double> etaValues00; vector<double> phiValues00;
    vector<double> etaValues03; vector<double> phiValues03;
    vector<double> etaValues05; vector<double> phiValues05;
    
    const char* filenameNtrg = "Ntrg.dat";

    // Open the file for reading
    std::ifstream inputFileNtrg(filenameNtrg);

    if (!inputFileNtrg.is_open()) {
        std::cerr << "Error: Unable to open file " << filenameNtrg << std::endl;
        return 1;
    }

    // Read the data from the file
    int Nevent;
    inputFileNtrg >> Nevent;
    if (Nmix > Nevent) {
        cout << "Need more events to calculate the background!" << endl;
        exit(1);
    }
    // Close the file
    inputFileNtrg.close();

    int etaBin, phiBin;
    for (int iev=0; iev<Nevent; iev++) {
        int Nchtemp = 0;
        etaValues00.clear(); etaValues03.clear(); etaValues05.clear();
        phiValues00.clear(); phiValues03.clear(); phiValues05.clear();
            
        if (iev < (Nevent - Nmix)) {
            for (int imix=0; imix < Nmix; imix++) {
                // Start one event
                InStream.read(reinterpret_cast<char*>(&total_number_of_particles), sizeof(int));
                for (int i=0; i<total_number_of_particles; i++) {
                     InStream.read(reinterpret_cast<char*>(&pid), sizeof(int));
                     InStream.read(reinterpret_cast<char*>(&temp), sizeof(float));
                     mass = temp;
                     InStream.read(reinterpret_cast<char*>(&temp), sizeof(float));
                     pT = temp;
                     InStream.read(reinterpret_cast<char*>(&temp), sizeof(float));
                     phi = temp;
                     InStream.read(reinterpret_cast<char*>(&temp), sizeof(float));
                     eta = temp;
                     // Calculate the dN/dDeltaphi_dDeltaeta
                     if (abs(pid)==211 || abs(pid) == 321 || abs(pid) == 2212) {
                         Nchtemp ++;
                         if (pT < 3.0) {
                             // pT > 0.5
                             if (pT > 0.5) {
                                 etaValues05.push_back(eta); phiValues05.push_back(phi);
                             }
                             // pT > 0.3
                             if (pT > 0.3) {
                                 etaValues03.push_back(eta); phiValues03.push_back(phi);
                             }
                             // pT > 0.0
                                 etaValues00.push_back(eta); phiValues00.push_back(phi);
                         }
                     } 
                }
            }
        }
        
        // Compute the differences in eta and phi and fill the histogram
        int binDeltaEta, binDeltaPhi;
        if (phiValues00.size() > 2) {
            for (int inch=0; inch<1; inch++) {
                for (int i = 0; i < phiValues00.size(); ++i) {
                    for (int j = i + 1; j < phiValues00.size(); ++j) {
                        double deltaEta = etaValues00[j] - etaValues00[i];
                        double deltaPhiVal = abs(phiValues00[j] - phiValues00[i]);
                        if (deltaPhiVal >= M_PI) deltaPhiVal = 2 * M_PI - deltaPhiVal;

                        // Map deltaEta and deltaPhi to histogram bins
                        binDeltaEta = static_cast<int>((deltaEta - minEta) / (maxEta - minEta) * nBinsEta);
                        binDeltaPhi = static_cast<int>((deltaPhiVal - minPhi) / (maxPhi - minPhi) * nBinsPhi);
                        // Check bounds
                        if (binDeltaEta < nBinsEta && binDeltaEta >-1) {
                            Shistogram00[inch][binDeltaEta][binDeltaPhi]++;
                        }
                        
                    }
                }
                
                for (size_t i = 0; i < etaValues03.size(); ++i) {
                    for (size_t j = i + 1; j < etaValues03.size(); ++j) {
                        double deltaEta = etaValues03[j] - etaValues03[i];
                        double deltaPhiVal = abs(phiValues03[j] - phiValues03[i]);
                        if (deltaPhiVal >= M_PI) deltaPhiVal = 2 * M_PI - deltaPhiVal;

                        // Map deltaEta and deltaPhi to histogram bins
                        int binDeltaEta = static_cast<int>((deltaEta - minEta) / (maxEta - minEta) * nBinsEta);
                        int binDeltaPhi = static_cast<int>((deltaPhiVal - minPhi) / (maxPhi - minPhi) * nBinsPhi);
        
                        // Check bounds
                        if (binDeltaEta < nBinsEta && binDeltaEta >-1 ) {
                            Shistogram03[inch][binDeltaEta][binDeltaPhi]++;
                        }
                    }
                }
                
                for (size_t i = 0; i < etaValues05.size(); ++i) {
                    for (size_t j = i + 1; j < etaValues05.size(); ++j) {
                        double deltaEta = etaValues05[j] - etaValues05[i];
                        double deltaPhiVal = abs(phiValues05[j] - phiValues05[i]);
                        if (deltaPhiVal >= M_PI) deltaPhiVal = 2 * M_PI - deltaPhiVal;
    
                        // Map deltaEta and deltaPhi to histogram bins
                        int binDeltaEta = static_cast<int>((deltaEta - minEta) / (maxEta - minEta) * nBinsEta);
                        int binDeltaPhi = static_cast<int>((deltaPhiVal - minPhi) / (maxPhi - minPhi) * nBinsPhi);
    
                        // Check bounds
                        if (binDeltaEta < nBinsEta && binDeltaEta >-1) {
                            Shistogram05[inch][binDeltaEta][binDeltaPhi]++;
                        }
                    }
                }
            }
        }
    }
    
    /*
    // Print the 2D histogram
    for (int i = 0; i < nBinsEta; ++i) {
        for (int j = 0; j < nBinsPhi; ++j) {
            cout << Shistogram03[i][j] << "\t";
        }
        cout << endl;
    }
    */
    
    // Output file
    char output_filename2[256]; 
    const char* output_filename_format;
    for (int inch=0; inch<1; inch++) {
        output_filename_format = "d2N_dDeltaeta_dDeltaphi_bg_00_3_%d.dat";
        std::sprintf(output_filename2, output_filename_format, inch);
        std::ofstream output00(output_filename2);
        if (!output00.is_open() ) {
            cout << "cannot open output file:"<< endl
                 << output_filename2 << endl;
            return -1;
        }
        for (auto i = 0; i < nBinsEta; ++i) {
            for (int j = 0; j < nBinsPhi; ++j) {
                output00  << Shistogram00[inch][i][j] << "\t";
            }
            output00 << endl;
        }
        output00.close();
        
        output_filename_format = "d2N_dDeltaeta_dDeltaphi_bg_03_3_%d.dat";
        std::sprintf(output_filename2, output_filename_format, inch);
        std::ofstream output2(output_filename2);
        if (!output2.is_open() ) {
            cout << "cannot open output file:"<< endl
                 << output_filename2 << endl;
            return -1;
        }
        for (auto i = 0; i < nBinsEta; ++i) {
            //output2 << i * (maxEta - minEta) / nBinsEta * 1.0 + minEta << "\t";
            for (int j = 0; j < nBinsPhi; ++j) {
                output2  << Shistogram03[inch][i][j] << "\t";
            }
            output2 << endl;
        }
        output2.close();
        
        output_filename_format = "d2N_dDeltaeta_dDeltaphi_bg_05_3_%d.dat";
        std::sprintf(output_filename2, output_filename_format, inch);
        std::ofstream output3(output_filename2);
        if (!output3.is_open() ) {
            cout << "cannot open output file:"<< endl
                 << output_filename2 << endl;
            return -1;
        }
        for (auto i = 0; i < nBinsEta; ++i) {
            for (int j = 0; j < nBinsPhi; ++j) {
                output3  << Shistogram05[inch][i][j] << "\t";
            }
            output3 << endl;
        }
        output3.close();
        
    }
    InStream.close();
}

