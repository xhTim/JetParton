#include<iostream>
#include <math.h>
#include<fstream>
#include<string>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include <vector>
#include <algorithm>

// *** The code to calculate observables output from jat+parton cascade framework *** //
// Copyright @  Wenbin Zhao, 2023 //

using namespace std;
int main(int argc, char* argv[] )
{
    int pid, total_number_of_particles;
    // Output file
    string output_filename2;
    output_filename2 = "QAQB.dat";// output files of final hadrons
    string output_filenamench;
    output_filenamench = "Nch.dat";// output files of final hadrons
    
    string output_filenamencheta;
    output_filenamencheta = "dNchdeta.dat";// output files of final hadrons
    
    string output_filenamepnch;
    output_filenamepnch = "pNch.dat";// output files of final hadrons
    

    // Read the binary output
       const char* filename = argv[1];
       // std::string fullPath = "..//" + std::string(relativePath);
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
    const int roweta = 100; const double deta = 0.1; int Deltaeta = 2.0/deta;
    const double minEta = 0.0; const double maxEta = 10.;

    int Nchtemp;
    double mass, pT, eta;
    double phi = 9999.;
    float temp, tempphi;
    double pTlab, philab, etalab;
    double eta_dis[11][90]= {0.}; int eta_nevent_count[90] = {0};
    double QnAB_0p0_3[11][4] = {0.0}; double QnAB_0p3_3[11][4] = {0.0}; double QnAB_0p5_3[11][4] = {0.0};
    double bQnAB_0p0_3[11][4] = {0.0}; double bQnAB_0p3_3[11][4] = {0.0}; double bQnAB_0p5_3[11][4] = {0.0};
    double cQnAB_0p0_3[11][4] = {0.0}; double cQnAB_0p3_3[11][4] = {0.0}; double cQnAB_0p5_3[11][4] = {0.0};
    
    double AQnAB_0p0_3[11][4] = {0.0}; double AQnAB_0p3_3[11][4] = {0.0}; double AQnAB_0p5_3[11][4] = {0.0};
    double AbQnAB_0p0_3[11][4] = {0.0}; double AbQnAB_0p3_3[11][4] = {0.0}; double AbQnAB_0p5_3[11][4] = {0.0};
    double AcQnAB_0p0_3[11][4] = {0.0}; double AcQnAB_0p3_3[11][4] = {0.0}; double AcQnAB_0p5_3[11][4] = {0.0};
    
    double meanpT[11][4] = {0.}; int meanpT_event[11][4] = {0};
    int Nchsum[11] = {0}; int Ncheventcount[11] = {0};
    int pNch[50] = {0};
    //int MAMB_0_3[10] = {0}; int MAMB_0p3_3[10] = {0}; int MAMB_0p5_3[10] = {0}; 
    std::vector<int> NchVector; 
    //double Nchbins[12] = {-1., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 500.,}; 
        double Nchbins[12] = {-1., 10., 20., 30., 40., 47., 55., 65., 75., 85., 95., 500.,}; 
    //int Nchtemp;
    /*
    double QnA_0p0_3_real[roweta][4] = {0.}; double QnA_0p3_3_real[roweta][4] = {0.}; double QnA_0p5_3_real[roweta][4] = {0.};
        double QnB_0p0_3_real[roweta][4] = {0.}; double QnB_0p3_3_real[roweta][4] = {0.}; double QnB_0p5_3_real[roweta][4] = {0.};
        double QnA_0p0_3_imag[roweta][4] = {0.}; double QnA_0p3_3_imag[roweta][4] = {0.}; double QnA_0p5_3_imag[roweta][4] = {0.};
        double QnB_0p0_3_imag[roweta][4] = {0.}; double QnB_0p3_3_imag[roweta][4] = {0.}; double QnB_0p5_3_imag[roweta][4] = {0.};
        double meanpT_temp[4] = {0.};    int meanpT_event_temp[4] = {0};
    */

    while( InStream.read(reinterpret_cast<char*>(&total_number_of_particles), sizeof(int))) {
        Nchtemp = 0;
        double QnA_0p0_3_real[roweta][4] = {0.}; double QnA_0p3_3_real[roweta][4] = {0.}; double QnA_0p5_3_real[roweta][4] = {0.}; 
        double QnB_0p0_3_real[roweta][4] = {0.};double QnB_0p3_3_real[roweta][4] = {0.}; double QnB_0p5_3_real[roweta][4] = {0.}; 
        double QnA_0p0_3_imag[roweta][4] = {0.}; double QnA_0p3_3_imag[roweta][4] = {0.};double QnA_0p5_3_imag[roweta][4] = {0.}; 
        double QnB_0p0_3_imag[roweta][4] = {0.}; double QnB_0p3_3_imag[roweta][4] = {0.}; double QnB_0p5_3_imag[roweta][4] = {0.}; 
        double meanpT_temp[4] = {0.};  double meanpT_event_temp[4] = {0};
        
        double eta_dis_temp[90]= {0.};
        
        // Start one event
        for (auto i=0; i<total_number_of_particles; i++) {
            InStream.read(reinterpret_cast<char*>(&pid), sizeof(int));
            InStream.read(reinterpret_cast<char*>(&temp), sizeof(float));
            mass = temp;
            InStream.read(reinterpret_cast<char*>(&temp), sizeof(float));
            pT = temp;
            InStream.read(reinterpret_cast<char*>(&tempphi), sizeof(float));
            if (tempphi != phi) {
            phi = tempphi;
            InStream.read(reinterpret_cast<char*>(&temp), sizeof(float));
            eta = temp;
            /*
            InStream.read(reinterpret_cast<char*>(&temp), sizeof(float));
            pTlab = temp;
            InStream.read(reinterpret_cast<char*>(&temp), sizeof(float));
            philab = temp;
            InStream.read(reinterpret_cast<char*>(&temp), sizeof(float));
            etalab = temp;
            */
           // cout <<total_number_of_particles << " " <<  pid << "  " << mass << "  " << pT << "  " << phi << "  " << eta << endl;
            if (abs(pid) == 211) {
                meanpT_event_temp[1] = meanpT_event_temp[1] +1;
                meanpT_temp[1] = meanpT_temp[1] + pT;
            }
            if (abs(pid) == 321) {
                meanpT_event_temp[2] = meanpT_event_temp[2] +1;
                meanpT_temp[2] = meanpT_temp[2] + pT;
            }
            if (abs(pid) == 2212) {
                meanpT_event_temp[3] = meanpT_event_temp[3] +1;
                meanpT_temp[3] = meanpT_temp[3] + pT;
            }
            
            if (abs(pid)==211 || abs(pid) == 321 || abs(pid) == 2212) {
                meanpT_event_temp[0] = meanpT_event_temp[0] +1;
                meanpT_temp[0] = meanpT_temp[0] + pT;
                Nchtemp =Nchtemp +1;
                // Calculate the dNch/deta*
                int binEta = static_cast<int>((eta - 0) / (9.0 - 0.0) * 90);
                eta_dis_temp[binEta] = eta_dis_temp[binEta] + 1.0/0.1;
                if (pT < 3.0) {
                    // pT > 0.5
                    int binDeltaEta = static_cast<int>((eta - minEta) / (maxEta - minEta) * roweta);
                    if (binDeltaEta > roweta) continue;
                    if (binDeltaEta < 0) continue;
                    if (pT > 0.5) {
                        for (auto iorder = 0; iorder < 4; iorder++) {
                            //if (eta > 1.0) {
                                QnA_0p5_3_real[binDeltaEta][iorder] = QnA_0p5_3_real[binDeltaEta][iorder] + cos(iorder * 1. * phi);
                                QnA_0p5_3_imag[binDeltaEta][iorder] = QnA_0p5_3_imag[binDeltaEta][iorder] + sin(iorder * 1. * phi);
                            //}
                            //if (eta < -1.0) {
                                QnB_0p5_3_real[binDeltaEta][iorder] = QnB_0p5_3_real[binDeltaEta][iorder] + cos(iorder * 1. * phi);
                                QnB_0p5_3_imag[binDeltaEta][iorder] = QnB_0p5_3_imag[binDeltaEta][iorder] + sin(iorder * 1. * phi);
                            //}
                        }
                    }
                    // pT > 0.3
                    if (pT > 0.3) {
                        for (auto iorder = 0; iorder < 4; iorder++) {
                            //if (eta > 1.0) {
                                QnA_0p3_3_real[binDeltaEta][iorder] = QnA_0p3_3_real[binDeltaEta][iorder] + cos(iorder * 1. * phi);
                                QnA_0p3_3_imag[binDeltaEta][iorder] = QnA_0p3_3_imag[binDeltaEta][iorder] + sin(iorder * 1. * phi);
                            //}
                            // if (eta < -1.0) {
                                QnB_0p3_3_real[binDeltaEta][iorder] = QnB_0p3_3_real[binDeltaEta][iorder] + cos(iorder * 1. * phi);
                                QnB_0p3_3_imag[binDeltaEta][iorder] = QnB_0p3_3_imag[binDeltaEta][iorder] + sin(iorder * 1. * phi);
                            // }
                        }
                    }
                    // pT > 0.0
                        for (auto iorder = 0; iorder < 4; iorder++) {
                            //if (eta > 1.0) {
                                QnA_0p0_3_real[binDeltaEta][iorder] = QnA_0p0_3_real[binDeltaEta][iorder] + cos(iorder * 1. * phi);
                                QnA_0p0_3_imag[binDeltaEta][iorder] = QnA_0p0_3_imag[binDeltaEta][iorder] + sin(iorder * 1. * phi);
                            //}
                            //if (eta < -1.0) {
                                QnB_0p0_3_real[binDeltaEta][iorder] = QnB_0p0_3_real[binDeltaEta][iorder] + cos(iorder * 1. * phi);
                                QnB_0p0_3_imag[binDeltaEta][iorder] = QnB_0p0_3_imag[binDeltaEta][iorder] + sin(iorder * 1. * phi);
                            //}
                        }
                }
            }
            }
        }
        
        
        if (Nchtemp <= 1) continue;
        //cout << Nchtemp << endl;
        if (Nchtemp < 500.) NchVector.push_back(Nchtemp);
        // Then Calculate the QAQB
        int binNch = static_cast<int>(Nchtemp / 3);
        if (binNch<50) pNch[binNch]++;
        for (auto inch=0; inch<11; inch++) {
            if (Nchtemp > Nchbins[inch] && Nchtemp <= Nchbins[inch+1]) {
                Ncheventcount[inch]++;
                Nchsum[inch] = Nchsum[inch] + Nchtemp;
                
                for (int ieta=0; ieta<90; ieta++) {
                    eta_dis[inch][ieta] = eta_dis[inch][ieta] + eta_dis_temp[ieta];
                }
                for (int iorder=0; iorder<4; iorder++) {
                    // Delta\eta = 2.0
                    double sumtep00 = 0.; double sumtep03 = 0.0; double sumtep05 = 0.;
                    for (int ieta=0; ieta<(roweta - Deltaeta); ieta++) {
                        for (int ietab=(ieta+Deltaeta); ietab<roweta; ietab++) {
                            sumtep00 = sumtep00 + QnA_0p0_3_real[ieta][iorder] * QnB_0p0_3_real[ietab][iorder] 
                                                + QnA_0p0_3_imag[ieta][iorder] * QnB_0p0_3_imag[ietab][iorder];
                            sumtep03 = sumtep03 + QnA_0p3_3_real[ieta][iorder] * QnB_0p3_3_real[ietab][iorder] 
                                                + QnA_0p3_3_imag[ieta][iorder] * QnB_0p3_3_imag[ietab][iorder];
                            sumtep05 = sumtep05 + QnA_0p5_3_real[ieta][iorder] * QnB_0p5_3_real[ietab][iorder] 
                                                + QnA_0p5_3_imag[ieta][iorder] * QnB_0p5_3_imag[ietab][iorder];
                                                
                        }
                    }
                    
                    QnAB_0p0_3[inch][iorder] = QnAB_0p0_3[inch][iorder] + sumtep00;
                    QnAB_0p3_3[inch][iorder] = QnAB_0p3_3[inch][iorder] + sumtep03;
                    QnAB_0p5_3[inch][iorder] = QnAB_0p5_3[inch][iorder] + sumtep05;
                    
                    // Delta\eta = 1.0
                    sumtep00 = 0.; sumtep03 = 0.0; sumtep05 = 0.;
                    for (int ieta=0; ieta<(roweta - Deltaeta+5); ieta++) {
                        for (int ietab=(ieta+Deltaeta-5); ietab<roweta; ietab++) {
                            sumtep00 = sumtep00 + QnA_0p0_3_real[ieta][iorder] * QnB_0p0_3_real[ietab][iorder] 
                                                + QnA_0p0_3_imag[ieta][iorder] * QnB_0p0_3_imag[ietab][iorder];
                            sumtep03 = sumtep03 + QnA_0p3_3_real[ieta][iorder] * QnB_0p3_3_real[ietab][iorder] 
                                                + QnA_0p3_3_imag[ieta][iorder] * QnB_0p3_3_imag[ietab][iorder];
                            sumtep05 = sumtep05 + QnA_0p5_3_real[ieta][iorder] * QnB_0p5_3_real[ietab][iorder] 
                                                + QnA_0p5_3_imag[ieta][iorder] * QnB_0p5_3_imag[ietab][iorder];
                                                
                        }
                    }
                    
                    bQnAB_0p0_3[inch][iorder] = bQnAB_0p0_3[inch][iorder] + sumtep00;
                    bQnAB_0p3_3[inch][iorder] = bQnAB_0p3_3[inch][iorder] + sumtep03;
                    bQnAB_0p5_3[inch][iorder] = bQnAB_0p5_3[inch][iorder] + sumtep05;
                    
                    // Delta\eta = 3.0
                    sumtep00 = 0.; sumtep03 = 0.0; sumtep05 = 0.;
                    for (int ieta=0; ieta<(roweta - Deltaeta-5); ieta++) {
                        for (int ietab=(ieta+Deltaeta+5); ietab<roweta; ietab++) {
                            sumtep00 = sumtep00 + QnA_0p0_3_real[ieta][iorder] * QnB_0p0_3_real[ietab][iorder] 
                                                + QnA_0p0_3_imag[ieta][iorder] * QnB_0p0_3_imag[ietab][iorder];
                            sumtep03 = sumtep03 + QnA_0p3_3_real[ieta][iorder] * QnB_0p3_3_real[ietab][iorder] 
                                                + QnA_0p3_3_imag[ieta][iorder] * QnB_0p3_3_imag[ietab][iorder];
                            sumtep05 = sumtep05 + QnA_0p5_3_real[ieta][iorder] * QnB_0p5_3_real[ietab][iorder] 
                                                + QnA_0p5_3_imag[ieta][iorder] * QnB_0p5_3_imag[ietab][iorder];
                                                
                        }
                    }
                    
                    cQnAB_0p0_3[inch][iorder] = cQnAB_0p0_3[inch][iorder] + sumtep00;
                    cQnAB_0p3_3[inch][iorder] = cQnAB_0p3_3[inch][iorder] + sumtep03;
                    cQnAB_0p5_3[inch][iorder] = cQnAB_0p5_3[inch][iorder] + sumtep05;
                    
                    
                    meanpT[inch][iorder]           = meanpT[inch][iorder] + meanpT_temp[iorder];
                    meanpT_event[inch][iorder]     = meanpT_event[inch][iorder] + meanpT_event_temp[iorder];
                }
                // Add +1 for bin effects
                for (int iorder=0; iorder<4; iorder++) {
                    // Delta\eta = 2.0
                    double sumtep00 = 0.; double sumtep03 = 0.0; double sumtep05 = 0.;
                    for (int ieta=0; ieta<(roweta - Deltaeta-1); ieta++) {
                        for (int ietab=(ieta+Deltaeta+1); ietab<roweta; ietab++) {
                            sumtep00 = sumtep00 + QnA_0p0_3_real[ieta][iorder] * QnB_0p0_3_real[ietab][iorder] 
                                                + QnA_0p0_3_imag[ieta][iorder] * QnB_0p0_3_imag[ietab][iorder];
                            sumtep03 = sumtep03 + QnA_0p3_3_real[ieta][iorder] * QnB_0p3_3_real[ietab][iorder] 
                                                + QnA_0p3_3_imag[ieta][iorder] * QnB_0p3_3_imag[ietab][iorder];
                            sumtep05 = sumtep05 + QnA_0p5_3_real[ieta][iorder] * QnB_0p5_3_real[ietab][iorder] 
                                                + QnA_0p5_3_imag[ieta][iorder] * QnB_0p5_3_imag[ietab][iorder];
                                                
                        }
                    }
                    
                    AQnAB_0p0_3[inch][iorder] = AQnAB_0p0_3[inch][iorder] + sumtep00;
                    AQnAB_0p3_3[inch][iorder] = AQnAB_0p3_3[inch][iorder] + sumtep03;
                    AQnAB_0p5_3[inch][iorder] = AQnAB_0p5_3[inch][iorder] + sumtep05;
                    
                    // Delta\eta = 1.0
                    sumtep00 = 0.; sumtep03 = 0.0; sumtep05 = 0.;
                    for (int ieta=0; ieta<(roweta - Deltaeta-1+5); ieta++) {
                        for (int ietab=(ieta+Deltaeta+1-5); ietab<roweta; ietab++) {
                            sumtep00 = sumtep00 + QnA_0p0_3_real[ieta][iorder] * QnB_0p0_3_real[ietab][iorder] 
                                                + QnA_0p0_3_imag[ieta][iorder] * QnB_0p0_3_imag[ietab][iorder];
                            sumtep03 = sumtep03 + QnA_0p3_3_real[ieta][iorder] * QnB_0p3_3_real[ietab][iorder] 
                                                + QnA_0p3_3_imag[ieta][iorder] * QnB_0p3_3_imag[ietab][iorder];
                            sumtep05 = sumtep05 + QnA_0p5_3_real[ieta][iorder] * QnB_0p5_3_real[ietab][iorder] 
                                                + QnA_0p5_3_imag[ieta][iorder] * QnB_0p5_3_imag[ietab][iorder];
                                                
                        }
                    }
                    
                    AbQnAB_0p0_3[inch][iorder] = AbQnAB_0p0_3[inch][iorder] + sumtep00;
                    AbQnAB_0p3_3[inch][iorder] = AbQnAB_0p3_3[inch][iorder] + sumtep03;
                    AbQnAB_0p5_3[inch][iorder] = AbQnAB_0p5_3[inch][iorder] + sumtep05;
                    
                    // Delta\eta = 3.0
                    sumtep00 = 0.; sumtep03 = 0.0; sumtep05 = 0.;
                    for (int ieta=0; ieta<(roweta - Deltaeta-1-5); ieta++) {
                        for (int ietab=(ieta+Deltaeta+1+5); ietab<roweta; ietab++) {
                            sumtep00 = sumtep00 + QnA_0p0_3_real[ieta][iorder] * QnB_0p0_3_real[ietab][iorder] 
                                                + QnA_0p0_3_imag[ieta][iorder] * QnB_0p0_3_imag[ietab][iorder];
                            sumtep03 = sumtep03 + QnA_0p3_3_real[ieta][iorder] * QnB_0p3_3_real[ietab][iorder] 
                                                + QnA_0p3_3_imag[ieta][iorder] * QnB_0p3_3_imag[ietab][iorder];
                            sumtep05 = sumtep05 + QnA_0p5_3_real[ieta][iorder] * QnB_0p5_3_real[ietab][iorder] 
                                                + QnA_0p5_3_imag[ieta][iorder] * QnB_0p5_3_imag[ietab][iorder];
                                                
                        }
                    }
                    
                    AcQnAB_0p0_3[inch][iorder] = AcQnAB_0p0_3[inch][iorder] + sumtep00;
                    AcQnAB_0p3_3[inch][iorder] = AcQnAB_0p3_3[inch][iorder] + sumtep03;
                    AcQnAB_0p5_3[inch][iorder] = AcQnAB_0p5_3[inch][iorder] + sumtep05;
                }
            }
        }
    }
    // Then output the Qn and NchVector
    
    ofstream output2(output_filename2.c_str()); 
    if (!output2.is_open() ) {
        cout << "cannot open output file:"<< endl
         << output_filename2 << endl;
        return -1;
    }
    output2 << " # Qn ( n = 0,1,2,3) of 0.0<pT<3.0, 0.3<pT<3.0, 0.5<pT<3.0 total_pT event_of_total_pT (charged, pi, k, p); Then the Qn with Delta_eta = 1.0 and 3.0" << endl;
    for (auto inch=0; inch<11; inch++) {
        for (int iorder=0; iorder<4; iorder++) {
            output2 << QnAB_0p0_3[inch][iorder] << "  "
                    << QnAB_0p3_3[inch][iorder] << "  "
                    << QnAB_0p5_3[inch][iorder] << "  "
                    << meanpT[inch][iorder] << "  "
                    << meanpT_event[inch][iorder] << "  ";
        }
        
        for (int iorder=0; iorder<4; iorder++) {
            output2 << bQnAB_0p0_3[inch][iorder] << "  "
                    << bQnAB_0p3_3[inch][iorder] << "  "
                    << bQnAB_0p5_3[inch][iorder] << "  ";
        }
        
        for (int iorder=0; iorder<4; iorder++) {
            output2 << cQnAB_0p0_3[inch][iorder] << "  "
                    << cQnAB_0p3_3[inch][iorder] << "  "
                    << cQnAB_0p5_3[inch][iorder] << "  ";
        }
        
        for (int iorder=0; iorder<4; iorder++) {
            output2 << AQnAB_0p0_3[inch][iorder] << "  "
                    << AQnAB_0p3_3[inch][iorder] << "  "
                    << AQnAB_0p5_3[inch][iorder] << "  ";
        }
        
        for (int iorder=0; iorder<4; iorder++) {
            output2 << AbQnAB_0p0_3[inch][iorder] << "  "
                    << AbQnAB_0p3_3[inch][iorder] << "  "
                    << AbQnAB_0p5_3[inch][iorder] << "  ";
        }
        
        for (int iorder=0; iorder<4; iorder++) {
            output2 << AcQnAB_0p0_3[inch][iorder] << "  "
                    << AcQnAB_0p3_3[inch][iorder] << "  "
                    << AcQnAB_0p5_3[inch][iorder] << "  ";
        }
        
        
        output2 << Ncheventcount[inch] << "  " << Nchsum[inch];
        output2 << endl;
    }
        
    ofstream outputnch(output_filenamench.c_str()); 
    if (!outputnch.is_open() ) {
        cout << "cannot open output file:"<< endl
         << output_filenamench << endl;
        return -1;
    }
    for (int ii=0; ii<NchVector.size(); ii++) {
        outputnch << NchVector[ii] << endl;
    }
    
    ofstream outputncheta(output_filenamencheta.c_str()); 
    if (!outputncheta.is_open() ) {
        cout << "cannot open output file:"<< endl
         << output_filenamencheta << endl;
        return -1;
    }
    
    for (int ii=0; ii<11; ii++) {
        outputncheta << Ncheventcount[ii] << "  ";
    }
    outputncheta << endl;
        
    ofstream outputpnch(output_filenamepnch.c_str()); 
    if (!outputpnch.is_open() ) {
        cout << "cannot open output file:"<< endl
         << output_filenamepnch << endl;
        return -1;
    }
    for (int ii=0; ii<50; ii++) {
        outputpnch << ii*3 << " " <<  ii*3 +3 << "  " << pNch[ii] << "  " << endl;
    }
    outputpnch << endl;
    
    for (int ieta=0; ieta<90; ieta++) {
        for (int ii=0; ii<11; ii++) {
            outputncheta << eta_dis[ii][ieta] << "  ";
        }
        outputncheta << endl;
    }
    
    
    NchVector.clear(); 
    output2.close();
    outputnch.close();
    InStream.close();
    outputncheta.close();
    outputpnch.close();
}

