#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh" 

/*
    The script to calculate the jet pT and save the data. 
*/

using namespace std;
using namespace fastjet;

int main(int argv, char* argc[])
{
    int Nevent = atoi(argc[1]);
    // selection for final particles which are used to reconstruct jet
    double absetamax = 2.4;
    double particle_ptmin = 0.3; // CMS cut, CMS PAS HIN-21-013
    // parameter setting
    const double R_jet = 0.8; // CMS cut, CMS PAS HIN-21-013
    const double jet_ptmin = 500.0;
    double jet_absetamax = 1.6;
    vector<fastjet::PseudoJet> input_particles;
    char inputfile[128];
    sprintf(inputfile, "particle_list.dat");
    FILE* infile;
    infile = fopen(inputfile,"r");
    char stemp1[100];
    char** stemp2;
    int total_number_of_particles, pid, event_id, int_temp, status;
    double px, py, pz, energy, mass, dummpx, dummpy, dummpz, dummpt, weight;
    int event_loop_flag = 1;
    int count_event_number = 0;
    /*
    //output the results to file
    char output_filename[128];
    sprintf(output_filename,"jet_hardons");
    ofstream output(output_filename);
    if( ! output.is_open() ) {
        cout << "cannot open output file:"<< endl
             << output_filename << endl;
        return -1;
    }
    */
    // open file for output
    std::string binary_output_filename = "final_state_hard_hadrons.bin";
    remove(binary_output_filename.c_str());
    FILE *outbin = NULL;
    outbin = fopen(binary_output_filename.c_str(), "wb");
    
    int njetevent_count = 0;
    for (int iev = 0; iev < Nevent; iev ++) {
        if(feof(infile)) {
            event_loop_flag = 0;
            cout << " End the event loop ~~~ " << endl;
            break;
        }
        fscanf(infile,"%s %d\n",stemp1, &total_number_of_particles);
        input_particles.clear();
        for (auto i=0; i<total_number_of_particles; i++) {
            if(feof(infile)) {
                event_loop_flag = 0;
                cout << " End the event loop, and drop last event ~~~ " << endl;
                break;
            }
            fscanf(infile,"%d %lf %lf %lf %lf %lf\n", &pid, &mass,
                   &energy, &px, &py, &pz);
            if (isnan(energy) || isnan(px) || isnan(py) || isnan(pz)) continue;
            fastjet::PseudoJet particle = PseudoJet(px, py, pz, energy);
            particle.set_user_index(pid);
            input_particles.push_back(particle);
        }
        if(event_loop_flag == 0) {
            cout << " End the event loop and drop last event ~~~ " << endl;
            break;
        }
        count_event_number++;

        // Then do the jet finding
        fastjet::Selector particle_selector = fastjet::SelectorAbsEtaMax(absetamax) && fastjet::SelectorPtMin( particle_ptmin );
        fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R_jet);
        // select jet
        fastjet::Selector jet_selector = fastjet::SelectorAbsEtaMax( jet_absetamax ) && fastjet::SelectorPtMin( jet_ptmin );
        input_particles = particle_selector(input_particles);
        fastjet::ClusterSequence clust_seq(input_particles, jet_def);
        // get the resulting jets ordered in pt
        vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt( jet_selector(clust_seq.inclusive_jets()) );

        // Select the jet pT and output the selected events, rotate the jet at pz direction
        for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
            // Select the jet with jet pT > 550 GeV/c.
            if (inclusive_jets[i].pt() > 550.) {
                njetevent_count++;
                double jx = inclusive_jets[i].px(); double jy = inclusive_jets[i].py(); double jz = inclusive_jets[i].pz();
                double theta = acos(jz/sqrt(jx*jx + jy*jy + jz*jz));
                double costheta = cos(theta);
                double sintheta = sin(theta);
                double rx = jy/sqrt(jx*jx + jy*jy); double ry =  -jx/sqrt(jx*jx + jy*jy); double rz =  0.;

                vector<fastjet::PseudoJet> constituents = inclusive_jets[i].constituents();
                //output << njetevent_count << "  " << constituents.size() << endl; 
                int Nparticle = 0;
                for (unsigned int jj=0; jj<constituents.size(); jj++){
                    if (constituents[jj].pt() > 0.3 && abs(constituents[jj].eta()) < 2.4 ) Nparticle ++;
                }
                //int size_temp = constituents.size();
                fwrite(&Nparticle, sizeof(int), 1, outbin);
                for (unsigned int jj=0; jj<constituents.size(); jj++){
                    if (constituents[jj].pt() <= 0.3 || abs(constituents[jj].eta()) >= 2.4 ) continue;
                    // rotate into jet going direction into the z-axis (0, 0, 1)
                    double ppx = constituents[jj].px() * (costheta + (1.-costheta)*rx*rx)  
                               + constituents[jj].py() * ((1.-costheta)*rx*ry - sintheta*rz) 
                               + constituents[jj].pz() * ((1.-costheta)*rx*rz+sintheta*ry); 
                    double ppy = constituents[jj].px() * ((1.-costheta)*rx*ry + sintheta*rz) 
                               + constituents[jj].py() * (costheta + (1.-costheta)*ry*ry)   
                               + constituents[jj].pz() * ((1.-costheta)*ry*rz - sintheta*rx);
                    double ppz = constituents[jj].px() * ((1.-costheta)*rx*rz - sintheta*ry) 
                               + constituents[jj].py() * ((1.-costheta)*ry*rz + sintheta*rx) 
                               + constituents[jj].pz() * (costheta + (1.-costheta)*rz*rz);
                    double pmag = sqrt(ppx*ppx + ppy*ppy + ppz*ppz);

                    // Output the constituents inside the jet with the jet going to +z direction
                    int pidtemp = constituents[jj].user_index();
                    fwrite(&pidtemp, sizeof(int), 1, outbin);
                    float array[] = {
                        static_cast<float>(constituents[jj].m()),
                        static_cast<float>(sqrt(ppx*ppx + ppy*ppy)), static_cast<float>(atan2(ppy, ppx)),
                        static_cast<float>(0.5*log((pmag+ppz)/(pmag-ppz))),
                    };
                    fwrite(array, sizeof(float), 4, outbin);
                    /*
                    output << constituents[jj].user_index() << "  " << constituents[jj].m() << "  " 
                           << sqrt(ppx*ppx + ppy*ppy) << "  " << atan2(ppy, ppx) << "  " 
                           << 0.5*log((pmag+ppz)/(pmag-ppz))
                           << endl;
                   */
                }
            }
        }
    }

    fclose(infile);
    //output.close();
    fclose(outbin);
    return 0;
}

