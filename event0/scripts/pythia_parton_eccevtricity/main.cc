#include<iostream>
#include<cmath>
#include<fstream>
#include<string>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh" 

using namespace fastjet;

using namespace std;


struct Particle {
    int pdgid;
    double px, py, pz, e, x, y, z, t;
};

int main(int argc, char* argv[] )
{

    // parameter setting
    // selection for final particles which are used to reconstruct jet
    double absetamax = 2.4;
    double particle_ptmin = 0.3; // CMS cut, CMS PAS HIN-21-013
    // parameter setting
    const double R_jet = 0.8; // CMS cut, CMS PAS HIN-21-013
    const double jet_ptmin = 500.0;
    double jet_absetamax = 1.6;
    
    double Eng[261][261];

    int partons;
    int temp;
    double bcoll;
    double DX = 0.1; double DY = 0.1;

    int itemp,ip,pId, mid1, mid2;
    char sch;
    double fx,fy,fz,ft,fetas;
    
    //fastjet setting
    // create a jet definition: 
    // a jet algorithm with a given radius parameter
   // select jet
    fastjet::Selector jet_selector = fastjet::SelectorAbsEtaMax( jet_absetamax ) && fastjet::SelectorPtMin( jet_ptmin );
    
    // selection for final particles which are used to reconstruct jet
    fastjet::Selector particle_selector = fastjet::SelectorAbsEtaMax(absetamax) && fastjet::SelectorPtMin( particle_ptmin );
    
    // Open the input file
    ifstream inputFile("parton_info.dat");
    if (!inputFile.is_open()) {
        cerr << "Error opening the file!" << endl;
        return 1;
    }
    
    char outfilename[128];
    sprintf(outfilename, "eccentricity.dat");
    FILE* outfile;
    outfile = fopen(outfilename,"w");

    vector<Particle> particles;
    string line;
    string line2;

    // Skip the first two lines (headers)
    getline(inputFile, line); // Read and ignore the first line
    for (int iev = 0; iev <20000; iev ++) {
        particles.clear();
        int numParticles;
        vector<fastjet::PseudoJet> input_particles;
        if (inputFile.eof()) {
           cout << "end the file" << endl;
           break;
        }
        
        std::getline(inputFile, line2);
        std::istringstream iss(line2);

        if (line2.empty() || line2[0] == '#') {
            // Skip comments and empty lines
            continue;
        }
        iss >> numParticles;
        // Read each particle data
        for (int i = 0; i < numParticles; ++i) {
            Particle particle;
            inputFile >> particle.pdgid >> particle.px >> particle.py >> particle.pz
                      >> particle.e >> particle.x >> particle.y >> particle.z >> particle.t >> mid1 >> mid2;
            particles.push_back(particle);
            fastjet::PseudoJet particlefastjet = PseudoJet(particle.px, particle.py, particle.pz, particle.e);
            particlefastjet.set_user_index(particle.pdgid);
            input_particles.push_back(particlefastjet);
        }
        input_particles = particle_selector(input_particles);
        
        
        
        // Then do the jet finding
        fastjet::Selector particle_selector = fastjet::SelectorAbsEtaMax(absetamax) && fastjet::SelectorPtMin( particle_ptmin );
        fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R_jet);
        // select jet
        fastjet::Selector jet_selector = fastjet::SelectorAbsEtaMax( jet_absetamax ) && fastjet::SelectorPtMin( jet_ptmin );
        input_particles = particle_selector(input_particles);
        fastjet::ClusterSequence clust_seq(input_particles, jet_def);
        // get the resulting jets ordered in pt
        vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt( jet_selector(clust_seq.inclusive_jets()) );
        
        if (inclusive_jets.size() == 0 || inclusive_jets[0].pt() < 550.){
            continue;
        }        
        // Select the jet pT and output the selected events, rotate the jet at pz direction
        for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
            // Select the jet with jet pT > 550 GeV/c.
            if (inclusive_jets[i].pt() > 250.) {
                //njetevent_count++;
                double jx = inclusive_jets[i].px(); double jy = inclusive_jets[i].py(); double jz = inclusive_jets[i].pz();
                double theta = acos(jz/sqrt(jx*jx + jy*jy + jz*jz));
                double costheta = cos(theta);
                double sintheta = sin(theta);
                double rx = jy/sqrt(jx*jx + jy*jy); double ry =  -jx/sqrt(jx*jx + jy*jy); double rz =  0.;

                vector<fastjet::PseudoJet> constituents = inclusive_jets[i].constituents();
                //output << njetevent_count << "  " << constituents.size() << endl; 
                int Nparticle = 0;
                for (unsigned int jj=0; jj<constituents.size(); jj++) {
                    if (constituents[jj].pt() > 0.3 && abs(constituents[jj].eta()) < 2.4 ) Nparticle ++;
                }
                
                double* particleEng;
                double* particleSx;
                double* particleSy;
                particleEng = new double[constituents.size()+1];
                particleSx = new double[constituents.size()+1];
                particleSy = new double[constituents.size()+1];
                int count = 0;
    
                //int size_temp = constituents.size();
                fprintf(outfile,"%20ld  ",constituents.size());
                for (unsigned int jj=0; jj<constituents.size(); jj++) {
                    if (constituents[jj].pt() <= 0.3 || abs(constituents[jj].eta()) >= 2.4 ) continue;
                    // rotate into jet going direction into the z-axis (0, 0, 1)
                    /*
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
                    */
                    //int pidtemp = constituents[jj].user_index();
                    // find the coordinte information
                    int which_particle = 0;
                    for (int ipart=0; ipart<numParticles; ipart++ ) {
                        if (abs((constituents[jj].px() - particles[ipart].px ) /constituents[jj].px()) < 0.001 &&
                            abs((constituents[jj].pz() - particles[ipart].pz ) /constituents[jj].pz()) < 0.001 ) {
                            which_particle = ipart;
                            break;
                        }
                    }
                    fx = particles[which_particle].x;
                    fy = particles[which_particle].y;
                    double energy = particles[which_particle].e;
                    particleEng[count] = energy;//(energy*ft - particles[which_particle].e*fz)/tau0;
                    particleSx[count] = fx;
                    particleSy[count] = fy;
                    count++;
                }
                // Then calculate the en
                double totalex = 0.; double totaley = 0.; double totale = 0.;
                for(int i = 0; i!= count; i++) {
                    double sx = particleSx[i];
                    double sy = particleSy[i];
                    double energy = particleEng[i];
                    totalex += energy * sx;
                    totaley += energy * sy;
                    totale  += energy;
                }
                double xmean = totalex/totale;
                double ymean = totaley/totale;
                double ecc_real[10];
                double ecc_imag[10];
                for(int i=0;i!=10;i++) {
                    ecc_real[i] = 0.;
                    ecc_imag[i] = 0.;
                }
                for (int ipart=0; ipart<count; ipart++) {
                    double posx = particleSx[ipart] - xmean;
                    double posy = particleSy[ipart] - ymean;
                    double phi =atan2(posy,posx);
                    double rT = sqrt(posx*posx + posy*posy);
                    for(int n=2;n!=4;n++) {
                        ecc_real[n] += particleEng[ipart]*pow(rT,n)*cos(n*phi);
                        ecc_imag[n] += particleEng[ipart]*pow(rT,n)*sin(n*phi);
                    }
                    ecc_real[0] += particleEng[ipart]*pow(rT,2);
                    ecc_imag[0] += particleEng[ipart]*pow(rT,3);
                }
                fprintf(outfile,"%8d   %20.6lf  %20.6lf  ", 0, ecc_real[0], ecc_imag[0]);
                for(int n=2;n!=4;n++) {
                    fprintf(outfile,"%8d   %20.6lf  %20.6lf  ", n, ecc_real[n], ecc_imag[n]);
                }
                fprintf(outfile,"\n");
                 delete particleEng;
                 delete particleSx;
                 delete particleSy;
            }
        }
    }
    inputFile.close();
    fclose(outfile);
	return 0;
}



