#include <iostream>    //....C++ headers
#include <fstream>
#include <string>
#include <ctime>

#include <unistd.h>   // access() funciton
#include <sys/stat.h>  // mkdir() function 

#include "Pythia8/Pythia.h"    //....PYTHIA8 headers
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh" 

using namespace Pythia8;
using namespace std;
using namespace fastjet;

#define pi 3.1415926

//////////////////////////////////////////
float ran33(long *idum);
void rotate(double px,double py,double pz,double pr[4],int icc);
void get_coordinate(vector<double>& parton_x, vector<double>& parton_y, const string& pathin);

// Function to recursively find all particles and their second-generation decay products
void findDecayProducts(const Event& event, int particleIndex, vector<int>& decayedParticles);
////////////////////////////////////////////

int main(int argv, char* argc[])
{
    string random_str = string(argc[2]);
    int totalevent = atoi(argc[1]);
    // Set up the Pythia8 configuration.
    Pythia pythia;
    
    pythia.readString("Random:setSeed = on");
    //pythia.readString("Random:seed = 0");  //92549 792 458 871
    pythia.readString("Random:seed = " + random_str);
    pythia.readString("Beams:eCM = 13000.0");//5020,2760 GeV
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212"); //2212 is proton; 1000791970 for 197Au; 1000822080 for 208Pb;
    // Standard settings
    pythia.readString("HardQCD:all = on");
    //pythia.readString("Tune:pp=18");
    pythia.readString("PhaseSpace:pTHatMin = 500.0");
    pythia.readString("PhaseSpace:pTHatMax = -1.");

    pythia.readString("Tune:pp=14");
    pythia.readString("PDF:pSet = 20");
    pythia.readString("Tune:ee=7");
    pythia.readString("MultipartonInteractions:ecmPow=0.03");
    pythia.readString("MultipartonInteractions:bProfile=2");
    pythia.readString("MultipartonInteractions:pT0Ref=1.41");
    pythia.readString("MultipartonInteractions:coreRadius=0.76");
    pythia.readString("MultipartonInteractions:coreFraction=0.63");
    pythia.readString("ColourReconnection:range=5.18");
    pythia.readString("SigmaTotal:zeroAXB=off");
    pythia.readString("SpaceShower:rapidityOrder=on");
    pythia.readString("SpaceShower:alphaSorder=2");
    pythia.readString("SpaceShower:alphaSvalue=0.118");
    pythia.readString("SigmaProcess:alphaSvalue=0.118");
    pythia.readString("SigmaProcess:alphaSorder=2");
    pythia.readString("MultipartonInteractions:alphaSvalue=0.118");
    pythia.readString("MultipartonInteractions:alphaSorder=2");
    pythia.readString("TimeShower:alphaSorder=2");
    pythia.readString("TimeShower:alphaSvalue=0.118");
    pythia.readString("SigmaTotal:mode = 0");
    pythia.readString("SigmaTotal:sigmaEl = 21.89");
    pythia.readString("SigmaTotal:sigmaTot = 100.309");
    //pythia.readString("PDF:pSet=LHAPDF6:NNPDF31_nnlo_as_0118");
    
    /*
    pythia.readString("PDF:useHard = on");
    pythia.readString("PDF:nPDFBeamB=100791970");
    pythia.readString("PDF:useHardNPDFB = on");
    pythia.readString("PDF:nPDFSetB=3");
    pythia.readString("PDF:nPDFBeamA=1000791970");
    pythia.readString("PDF:useHardNPDFA = on");
    pythia.readString("PDF:nPDFSetA=3");
    pythia.readString("TimeShower:pTmin=2.0"); // the mimum value of showed parton Wenbin
    pythia.readString("TimeShower:nGluonToQuark=3");
    */
    pythia.readString("HadronLevel:Hadronize = on");
    // Initialize Pythia
    if (!pythia.init()) {
        std::cerr << "Pythia initialization failed!" << std::endl;
        return 1;
    }

    // Output the parton information
    /*
    string output_filename;
    output_filename = "hadorn_info.dat";
    ofstream output_parton(output_filename.c_str());  //-----zhao

    if( ! output_parton.is_open() ) {
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
    
    // parameter setting
    // selection for final particles which are used to reconstruct jet
    double absetamax = 2.4;
    double particle_ptmin = 0.3; // CMS cut, CMS PAS HIN-21-013
    // parameter setting
    const double R_jet = 0.8; // CMS cut, CMS PAS HIN-21-013
    const double jet_ptmin = 500.0;
    double jet_absetamax = 1.6;
    
    //fastjet setting
    // create a jet definition: 
    // a jet algorithm with a given radius parameter
   // select jet
    fastjet::Selector jet_selector = fastjet::SelectorAbsEtaMax( jet_absetamax ) && fastjet::SelectorPtMin( jet_ptmin );
    
    // selection for final particles which are used to reconstruct jet
    fastjet::Selector particle_selector = fastjet::SelectorAbsEtaMax(absetamax) && fastjet::SelectorPtMin( particle_ptmin );
    
    // Number of events to generate
    int numEvent = 0;
    // Loop over events
    for (int iEvent = 0; ; ++iEvent) {
        if (!pythia.next()) continue;
        vector<fastjet::PseudoJet> input_particles;
        // Access the event record
        Event& event = pythia.event;
        int Nparton = 0;
        
        // First use the fastjet to pre-select the event with jet pT > 500 GeV at parton level
        for (int j = 0; j < event.size(); ++j) {
            Particle& particle = event[j];
            if (particle.isFinal()) {
                Nparton++;
                fastjet::PseudoJet particlefastjet = PseudoJet(particle.px(), particle.py(), particle.pz(), particle.e());
                particlefastjet.set_user_index(particle.id());
                input_particles.push_back(particlefastjet);
            }
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
        } else {
            numEvent +=1;
            if ( numEvent > totalevent) break;
        }
        
                // Select the jet pT and output the selected events, rotate the jet at pz direction
        for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
            // Select the jet with jet pT > 550 GeV/c.
            if (inclusive_jets[i].pt() > 550.) {
                //njetevent_count++;
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
                        //static_cast<float>(constituents[jj].pt()),
                        //static_cast<float>(atan2(constituents[jj].py(), constituents[jj].px())),
                        //static_cast<float>(constituents[jj].eta()),
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

    // End Pythia8
    pythia.stat();
    fclose(outbin);
    return 0;
}

// *** define some needed functions ***

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran33(long *idum) {
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0) { 
	if (-(*idum) < 1) *idum=1; 
	else *idum = -(*idum);
	for (j=NTAB+7;j>=0;j--) { 
	  k=(*idum)/IQ1;
	  *idum=IA1*(*idum-k*IQ1)-k*IR1;
	  if (*idum < 0) *idum += IM1;
	  if (j < NTAB) iv[j] = *idum;
	}
	iy=iv[0];
  }
  k=(*idum)/IQ1; 
  *idum=IA1*(*idum-k*IQ1)-k*IR1; 
  if (*idum < 0) *idum += IM1; 
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; 
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV; 
  iy=iv[j]-idum2; 
  iv[j] = *idum; 
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; 
  else return temp;
}






void rotate(double px,double py,double pz,double pr[4], int icc) {
    //     input:  (px,py,pz), (wx,wy,wz), argument (i)
    //     output: new (wx,wy,wz)
    //     if i=1, turn (wx,wy,wz) in the direction (px,py,pz)=>(0,0,E)
    //     if i=-1, turn (wx,wy,wz) in the direction (0,0,E)=>(px,py,pz)

	   
    double wx,wy,wz,E,pt, cosa,sina,cosb,sinb;
    double wx1,wy1,wz1;	   	   

    wx=pr[1];
    wy=pr[2];
    wz=pr[3];

    E=sqrt(px*px+py*py+pz*pz);
    pt=sqrt(px*px+py*py);

    //w = sqrt(wx*wx+wy*wy+wz*wz);

    //  if(pt==0)
    if (pt<1e-6) {
        cosa=1;
        sina=0;
    } 
    else {
        cosa=px/pt;
        sina=py/pt;
    }

    if (E>1e-6) {
        cosb=pz/E;
        sinb=pt/E;
        if (icc==1) {
            wx1=wx*cosb*cosa+wy*cosb*sina-wz*sinb;
    	      wy1=-wx*sina+wy*cosa;
    	      wz1=wx*sinb*cosa+wy*sinb*sina+wz*cosb;
          } else {
    	      wx1=wx*cosa*cosb-wy*sina+wz*cosa*sinb;
    	      wy1=wx*sina*cosb+wy*cosa+wz*sina*sinb;
    	      wz1=-wx*sinb+wz*cosb;
          }
        wx=wx1;
        wy=wy1;
        wz=wz1;
    } else {
        cout << "warning: small E in rotation" << endl;
    }

    pr[1]=wx;
    pr[2]=wy;
    pr[3]=wz;      
}



// get the random position!! --> add by piduan
void get_coordinate(vector<double>& parton_x, vector<double>& parton_y, const string& pathin) {

    ifstream finput(pathin,ifstream::in);
    double x_tep,y_tep;
    if(finput.is_open()) {
      while(finput >> x_tep) {
        finput >> y_tep;

        parton_x.push_back(x_tep);
        parton_y.push_back(y_tep);
      }

      //cout << parton_x.size() << endl;

    }
    else { 
      cout <<" No FOUND coordinate file!" << endl;
    }

}

// Function to recursively find all particles and their second-generation decay products
void findDecayProducts(const Event& event, int particleIndex, vector<int>& decayedParticles) {
    const Particle& particle = event[particleIndex];

    // Add the current particle to the list
    decayedParticles.push_back(particleIndex);

    if (particle.daughter1() == particle.daughter2()) {
        // Recursively find daughters
        if (particle.daughter1() != 0) {
            // Recursively find daughters
            findDecayProducts(event, particle.daughter1(), decayedParticles);
        }
    } else {
        if (particle.daughter1() > particle.daughter2()) {
            // Check if the particle has daughters
            if (particle.daughter1() != 0) {
                // Recursively find daughters
                findDecayProducts(event, particle.daughter1(), decayedParticles);
            }
            if (particle.daughter2() != 0) {
                // Recursively find daughters
                findDecayProducts(event, particle.daughter2(), decayedParticles);
            }
        } else {
            for (int kk = particle.daughter1(); kk <=particle.daughter2(); kk++) {
                if (kk != 0) findDecayProducts(event, kk, decayedParticles);
            }
        }
    }
}

