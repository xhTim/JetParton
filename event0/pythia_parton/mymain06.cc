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
    //string random_str = string(argc[1]);
    int Spatial_mode = 1; // 0: use the accumulant method; 1: uses the free-streaming to formation time; 
    // Set up the Pythia8 configuration.
    Pythia pythia;
    string random_str = string(argc[2]);
    int totalevent = atoi(argc[1]);
    /*time_t now = time(NULL);
    int offset = atoi(argc[2]);
    int seed = static_cast<int>(now) + offset;*/
    pythia.readString("Random:setSeed = on");
    //pythia.readString("Random:seed = 0");  //92549 792 458 871
    pythia.readString("Random:seed = " + random_str);
    // pythia.readString("Random:seed = " + seed);
    pythia.readString("Beams:eCM = 13000.0");//5020,2760 GeV
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212"); //2212 is proton; 1000791970 for 197Au; 1000822080 for 208Pb;
    // Standard settings
    pythia.readString("HardQCD:all = on");
    //pythia.readString("Tune:pp=18");
    pythia.readString("PhaseSpace:pTHatMin = 500.0");
    pythia.readString("PhaseSpace:pTHatMax = -1.");

    pythia.readString("Tune:pp=14");
    //pythia.readString("PDF:pSet = 20");
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
    pythia.readString("SigmaTotal:sigmaEl = 22.08");
    pythia.readString("SigmaTotal:sigmaTot = 101.037");
    pythia.readString("PDF:pSet=LHAPDF6:NNPDF31_nnlo_as_0118");
    
    pythia.readString("HadronLevel:Hadronize = off");
    // Initialize Pythia
    if (!pythia.init()) {
        std::cerr << "Pythia initialization failed!" << std::endl;
        return 1;
    }

    // Output the parton information
    string output_filename;
    output_filename = "parton_info.dat";
    ofstream output_parton(output_filename.c_str());  //-----zhao

    if( ! output_parton.is_open() ) {
        cout << "cannot open output file:"<< endl
             << output_filename << endl;
        return -1;
    }
    
    // parameter setting
    const double jet_absetamax = 1.6;
    const double jet_ptmin = 450.0;
    
    //fastjet setting
    // create a jet definition: 
    // a jet algorithm with a given radius parameter
    double R = 0.8;
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);

   // select jet
    fastjet::Selector jet_selector = fastjet::SelectorAbsEtaMax( jet_absetamax ) && fastjet::SelectorPtMin( jet_ptmin );
    
    // selection for final particles which are used to reconstruct jet
    double absetamax = R+jet_absetamax;
    double particle_ptmin = 0.3;
    fastjet::Selector particle_selector = fastjet::SelectorAbsEtaMax(absetamax) && fastjet::SelectorPtMin( particle_ptmin );
    
    // Number of events to generate
    int numEvent = 0;
    output_parton << "# pdgid  px  py  pz  energy  x  y  z  t" << std::endl;
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
            //if (particle.isFinal() &&
            //    (fabs(particle.id()) == 1 || fabs(particle.id()) == 2 ||  fabs(particle.id()) == 3 || particle.id() == 21)) {
            if (particle.isFinal() && particle.isParton()) {
                Nparton++;
                fastjet::PseudoJet particlefastjet = PseudoJet(particle.px(), particle.py(), particle.pz(), particle.e());
                particlefastjet.set_user_index(particle.id());
                input_particles.push_back(particlefastjet);
            }
        }
        input_particles = particle_selector(input_particles);
        // run the jet clustering with the above jet definition
        //----------------------------------------------------------
        fastjet::ClusterSequence clust_seq(input_particles, jet_def);

        // get the resulting jets ordered in pt
        //----------------------------------------------------------
        vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt( jet_selector(clust_seq.inclusive_jets()) );
        
        if (inclusive_jets.size() == 0 || inclusive_jets[0].pt() < 500.){
            continue;
        } else {
            numEvent +=1;
            if ( numEvent > totalevent) break;
        }
        
        output_parton << "# event id " << numEvent  <<  ", Number_of_parton  "  << std::endl;
        output_parton << Nparton << std::endl;
        for (int j = 0; j < event.size(); ++j) {
            Particle& particle = event[j];
            //if (particle.isFinal() &&
            //    (fabs(particle.id()) == 1 || fabs(particle.id()) == 2 ||  fabs(particle.id()) == 3 || particle.id() == 21)) {
            if ( particle.isFinal() && particle.isParton() ) {
                double position[3] = {0.0}; // spatial information of partons (x, y, z)
                position[0] = 0.0;
                position[1] = 0.0;
                position[2] = 0.0;
                //......formation time + ......
                //int ishower = 0;
                double p0[4] = {0.0}; // particle's four momentum
                double p4[4] = {0.0}; // mother's four momentum
		double qt, time_step;
		double timeplus = 0.0;
                int IDmom1, IDmom2;
                int timebreaker = 0;
                int IDmom0 = j;
                //... start to calculate the formation and the spatial information of partons ...
                while (timebreaker == 0) {
                    int IDiii = IDmom0;
                    if (abs(pythia.event[IDiii].status())==23 || abs(pythia.event[IDiii].status())==21 ||
                        abs(pythia.event[IDiii].status())==12) timebreaker=1;			
                    IDmom1 = pythia.event[IDiii].mother1();
                    IDmom2 = pythia.event[IDiii].mother2();
                    if (IDmom1==IDmom2 && IDmom1==0) timebreaker=1;
                    if (IDmom1==IDmom2 && IDmom1>0) IDmom0=IDmom1;
                    if (IDmom1>0 && IDmom2==0) {
                        IDmom0=IDmom1;
                        double IDdaughter1 = pythia.event[IDmom0].daughter1();
                        double IDdaughter2 = pythia.event[IDmom0].daughter2();
                        if (IDdaughter1 != IDdaughter2 && IDdaughter1>0 && IDdaughter2>0) {
                            p4[0] = pythia.event[IDdaughter1].e()+pythia.event[IDdaughter2].e();
                            p4[1] = pythia.event[IDdaughter1].px()+pythia.event[IDdaughter2].px();
                            p4[2] = pythia.event[IDdaughter1].py()+pythia.event[IDdaughter2].py();
                            p4[3] = pythia.event[IDdaughter1].pz()+pythia.event[IDdaughter2].pz();
                            double x_split = pythia.event[IDiii].e()/p4[0];
                            if (x_split>1) x_split=1.0/x_split;
                            p0[0] = pythia.event[IDiii].e();
                            p0[1] = pythia.event[IDiii].px();
                            p0[2] = pythia.event[IDiii].py();
                            p0[3] = pythia.event[IDiii].pz(); 
                            //double pt_daughter = sqrt(pow(p0[1],2)+pow(p0[2],2));
                            //double pt_mother = sqrt(pow(p4[1],2)+pow(p4[2],2));			  

                            rotate(p4[1],p4[2],p4[3],p0,1);
                            qt = sqrt(pow(p0[1],2)+pow(p0[2],2));
                            rotate(p4[1],p4[2],p4[3],p0,-1);
                            double kt_daughter=qt;
                            if (x_split<0.5) {
                                if(kt_daughter > 0.0001) {
                                    time_step = 2.0*p4[0]*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1]/p4[0];
                                    position[1] = position[1] + time_step * p4[2]/p4[0];
                                    position[2] = position[2] + time_step * p4[3]/p4[0];
                                }
                            } else {
                                if(kt_daughter > 0.0001) {
                                    time_step = 1/p4[0];
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1]/p4[0];
                                    position[1] = position[1] + time_step * p4[2]/p4[0];
                                    position[2] = position[2] + time_step * p4[3]/p4[0];
                                }
                            }
                        }
                    }//if(IDmom1>0 && IDmom2==0)

                    if (IDmom1<IDmom2 && IDmom1>0 && IDmom2>0) {
                        if(pythia.event[IDmom1].e()>pythia.event[IDmom2].e()) {
                          IDmom0=IDmom1;			
                        }
                        if(pythia.event[IDmom1].e()<=pythia.event[IDmom2].e()) {
                          IDmom0=IDmom2;			
                        }
                        double IDdaughter1=pythia.event[IDmom0].daughter1();
                        double IDdaughter2=pythia.event[IDmom0].daughter2();
                        if (IDdaughter1 != IDdaughter2 && IDdaughter1>0 && IDdaughter2>0) {
                            p4[0]=pythia.event[IDdaughter1].e()+pythia.event[IDdaughter2].e();
                            p4[1]=pythia.event[IDdaughter1].px()+pythia.event[IDdaughter2].px();
                            p4[2]=pythia.event[IDdaughter1].py()+pythia.event[IDdaughter2].py();
                            p4[3]=pythia.event[IDdaughter1].pz()+pythia.event[IDdaughter2].pz();
                            double x_split = pythia.event[IDiii].e()/p4[0];
                            if(x_split>1) x_split=1.0/x_split;

                            p0[0]=pythia.event[IDiii].e();
                            p0[1]=pythia.event[IDiii].px();
                            p0[2]=pythia.event[IDiii].py();
                            p0[3]=pythia.event[IDiii].pz(); 
                            rotate(p4[1],p4[2],p4[3],p0,1);
                            qt=sqrt(pow(p0[1],2)+pow(p0[2],2));
                            rotate(p4[1],p4[2],p4[3],p0,-1);
                            double kt_daughter=qt;

                            if (x_split<0.5) {
			        if (kt_daughter > 0.0001) {
                                    time_step = 2.0*p4[0]*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1]/p4[0];
                                    position[1] = position[1] + time_step * p4[2]/p4[0];
                                    position[2] = position[2] + time_step * p4[3]/p4[0];
                                    
                                }
                            } else {
                                if(kt_daughter > 0.0001) {
                                    time_step = 1/p4[0];
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1]/p4[0];
                                    position[1] = position[1] + time_step * p4[2]/p4[0];
                                    position[2] = position[2] + time_step * p4[3]/p4[0];
                                }
                            }
                        }
			
			if ( (IDdaughter1 > 0 && IDdaughter2 == 0) || (IDdaughter1 > 0 && IDdaughter2 == IDdaughter1)) {
                            p4[0] = event[IDdaughter1].e();
                            p4[1] = event[IDdaughter1].px();
                            p4[2] = event[IDdaughter1].py();
                            p4[3] = event[IDdaughter1].pz();
                            double x_split = event[IDiii].e()/p4[0];
                            if (x_split>1) x_split=1.0/x_split; // revise the mother and daughter
                            p0[0] = event[IDiii].e();
                            p0[1] = event[IDiii].px();
                            p0[2] = event[IDiii].py();
                            p0[3] = event[IDiii].pz(); 
                            //double pt_daughter=sqrt(pow(p0[1],2)+pow(p0[2],2));
                            //double pt_mother=sqrt(pow(p4[1],2)+pow(p4[2],2));			 
                            rotate(p4[1],p4[2],p4[3],p0,1); // rotate into the 
                            qt = sqrt(pow(p0[1],2)+pow(p0[2],2));
                            rotate(p4[1],p4[2],p4[3],p0,-1);
                            double kt_daughter = qt;
                            //double Q2 = 1.0/(x_split*(1-x_split)/pow(kt_daughter,2));
                            if(x_split<0.5){
			        if (kt_daughter > 0.0001) {						
                                    time_step = 2.0*p4[0]*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1]/p4[0];
                                    position[1] = position[1] + time_step * p4[2]/p4[0];
                                    position[2] = position[2] + time_step * p4[3]/p4[0];					
			        }
                            } else {
                                if(kt_daughter > 0.0001) {
                                    time_step = 1/p4[0];
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1]/p4[0];
                                    position[1] = position[1] + time_step * p4[2]/p4[0];
                                    position[2] = position[2] + time_step * p4[3]/p4[0];
                                }
                            }
                        }
                        
                    }//if(IDmom1<IDmom2 && IDmom1>0 && IDmom2>0)

                    if (IDmom1>IDmom2 && IDmom1>0 && IDmom2>0) {
                        if (pythia.event[IDmom1].e()>pythia.event[IDmom2].e()) {
                            IDmom0=IDmom1;			
                        }
                        if (pythia.event[IDmom1].e()<=pythia.event[IDmom2].e()) {
                            IDmom0=IDmom2;			
                        }			
                        double IDdaughter1=pythia.event[IDmom0].daughter1();
                        double IDdaughter2=pythia.event[IDmom0].daughter2();
                        if (IDdaughter1 != IDdaughter2 && IDdaughter1>0 && IDdaughter2>0) {
                            p4[0]=pythia.event[IDdaughter1].e()+pythia.event[IDdaughter2].e();
                            p4[1]=pythia.event[IDdaughter1].px()+pythia.event[IDdaughter2].px();
                            p4[2]=pythia.event[IDdaughter1].py()+pythia.event[IDdaughter2].py();
                            p4[3]=pythia.event[IDdaughter1].pz()+pythia.event[IDdaughter2].pz();
                            double x_split=pythia.event[IDiii].e()/p4[0];
                            if(x_split>1) x_split=1.0/x_split;
                            p0[0]=pythia.event[IDiii].e();
                            p0[1]=pythia.event[IDiii].px();
                            p0[2]=pythia.event[IDiii].py();
                            p0[3]=pythia.event[IDiii].pz(); 
                            rotate(p4[1],p4[2],p4[3],p0,1);
                            qt=sqrt(pow(p0[1],2)+pow(p0[2],2));
                            rotate(p4[1],p4[2],p4[3],p0,-1);
                            double kt_daughter=qt;
                            if(x_split<0.5){
			        if (kt_daughter > 0.0001) {						
                                    time_step = 2.0*p4[0]*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1]/p4[0];
                                    position[1] = position[1] + time_step * p4[2]/p4[0];
                                    position[2] = position[2] + time_step * p4[3]/p4[0];					
			        }
                            } else {
                                if(kt_daughter > 0.0001) {
                                    time_step = 1/p4[0];
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1]/p4[0];
                                    position[1] = position[1] + time_step * p4[2]/p4[0];
                                    position[2] = position[2] + time_step * p4[3]/p4[0];
                                }
                            }
                        }
                        
                        if ( (IDdaughter1 > 0 && IDdaughter2 == 0) || (IDdaughter1 > 0 && IDdaughter2 == IDdaughter1)) {
                            p4[0] = event[IDdaughter1].e();
                            p4[1] = event[IDdaughter1].px();
                            p4[2] = event[IDdaughter1].py();
                            p4[3] = event[IDdaughter1].pz();
                            double x_split = event[IDiii].e()/p4[0];
                            if (x_split>1) x_split=1.0/x_split; // revise the mother and daughter
                            p0[0] = event[IDiii].e();
                            p0[1] = event[IDiii].px();
                            p0[2] = event[IDiii].py();
                            p0[3] = event[IDiii].pz(); 
                            //double pt_daughter=sqrt(pow(p0[1],2)+pow(p0[2],2));
                            //double pt_mother=sqrt(pow(p4[1],2)+pow(p4[2],2));			 
                            rotate(p4[1],p4[2],p4[3],p0,1); // rotate into the 
                            qt = sqrt(pow(p0[1],2)+pow(p0[2],2));
                            rotate(p4[1],p4[2],p4[3],p0,-1);
                            double kt_daughter = qt;
                            //double Q2 = 1.0/(x_split*(1-x_split)/pow(kt_daughter,2));
                            if(x_split<0.5){
			        if (kt_daughter > 0.0001) {						
                                    time_step = 2.0*p4[0]*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1]/p4[0];
                                    position[1] = position[1] + time_step * p4[2]/p4[0];
                                    position[2] = position[2] + time_step * p4[3]/p4[0];					
			        }
                            } else {
                                if(kt_daughter > 0.0001) {
                                    time_step = 1/p4[0];
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1]/p4[0];
                                    position[1] = position[1] + time_step * p4[2]/p4[0];
                                    position[2] = position[2] + time_step * p4[3]/p4[0];
                                }
                            }
                        }
                        
                    }//if(IDmom1>IDmom2 && IDmom1>0 && IDmom2>0)
                    
                    if (IDmom1==IDmom2 && IDmom1>0) {
                        IDmom0=IDmom1;
                        int IDdaughter1 = event[IDmom0].daughter1();
                        int IDdaughter2 = event[IDmom0].daughter2();
                        
                        if (IDdaughter1 != IDdaughter2 && IDdaughter1>0 && IDdaughter2>0) {
                            p4[0] = event[IDdaughter1].e()+event[IDdaughter2].e();
                            p4[1] = event[IDdaughter1].px()+event[IDdaughter2].px();
                            p4[2] = event[IDdaughter1].py()+event[IDdaughter2].py();
                            p4[3] = event[IDdaughter1].pz()+event[IDdaughter2].pz();
                            double x_split = event[IDiii].e()/p4[0];
                            if (x_split>1) x_split=1.0/x_split; // revise the mother and daughter
                            p0[0] = event[IDiii].e();
                            p0[1] = event[IDiii].px();
                            p0[2] = event[IDiii].py();
                            p0[3] = event[IDiii].pz(); 
                            rotate(p4[1],p4[2],p4[3],p0,1); // rotate into the 
                            qt = sqrt(pow(p0[1],2)+pow(p0[2],2));
                            rotate(p4[1],p4[2],p4[3],p0,-1);
                            double kt_daughter = qt;
                            if(x_split<0.5){
			        if (kt_daughter > 0.0001) {						
                                    time_step = 2.0*p4[0]*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1]/p4[0];
                                    position[1] = position[1] + time_step * p4[2]/p4[0];
                                    position[2] = position[2] + time_step * p4[3]/p4[0];					
			        }
                            } else {
                                if(kt_daughter > 0.0001) {
                                    time_step = 1/p4[0];
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1]/p4[0];
                                    position[1] = position[1] + time_step * p4[2]/p4[0];
                                    position[2] = position[2] + time_step * p4[3]/p4[0];
                                }
                            }
                        }
                        
                        if ( (IDdaughter1 > 0 && IDdaughter2 == 0) || (IDdaughter1 > 0 && IDdaughter2 == IDdaughter1)) {
                            p4[0] = event[IDdaughter1].e();
                            p4[1] = event[IDdaughter1].px();
                            p4[2] = event[IDdaughter1].py();
                            p4[3] = event[IDdaughter1].pz();
                            double x_split = event[IDiii].e()/p4[0];
                            if (x_split>1) x_split=1.0/x_split; // revise the mother and daughter
                            p0[0] = event[IDiii].e();
                            p0[1] = event[IDiii].px();
                            p0[2] = event[IDiii].py();
                            p0[3] = event[IDiii].pz(); 
                            rotate(p4[1],p4[2],p4[3],p0,1); // rotate into the 
                            qt = sqrt(pow(p0[1],2)+pow(p0[2],2));
                            rotate(p4[1],p4[2],p4[3],p0,-1);
                            double kt_daughter = qt;
                            //if (kt_daughter <= 0.0001) kt_daughter = 0.0001;
                            if(x_split<0.5){
			        if (kt_daughter > 0.0001) {						
                                    time_step = 2.0*p4[0]*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1]/p4[0];
                                    position[1] = position[1] + time_step * p4[2]/p4[0];
                                    position[2] = position[2] + time_step * p4[3]/p4[0];					
			        }
                            } else {
                                if(kt_daughter > 0.0001) {
                                    time_step = 1/p4[0];
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1]/p4[0];
                                    position[1] = position[1] + time_step * p4[2]/p4[0];
                                    position[2] = position[2] + time_step * p4[3]/p4[0];
                                }
                            }
                        }
                    }
                    
                } //while(timebreaker == 0)

                // Extract relevant information
                int pdgId = particle.id();
                /*
                int status = particle.status();
                int mother1 = particle.mother1();
                int mother2 = particle.mother2();
                int daughter1 = particle.daughter1();
                int daughter2 = particle.daughter2();
                */
                double px = particle.px();
                double py = particle.py();
                double pz = particle.pz();
                double energy = particle.e();
                if (Spatial_mode == 1) {
                    position[0] = 0.0 + timeplus * px / energy;
                    position[1] = 0.0 + timeplus * py / energy;
                    position[2] = 0.0 + timeplus * pz / energy;
                }
                output_parton << pdgId << "  " 
                              << px << "  " << py << "  " << pz << "  " << energy <<"  "
                              << position[0] << "  "<< position[1] << "  " << position[2] << "  " << timeplus << "  "
                              << particle.col() << "  " << particle.acol()
                              << std::endl;
            }
        }
    }

    // End Pythia8
    pythia.stat();
    output_parton.close();
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


