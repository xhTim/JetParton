#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "Pythia8/Pythia.h"
#include "TMath.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TSystem.h"
#include "TInterpreter.h"

//standard cpp libraries
#include <vector>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>

//namespaces
using namespace std;
using namespace Pythia8;
using namespace fastjet;
Pythia pythia;
Pythia8::ParticleData &particleData = pythia.particleData;

int main(int argv, char* argc[])
{
    gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
    gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");

    int Nevent = atoi(argc[1]);
    int jobnumber = atoi(argc[2]);
    // selection for final particles which are used to reconstruct jet
    double absetamax = 2.4;
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
    
    //Save the parton info

    ifstream inputFile_p("parton_info.dat");
    if (!inputFile_p.is_open()) {
        cerr << "Error opening parton_info.dat!" << endl;
        return 1;
    }
    string line;
    string line2;
    getline(inputFile_p, line); // Read and ignore the first line

    ifstream inputFile_h("hadrons_frag_full.dat");
    if (!inputFile_h.is_open()) {
        cerr << "Error opening hadrons_frag_full.dat!" << endl;
        return 1;
    }
    string lineh;

    int par_pdgid;
    double par_px, par_py, par_pz, par_e, par_x, par_y, par_z, par_t, mid1, mid2;

    // parton stuff
    std::vector<int> parton_pid;
    std::vector<float> parton_px;
    std::vector<float> parton_py;
    std::vector<float> parton_pz;
    std::vector<float> parton_e;
    std::vector<float> parton_x;
    std::vector<float> parton_y;
    std::vector<float> parton_z;
    std::vector<float> parton_t;

    std::vector<float> genpx;
    std::vector<float> genpy;
    std::vector<float> genpz;
    std::vector<float> genm;
    std::vector<int> genpid;
    std::vector<int> genchg;

    // jet stuff
    std::vector<float> genJetPt;
    std::vector<float> genJetEta;
    std::vector<float> genJetPhi;
    std::vector<int> genJetChargedMultiplicity;
    std::vector<std::vector<int>> gendau_chg;
    std::vector<std::vector<int>> gendau_pid;
    std::vector<std::vector<float>> gendau_pt;
    std::vector<std::vector<float>> gendau_eta;
    std::vector<std::vector<float>> gendau_phi;

    TTree *trackTree = new TTree("trackTree", "v1");

    trackTree->Branch("par_pdgid", &parton_pid);
    trackTree->Branch("par_px", &parton_px);
    trackTree->Branch("par_py", &parton_py);
    trackTree->Branch("par_pz", &parton_pz);
    trackTree->Branch("par_e", &parton_e);
    trackTree->Branch("par_x", &parton_x);
    trackTree->Branch("par_y", &parton_y);
    trackTree->Branch("par_z", &parton_z);
    trackTree->Branch("par_t", &parton_t);

    trackTree->Branch("px", &genpx);
    trackTree->Branch("py", &genpy);
    trackTree->Branch("pz", &genpz);
    trackTree->Branch("m", &genm);
    trackTree->Branch("pid", &genpid);
    trackTree->Branch("chg", &genchg);

    trackTree->Branch("genJetEta", &genJetEta);
    trackTree->Branch("genJetPt", &genJetPt);
    trackTree->Branch("genJetPhi", &genJetPhi);
    trackTree->Branch("genJetChargedMultiplicity", &genJetChargedMultiplicity);
    trackTree->Branch("genDau_chg", &gendau_chg);
    trackTree->Branch("genDau_pid", &gendau_pid);
    trackTree->Branch("genDau_pt", &gendau_pt);
    trackTree->Branch("genDau_eta", &gendau_eta);
    trackTree->Branch("genDau_phi", &gendau_phi);

    int njetevent_count = 0;

    //*******************************START EVENT LOOP****************************************
    for (int iev = 0; iev < Nevent; iev++)
    {
	parton_pid.clear();
	parton_px.clear();
	parton_py.clear();
	parton_pz.clear();
	parton_e.clear();
	parton_x.clear();
	parton_y.clear();
	parton_z.clear();
	parton_t.clear();

	genpx.clear();
	genpy.clear();
	genpz.clear();
	genm.clear();
	genpid.clear();
	genchg.clear();

	genJetPt.clear();
	genJetEta.clear();
	genJetPhi.clear();
	genJetChargedMultiplicity.clear();
	gendau_chg.clear();
	gendau_pid.clear();
	gendau_pt.clear();
	gendau_eta.clear();
	gendau_phi.clear();

	// PARTON stuff

	int numPartons;
	if (inputFile_p.eof())
	{
		cout << "end the file" << endl;
		break;
	}
	while (getline(inputFile_p, line2))
	{
		if (line2.empty() || line2[0] == '#') continue;
		break;
	}
	istringstream iss(line2);
	iss >> numPartons;

	if (inputFile_h.eof())
	{
		cout << "end the hadron file" << endl;
		break;
	}
	int idx, nHadrons, z3, z4;
	while (getline(inputFile_h, lineh))
	{
		istringstream issh(lineh);
		issh >> idx >> nHadrons >> z3 >> z4;
		if (nHadrons > 0) break;

		// Skip 0-hadron events
		for (int i = 0; i < numPartons; ++i)
			inputFile_p >> par_pdgid >> par_px >> par_py >> par_pz >> par_e >> par_x >> par_y >> par_z >> par_t >> mid1 >> mid2;
		while (getline(inputFile_p, line2))
		{
			if (line2.empty() || line2[0] == '#') continue;
			break;
		}
		istringstream iss(line2);
		iss >> numPartons;

		if (inputFile_p.eof() || inputFile_h.eof()) 
            	{
                	cout << "Skip until EOF" << endl;
                	break;
            	}
	}

	for (int i = 0; i < numPartons; ++i)
	{
		inputFile_p >> par_pdgid >> par_px >> par_py >> par_pz >> par_e >> par_x >> par_y >> par_z >> par_t >> mid1 >> mid2;

		parton_pid.push_back(par_pdgid);
		parton_px.push_back(par_px);
		parton_py.push_back(par_py);
		parton_pz.push_back(par_pz);
		parton_e.push_back(par_e);
		parton_x.push_back(par_x);
		parton_y.push_back(par_y);
		parton_z.push_back(par_z);
		parton_t.push_back(par_t);
	}

	// PARTICLE stuff

	// test infile status
	if (feof(infile))
	{
		event_loop_flag = 0;
		cout << " End the event loop ~~~ " << endl;
		break;
	}

	fscanf(infile, "%s %d\n", stemp1, &total_number_of_particles);
	cout << "nHadrons: " << nHadrons << "\ttotal_number_of_particles: " << total_number_of_particles << endl;
	input_particles.clear();
	//******************START daughter particles LOOP********************************
	for (auto i = 0; i < total_number_of_particles; i++)
	{
		if (feof(infile))
		{
			event_loop_flag = 0;
			cout << " End the event loop, and drop last event ~~~ " << endl;
			break;
		}
		fscanf(infile, "%d %lf %lf %lf %lf %lf\n", &pid, &mass,
		       &energy, &px, &py, &pz);
		if (isnan(energy) || isnan(px) || isnan(py) || isnan(pz))
			continue;
		fastjet::PseudoJet particle = PseudoJet(px, py, pz, energy);
		particle.set_user_index(pid);
		input_particles.push_back(particle);

		genpx.push_back(px);
		genpy.push_back(py);
		genpz.push_back(pz);
		genm.push_back(mass);
		genpid.push_back(pid);
		genchg.push_back(particleData.charge(pid));
	} //**************************END daughter partcicles LOOP**********************************

	if (event_loop_flag == 0)
	{
		cout << " End the event loop and drop last event ~~~ " << endl;
		break;
	}

	count_event_number++;

	// Then do the jet finding
	// fastjet::Selector particle_selector = fastjet::SelectorAbsEtaMax(absetamax) && fastjet::SelectorPtMin( particle_ptmin );
	fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R_jet);
	// select jet
	// fastjet::Selector jet_selector = fastjet::SelectorAbsEtaMax( jet_absetamax ) && fastjet::SelectorPtMin( jet_ptmin );
	// input_particles = particle_selector(input_particles);
	fastjet::ClusterSequence clust_seq(input_particles, jet_def);
	// get the resulting jets ordered in pt
	vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets());

	// Select the jet pT and output the selected events, rotate the jet at pz direction
	//*************************************START JET LOOP*****************************************
	for (unsigned int i = 0; i < inclusive_jets.size(); i++)
	{
		if (inclusive_jets[i].pt() < jet_ptmin)
			continue;
		std::vector<float> tmp_pt;
		std::vector<float> tmp_eta;
		std::vector<float> tmp_phi;
		std::vector<int> tmp_chg;
		std::vector<int> tmp_pid;
		vector<PseudoJet> constituents = inclusive_jets[i].constituents();
		int chMult = 0;
		//***********************START constituents LOOP*********************************
		for (unsigned j = 0; j < constituents.size(); j++)
		{
			if (particleData.charge(constituents[j].user_index()))
			{
				chMult++;
			}
			tmp_pt.push_back(constituents[j].pt());
			tmp_eta.push_back(constituents[j].eta());
			tmp_phi.push_back(constituents[j].phi());
			tmp_chg.push_back(particleData.charge(constituents[j].user_index()));
			tmp_pid.push_back(constituents[j].user_index());
		}
		//***********************END constituents LOOP*********************************

		genJetPt.push_back(inclusive_jets[i].pt());
		genJetEta.push_back(inclusive_jets[i].eta());
		genJetPhi.push_back(inclusive_jets[i].phi());
		genJetChargedMultiplicity.push_back(chMult);
		gendau_pt.push_back(tmp_pt);
		gendau_eta.push_back(tmp_eta);
		gendau_phi.push_back(tmp_phi);
		gendau_chg.push_back(tmp_chg);
		gendau_pid.push_back(tmp_pid);
	} //************************************END JET LOOP******************************************
	    trackTree->Fill();




    }//****************************************END EVENT LOOP***************************************

    fclose(infile);
	inputFile_p.close();
    TFile * fout = TFile::Open( Form("/eos/cms/store/group/phys_heavyions/huangxi/PC/pp_parton_cascade_%d.root",jobnumber) ,"recreate");
    trackTree->Write();
    fout->Close();
    return 0;
}

