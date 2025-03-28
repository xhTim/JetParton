#include <assert.h>
#include <time.h>
#include <sstream>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string.h>
#include <algorithm>
#include <vector>
#include "Pythia8/Pythia.h" 
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh" 

#define PI 3.1415926
#define QQ0 1.0

//*****************************
//  The Pythia 8 string fragmentation
//  copywrite Wenbin Zhao, 2020, 12.24
//  When you use this code, please cite our papers:
//      W.~Zhao, etc. al.[arXiv:2103.14657 [hep-ph]].                                               
//      W.~Zhao, etc. al.Phys. Rev. Lett.125, (2020) no.7, 072301.
//*****************************

using namespace Pythia8;
using namespace fastjet;

int searchmin(double*p, int *q,int*s,int len);
int searchmax(double*p, int *q,int len);
int findMax(std::vector<int>& vec);
int searchmin2(double*p,int len);
bool findarray(int*p,int len,int val);
char infiles[128];

//===========================================================================================
int main(int argv, char* argc[])
{
    //string random_str = string(argc[1]);
    int n_event = atoi(argc[1]);
    int Nchcut = atoi(argc[2]);
    bool DO_Colorless_frag = false;
    string output_filename2;
    output_filename2 = "hadrons_frag1.dat";// output files of final hadrons
    //string ramdomseed_str = "Random:seed = "+random_str;
    ofstream output2(output_filename2.c_str()); 
    if (!output2.is_open() ) {
        cout << "cannot open output file:"<< endl
         << output_filename2 << endl;
        return -1;
    }
    string output_filename;
    output_filename = "selected_parton/parton_info.dat";// output files of final hadrons
    //string ramdomseed_str = "Random:seed = "+random_str;
    ofstream output(output_filename.c_str()); 
    if (!output.is_open() ) {
        cout << "cannot open output file:"<< endl
         << output_filename << endl;
        return -1;
    }
    output << "# pdgid  px  py  pz  energy  x  y  z  t" << std::endl;
     /*
	output2<<"OSC1997A"<<endl;
	output2<<"final_id_p_x"<<endl;
	output2<<" 3DHydro       1.1  (197,    79)+(197,    79)  eqsp  0.1000E+03         1"<<endl;
	*/
    FILE* infile1;
    infile1 = fopen("parton_info.dat","r");

    Pythia pythia;
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = 0");  
    pythia.readString("Beams:eCM = 5020");
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212"); 
    // Standard settings
    pythia.readString("HardQCD:all = on");
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("ColourReconnection:reconnect = on");
    pythia.readString("ColourReconnection:mode = 1");
    pythia.readString("MultipartonInteractions:pT0Ref = 2.15");     
    pythia.readString("ColourReconnection:allowDoubleJunRem = off");
    pythia.readString("ColourReconnection:junctionCorrection = 1.15");
    pythia.readString("ColourReconnection:timeDilationMode = 3");    
    pythia.readString("ColourReconnection:timeDilationPar = 0.18");  
    pythia.readString("ProcessLevel:all = off"); // The trick!
    // CMS CP5 setting
    pythia.readString("Tune:pp=14");
    pythia.readString("Tune:ee=7");
    pythia.readString("MultipartonInteractions:ecmPow=0.03344");
    pythia.readString("MultipartonInteractions:bProfile=2");
    pythia.readString("MultipartonInteractions:pT0Ref=1.407");
    pythia.readString("MultipartonInteractions:coreRadius=0.6671");
    pythia.readString("MultipartonInteractions:coreFraction=0.4281");
    pythia.readString("ColourReconnection:range=4.881");
    pythia.readString("SigmaTotal:zeroAXB=off");
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
    pythia.readString(" StringFlav:probStoUD=0.50");
    //general CMS settings
    pythia.readString("Check:epTolErr = 0.01");
    pythia.readString("Beams:setProductionScalesFromLHEF = off");
    pythia.readString("ParticleDecays:limitTau0 = on");
    pythia.readString("ParticleDecays:tau0Max = 10");
    pythia.readString("ParticleDecays:allowPhotonRadiation = on");
    pythia.readString("111:mayDecay = off");//pion0
    pythia.readString("HadronLevel:Decay = on");
    pythia.readString("HadronLevel:Hadronize = on");
    pythia.init();

    double hbarc = 0.19732;
    double c_px, c_py ,c_pz, c_e, c_energy,c_x,c_y,c_z,c_t;
    int c_id,tt, c_col, c_acol;
    int acol_ip[5000]={0};
    int pp_collision=0;    //used for total cross section
    int  ccbar_num=0;
    double cmeson_px,cmeson_py,cmeson_pz, cmeson_P,tempp;
    double cbar_meson_px,cbar_meson_py,cbar_meson_pz, cbar_meson_P,pt_square,cbar_meson_energy;
    double pt,mass,amid,Qmid;
    int mid,ie,II,status,col,acol,Npart,NN,Ntotal,PAosi,simble,mmaxindex;
    int mid1, mid2, mid3, mid4, midd5;
    double midd;
    int idp[10000]={0},idpo[10000]={0};
    double pxpo[10000]={0.0},pypo[10000]={0.0},pzpo[10000]={0.0},epo[10000]={0.0},ptpo[10000]={0.0},xxpo[10000]={0.0},yypo[10000]={0.0},zzpo[10000]={0.0},ttpo[10000]={0.0},phio[10000]={0.0},etao[10000]={0.0};//,mass[10000]={0.0};
    double pxp[10000]={0.0},pyp[10000]={0.0},pzp[10000]={0.0},ep[10000]={0.0},ptp[10000]={0.0},xxp[10000]={0.0},yyp[10000]={0.0},zzp[10000]={0.0},ttp[10000]={0.0},phi[10000]={0.0},distance[10000]={0.0},dsting[100][1000]={0.0},Qscale[10000]={0.0};//,mass[10000]={0.0};
    int nncol[10000]={0},aacol[10000]={0},index[10000]={0},qindex[10000]={0},aqindex[10000]={0},used[10000]={0},pair[1000]={0},apair[1000]={0},gindex[1000]={0},strings[100][1000]={0},Snum[1000]={0};//,nncol_mid[10000]={0},aacol_mid[10000]={0};
    // event loop
    char header[256];
    fgets(header, sizeof(header), infile1);
    std::vector<int> Nchevent;
    int event_loop_flag = 1;
    for (int iEvent=0; iEvent<n_event; iEvent++) {
        if(feof(infile1)) {
            event_loop_flag = 0;
            cout << " End the event loop ~~~ " << endl;
            break;
        }
        Nchevent.clear();     
        // Read and ignore the first line starting with '#'
        
        char header2[256];
        fgets(header2, sizeof(header2), infile1);
        
        fscanf(infile1,"%d\n",&Npart);
        cout << iEvent << " " << Npart << endl;
        if (Npart==0 ) {output2 << "         " << iEvent <<"          " << 0  << endl; continue;}
        int Nquark=0;
        int Naquark =0;
        int Ngluon=0;
        int Npair=0;
        for (int ll=0;ll<Npart;ll++) {
                if(event_loop_flag == 0) {
                    cout << " End the event loop and drop last event ~~~ " << endl;
                    break;
                }
                fscanf(infile1, "%d %lf %lf %lf %lf %lf %lf %lf %lf %d %d\n",
                       &c_id, &c_px, &c_py, &c_pz, &c_energy, &c_x, &c_y, &c_z, &c_t, &c_col, &c_acol);
                //cout << c_id << endl;
                Qmid = 0.0; // The scale for the parton shower. 
                idpo[ll]=c_id;
                if(c_id==21){gindex[Ngluon]=ll;Ngluon++;}
                pxpo[ll]=c_px;
                pypo[ll]=c_py;
                pzpo[ll]=c_pz;
                ptpo[ll]=c_px*c_px+c_py*c_py;
                double pmg=sqrt(c_px*c_px+c_py*c_py+c_pz*c_pz);
                epo[ll] = c_energy;//sqrt(pmg * pmg + c_energy*c_energy);
                xxpo[ll]=c_x;
                yypo[ll]=c_y;
                zzpo[ll]=c_z;
                ttpo[ll]=c_t;
                Qscale[ll]=sqrt(Qmid);
                double aamid=0.5*log((pmg+c_pz)/(pmg-c_pz));
                etao[ll]=aamid;
                phio[ll]=atan2(c_py,c_px);
                index[ll]=ll;
                used[ll]=0;
                nncol[ll] = c_col;
                aacol[ll] = c_acol;
        }
       if (DO_Colorless_frag) {
        mmaxindex=searchmax(ptpo,index,Npart);// get the leading parton
        //**** all partons are gluon ****
        // add fictive quark anti-quark
        if (Npart>0&&Ngluon==Npart) {
            idpo[Npart]=-3;pxpo[Npart]=0.10;pypo[Npart]=0.20;pzpo[Npart]=100000.0;
            epo[Npart]=sqrt(pzpo[Npart]*pzpo[Npart]+pypo[Npart]*pypo[Npart]+pxpo[Npart]*pxpo[Npart]+0.5*0.5);
            xxpo[Npart]=0.0;yypo[Npart]=0.0;zzpo[Npart]=0.0;ttpo[Npart]=0.0;phio[Npart]=atan2(pypo[Npart],pxpo[Npart]);index[Npart]=Npart;
            used[Npart]=0;
            etao[Npart]=10000.0;
            Npart=Npart+1;
            idpo[Npart]=3;pxpo[Npart]=0.10;pypo[Npart]=0.20;pzpo[Npart]=-100000.0;
            epo[Npart]=sqrt(pzpo[Npart]*pzpo[Npart]+pypo[Npart]*pypo[Npart]+pxpo[Npart]*pxpo[Npart]+0.5*0.5);
            xxpo[Npart]=0.0;yypo[Npart]=0.0;zzpo[Npart]=0.0;ttpo[Npart]=0.0;phio[Npart]=atan2(pypo[Npart],pxpo[Npart]);index[Npart]=Npart;
            used[Npart]=0;
            etao[Npart]=-10000.0;
            Npart=Npart+1;
        }

        // ****************** give the col and acol to gluon and quarks *********************//
        int Nmax=101;
        int add=0;
        int jj,Nmin;
        //*** find the quark pairs ***
        for (int ll=0;ll<Npart;ll++) {
            if ((idpo[ll]<=5)&((idpo[ll])>0)) {
                qindex[Nquark]=ll;Nquark++;
            }
            if ((idpo[ll]>-5)&((idpo[ll])<0)) {
                aqindex[Naquark]=ll;Naquark++;
            }
        }
        int excess=abs(Nquark-Naquark);
        if (excess>0) {
            int ppp=excess % 2;
            int ecc=excess/2;   
            if (Nquark>Naquark) {
                if (ecc>0) {
                    for (int gg=0; gg<ecc; gg++) {
                        idpo[qindex[Nquark-1]]=-1*idpo[qindex[Nquark-1]];
                        aqindex[Naquark]=qindex[Nquark-1];
                        Naquark++;
                        Nquark--;
                    }
                }
                if (ppp==1) {
                    idpo[Npart]=-3; pxpo[Npart]=0.20; pypo[Npart]=0.20;
                    pzpo[Npart]=10000.20;epo[Npart]=sqrt(pzpo[Npart]*pzpo[Npart]+pypo[Npart]*pypo[Npart]+pxpo[Npart]*pxpo[Npart]+0.5*0.5);
                    xxpo[Npart]=0.0;yypo[Npart]=0.0;zzpo[Npart]=0.0;ttpo[Npart]=0.0;phio[Npart]=atan2(pypo[Npart],pxpo[Npart]);index[Npart]=Npart; 
                    etao[Npart]=1.0;
                    used[Npart]=0;
                    aqindex[Naquark]=Npart;
                    Npart++;
                    Naquark++;
                }
            }
            if (Naquark>Nquark) {
                if (ecc>0) {
                    for (int gg=0; gg<ecc; gg++) {
                        idpo[aqindex[Naquark-1]]=-1*idpo[aqindex[Naquark-1]];
                        qindex[Nquark]=aqindex[Naquark-1];
                        Nquark++;
                        Naquark--;
                    }
                }
                if (ppp==1) {
                    idpo[Npart]=3; pxpo[Npart]=0.20; pypo[Npart]=0.20;
                    pzpo[Npart]=100000.20;epo[Npart]=sqrt(pzpo[Npart]*pzpo[Npart]+pypo[Npart]*pypo[Npart]+pxpo[Npart]*pxpo[Npart]+0.5*0.5);
                    xxpo[Npart]=0.0;yypo[Npart]=0.0;zzpo[Npart]=0.0;ttpo[Npart]=0.0;phio[Npart]=atan2(pypo[Npart],pxpo[Npart]);index[Npart]=Npart;
                    etao[Npart]=1.0;
                    used[Npart]=0;
                    qindex[Nquark]=Npart;
                    Npart++;
                    Nquark++;
                }
            }
        }
        // ************* find the quark-anti-quark pairs with smallest distance *************
       //Nquark>Naquark>0
        if (Nquark>Naquark) {
            // if the only one (anti-)quark
            if(Naquark==1){
                int used_q[1000]={0};
                for (int pp=0;pp<Nquark;pp++) {
                    distance[pp]=(phio[aqindex[0]]-phio[qindex[pp]])*(phio[aqindex[0]]-phio[qindex[pp]])+(etao[aqindex[0]]-etao[qindex[pp]])*(etao[aqindex[0]]-etao[qindex[pp]]);
                    used_q[pp]=used[qindex[pp]];
                }
                int minindex=searchmin(distance,qindex,used_q,Nquark);
                apair[Npair]=aqindex[0]; used[aqindex[0]]=1; pair[Npair]=minindex; used[minindex]=1; Npair++;// first pair
                for(int pp=0;pp<Nquark;pp++){
                    int iindex=qindex[pp];
                    if(used[iindex]==1)continue;
                    idpo[Npart]=-3;pxpo[Npart]=0.0;pypo[Npart]=0.0;pzpo[Npart]=100000.0;epo[Npart]=sqrt(pzpo[Npart]*pzpo[Npart]+pypo[Npart]*pypo[Npart]+pxpo[Npart]*pxpo[Npart]+0.5*0.5);
                    apair[Npair]=Npart; used[Npart]=1; pair[Npair]=iindex; used[iindex]=1; Npair++;//  pair
                    Npart++;
                }
            }else{
                //sort the input partons according to pT
                int temindex[1000]={0};
                for(int yy=0;yy<Naquark;yy++){
                     temindex[yy]=aqindex[yy];
                 }
                for (int i = 0; i < Naquark-1; i++){
                    for (int j = 0; j < Naquark - 1 - i; j++){
                        if (ptpo[temindex[j]] > ptpo[temindex[j + 1]]) {
                            tempp = ptpo[temindex[j]];
                            ptpo[temindex[j]] = ptpo[temindex[j + 1]];
                            ptpo[temindex[j+1]] = tempp;
                            tt=aqindex[j];
                            aqindex[j] = aqindex[j+1];
                            aqindex[j+1] = tt; 
                        }
                    }
                }
                // Assign the pairs
                for (int ii=0; ii<Naquark; ii++) {
                    int used_q[1000]={0};
                    for (int pp=0;pp<Nquark;pp++) {
                        distance[pp]=(phio[aqindex[ii]]-phio[qindex[pp]])*(phio[aqindex[ii]]-phio[qindex[pp]])+(etao[aqindex[ii]]-etao[qindex[pp]])*(etao[aqindex[ii]]-etao[qindex[pp]]);
                        used_q[pp]=used[qindex[pp]];
                    }
                    int minindex=searchmin(distance,qindex,used_q,Nquark);
                    apair[Npair]=aqindex[ii]; used[aqindex[ii]]=1; pair[Npair]=minindex; used[minindex]=1; Npair++;//  pair
                }
            }
            // the remnant quarks
            for(int yy=0;yy<Nquark;yy++){
                if(used[qindex[yy]]==1)continue;
                idpo[Npart]=-3;pxpo[Npart]=0.0;pypo[Npart]=0.0;pzpo[Npart]=100000.0;
                epo[Npart]=sqrt(pzpo[Npart]*pzpo[Npart]+pypo[Npart]*pypo[Npart]+pxpo[Npart]*pxpo[Npart]+0.5*0.5);
                xxpo[Npart]=0.0;yypo[Npart]=0.0;zzpo[Npart]=0.0;ttpo[Npart]=0.0;phio[Npart]=atan2(pypo[Npart],pxpo[Npart]);index[Npart]=Npart;
                apair[Npair]=Npart; used[Npart]=1; pair[Npair]=qindex[yy]; used[qindex[yy]]=1; Npair++;//  pair
                Npart++;
            }
        }
    
        //Naquark>Nquark>0
        if(Naquark>Nquark){
            // if the only one (anti-)quark
            if(Nquark==1){
                int used_q[1000]={0};
                for(int pp=0;pp<Naquark;pp++){
                    distance[pp]=(phio[qindex[0]]-phio[aqindex[pp]])*(phio[qindex[0]]-phio[aqindex[pp]])+(etao[qindex[0]]-etao[aqindex[pp]])*(etao[qindex[0]]-etao[aqindex[pp]]);
                    used_q[pp]=used[aqindex[pp]];
                }
                int minindex=searchmin(distance,aqindex,used_q,Naquark);
                pair[Npair]=qindex[0]; used[qindex[0]]=1; apair[Npair]=minindex; used[minindex]=1; Npair++;// first pair
                for(int pp=0;pp<Naquark;pp++){
                    int iindex=aqindex[pp];
                    if(used[iindex]==1)continue;
                    idpo[Npart]=3;pxpo[Npart]=0.0;pypo[Npart]=0.0;pzpo[Npart]=100000.0;epo[Npart]=sqrt(pzpo[Npart]*pzpo[Npart]+pypo[Npart]*pypo[Npart]+pxpo[Npart]*pxpo[Npart]+0.5*0.5);
                    pair[Npair]=Npart; used[Npart]=1; apair[Npair]=iindex; used[iindex]=1; Npair++;//  pair
                    Npart++;
                }
            }else{
                //sort the input partons according to pT// track the index
                int temindex[1000]={0};
                for(int yy=0;yy<Nquark;yy++){
                     temindex[yy]=qindex[yy];
                 }
                for (int i = 0; i < Nquark-1; i++){
                    for (int j = 0; j < Nquark - 1 - i; j++){
                        if (ptpo[temindex[j]] > ptpo[temindex[j + 1]]) {
                            tempp = ptpo[temindex[j]];
                            ptpo[temindex[j]] = ptpo[temindex[j + 1]];
                            ptpo[temindex[j+1]] = tempp;
                            tt=qindex[j];
                            qindex[j] = qindex[j+1];
                            qindex[j+1] = tt; 
                        }
                    }
                }
                // Assign the pairs
                for(int ii=0;ii<Nquark;ii++){
                    int used_q[1000]={0};
                    for(int pp=0;pp<Naquark;pp++){
                        distance[pp]=(phio[qindex[ii]]-phio[aqindex[pp]])*(phio[qindex[ii]]-phio[aqindex[pp]])+(etao[qindex[ii]]-etao[aqindex[pp]])*(etao[qindex[ii]]-etao[aqindex[pp]]);
                        used_q[pp]=used[aqindex[pp]];
                    }
                    int minindex=searchmin(distance,aqindex,used_q,Nquark);
                    pair[Npair]=qindex[ii]; used[qindex[ii]]=1; apair[Npair]=minindex; used[minindex]=1; Npair++;//  pair
                }
            }
            // the remnant ani-quarks
            for(int yy=0;yy<Naquark;yy++){
                if(used[aqindex[yy]]==1)continue;
                idpo[Npart]=3;pxpo[Npart]=0.0;pypo[Npart]=0.0;pzpo[Npart]=100000.0;
                epo[Npart]=sqrt(pzpo[Npart]*pzpo[Npart]+pypo[Npart]*pypo[Npart]+pxpo[Npart]*pxpo[Npart]+0.5*0.5);
                xxpo[Npart]=0.0;yypo[Npart]=0.0;zzpo[Npart]=0.0;ttpo[Npart]=0.0;phio[Npart]=atan2(pypo[Npart],pxpo[Npart]);index[Npart]=Npart;
                pair[Npair]=Npart; used[Npart]=1; apair[Npair]=aqindex[yy]; used[aqindex[yy]]=1; Npair++;//  pair
                Npart++;
            }
        }
    
        //Naquark=Nquark>0
        if(Naquark==Nquark){
            if (Nquark>1) {
                //sort the input partons according to Delta R
                for(int www=0;www<Nquark;www++){
                    double disquark[200][200]={0.0};
                    // calculate the distance between quark pairs
                    for(int ii=0;ii<Nquark;ii++){
                        for(int pp=0;pp<Naquark;pp++){
                            if((used[qindex[ii]]==0)&&(used[aqindex[pp]]==0)){
                                double sphi=phio[qindex[ii]]-phio[aqindex[pp]];
                                if(abs(sphi)>PI){sphi=2*PI-abs(sphi);}
                                disquark[ii][pp]=sphi*sphi+(etao[qindex[ii]]-etao[aqindex[pp]])*(etao[qindex[ii]]-etao[aqindex[pp]]);
                            }else{
                                disquark[ii][pp]=1000000000000000.0;
                            }
                        }
                    }
                    // get the index of the minum value of distance
                    int qnm=0,aqnm=0;
                    double mm=disquark[0][0];
                    for (int i = 0; i < Nquark; i++){
                        for(int j=0;j<Naquark;j++){
                            if (mm > disquark[i][j]){
                                mm = disquark[i][j];
                                qnm = i;
                                aqnm=j;
                            }
                        }
                    }
                    pair[Npair]=qindex[qnm]; used[qindex[qnm]]=1; apair[Npair]=aqindex[aqnm]; used[aqindex[aqnm]]=1; Npair++;//  pair
                }
            }
            if(Naquark==1){
                pair[0]=qindex[0];used[qindex[0]]=1; apair[0]=aqindex[0]; used[aqindex[0]]=1; Npair++;//  pair
            }
        }
        // *** connect gluon the quark anti-quark pairs
        for(int zz=0;zz<Npair;zz++){
            Snum[zz]=0;
        }
        if(Ngluon>0){
            for(int pp=0;pp<Ngluon;pp++){
                double disg[1000]={0.0},disq[1000]={0.0};
                for(int nn=0;nn<Npair;nn++){
                    double sphi=phio[gindex[pp]]-phio[pair[nn]];
                    if(abs(sphi)>PI){sphi=2*PI-abs(sphi);}
                    double d1=sphi*sphi+(etao[gindex[pp]]-etao[pair[nn]])*(etao[gindex[pp]]-etao[pair[nn]]);
                    double sphi2=phio[gindex[pp]]-phio[apair[nn]];
                    if(abs(sphi2)>PI){sphi2=2*PI-abs(sphi2);}
                    double d2=sphi2*sphi2+(etao[gindex[pp]]-etao[apair[nn]])*(etao[gindex[pp]]-etao[apair[nn]]);
                    disg[nn]=sqrt(d1)+sqrt(d2); disq[nn]=d1;
                }
                int minindex=searchmin2(disg,Npair);//minindex belongs to Npair
                    strings[minindex][Snum[minindex]]=gindex[pp];
                    dsting[minindex][Snum[minindex]]=disq[minindex];
                    Snum[minindex]++;
            }
        }
        // sort the gluons // track the index
        for(int kkk=0;kkk<Npair;kkk++){
            if(Snum[kkk]>=2){
                int temindex[1000]={0};
                     for(int yy=0;yy<Snum[kkk];yy++){
                         temindex[yy]=strings[kkk][yy];
                     }
                for (int i = 0; i < Snum[kkk]-1; i++){
                    for (int j = 0; j < Snum[kkk] - 1 - i; j++){
                        if (dsting[kkk][j] > dsting[kkk][j+1]) {
                            tempp = dsting[kkk][j];
                            dsting[kkk][j] = dsting[kkk][j+1];
                            dsting[kkk][j+1] = tempp;
                            tt=strings[kkk][j];
                            strings[kkk][j] = strings[kkk][j+1];
                            strings[kkk][j+1] = tt; 
                        }
                    }
                }
            }
        }
        //***** connect the partons by strings *****
        int Cocon=102; 
        //cout<<"=== The string connections ==="<<endl;
        for(int pp=0;pp<Npair;pp++){
            int Aocon=Cocon+Snum[pp];
            int Cocon2=Cocon;
            //cout<<"====   "<< Snum[pp] <<"   ======"<<endl;
            nncol[pair[pp]]=Cocon;  aacol[pair[pp]]=0;// quark
            //cout<< nncol[pair[pp]]<<" "<<aacol[pair[pp]]<<" "<<idpo[pair[pp]]<<endl;
            for(int yy=0;yy<Snum[pp];yy++){
                int I1=strings[pp][yy]; 
                aacol[I1]=Cocon; Cocon++;
                nncol[I1]=Cocon; //cout<< nncol[I1]<<" "<<idpo[I1]<<" ";cout<< aacol[I1]<<" "<<idpo[I1]<<endl;
            }
            nncol[apair[pp]]=0; aacol[apair[pp]]=Cocon;Cocon++;// anti-quark
            //cout<< nncol[apair[pp]]<<" "<<aacol[apair[pp]]<<" "<<idpo[apair[pp]]<<endl;
        }
        //cout<<"================================"<<endl;
       }
// ****************** append the partons into pythia event *********************//
        double m_str=0.0, x_str=0.0,y_str=0.0,z_str=0.0,t_str=0.0;
        double x_hadron,y_hadron,z_hadron,t_hadron,hmt;
        pythia.event.clear();
        double maxQ0 = 2.0;//maxQ0>=QQ0; wenbin 
        double minQ0 = 0.4;//minum Q0;
        pythia.event.reset();
        for (int tt=0;tt<Npart;tt++){
            if(idpo[tt]==21){mass=0.0;}
            else{
                //if(abs(idp[tt])<=2)mass=0.330;
                //if(abs(idp[tt])==3)mass=0.50;
                mass=sqrt(abs(epo[tt]*epo[tt]-pzpo[tt]*pzpo[tt]-pypo[tt]*pypo[tt]-pxpo[tt]*pxpo[tt]));
            }
            pythia.event.append(idpo[tt],62,nncol[tt],aacol[tt],pxpo[tt],pypo[tt],pzpo[tt],epo[tt],mass);
            if(maxQ0<Qscale[tt])Qscale[tt]=maxQ0;
            //if(minQ0>Qscale[tt])Qscale[tt]=minQ0;
            //pythia.event[tt].scale(Qscale[tt]);//QQ0 the initial scale of input partons 
            // get the center of mass of the strings and corresponding posistion 
            m_str=m_str+mass;
            x_str=x_str+mass*xxpo[tt];
            y_str=y_str+mass*yypo[tt];
            z_str=z_str+mass*zzpo[tt];
            t_str=t_str+mass*ttpo[tt];
        }
        if(m_str==0)m_str=0.10;
        
        if (!DO_Colorless_frag) {
            //first, find unpaired color and anticolor tags. 
            std::vector<int> cols;
            std::vector<int> acols;
            for (unsigned int ipart = 0; ipart < Npart; ++ipart) {
                if (idpo[ipart] == 22) continue;
                if (nncol[ipart] != 0) cols.push_back(nncol[ipart]);
                if (aacol[ipart] != 0) acols.push_back(aacol[ipart]);
            }
            //the outcomes are: 1-unpaired color tag, 2-unpaired anticolor tag, 3-both an unpaired color & anticolor tag, 4-no unpaired tags
            //1-add an antiquark, 2-add a quark, 3-add a gluon, 4-add nothing (possibly photon only event)
            int icol = 0;
            while (icol < cols.size()) {
                bool foundpair = false;
                for (int iacol = 0; iacol < acols.size(); ++iacol) {
                    if (cols[icol] == acols[iacol]) {
                        cols.erase(cols.begin() + icol);
                        acols.erase(acols.begin() + iacol);
                        foundpair = true;
                        continue;
                    }
                }
                if (!foundpair) ++icol;
            }
        
            double sign_added = -1.; double p_fake = 200.;
            while (cols.size() >0 || acols.size() > 0) {
                int pid = 0;
                int color = 0;
                int anti_color = 0;
                if ((cols.size() > 0) && (acols.size() > 0)) {
                    pid = 21;
                    color = cols[0];
                    anti_color = acols[0];
                    cols.erase(cols.begin() );
                    acols.erase(acols.begin());
                } else if ((cols.size() > 0) && (acols.size() == 0)) {
                    pid = -1;
                    color = cols[0];
                    anti_color = 0;
                    cols.erase(cols.begin() );
                } else if ((cols.size() == 0) && (acols.size() > 0)) {
                    pid = 1;
                    color = 0;
                    anti_color = acols[0];
                    acols.erase(acols.begin() );
                }
                if (pid != 0) {
                    p_fake = sign_added * 1. * p_fake;
                    pythia.event.append(pid, 62, anti_color, color, 0.2, 0.2, p_fake,
                                        sqrt(p_fake * p_fake + 0.08));
                    sign_added = sign_added * -1.;
                }
            }
        }
        
//****** fragment the remnant partons ***********
        //pythia.forceTimeShower(1,Npart,maxQ0);//Continue the FSR to the defaulted scale 
        pythia.forceHadronLevel();
        //pgd particle
        char infilepdg[128];
	int pdgid[600]={0};
        sprintf(infilepdg,"chosen_particles.dat");
        ifstream inhypdg(infilepdg);
	int lenght=251;
	for (int ll=0;ll<lenght;ll++) {
		inhypdg>>pdgid[ll];
	}
	inhypdg.close();

        int simble = 0;
        for(int i=0; i<pythia.event.size();i++) {
            if (pythia.event[i].isFinal() ) {
                c_id = pythia.event[i].id();
                bool sss=findarray(pdgid,lenght, c_id);
                if (sss) {
                    simble=simble+1;
                }
            }
        }
        //if(simble==0){output2 << iEvent+1<<" "<<simble << endl;}
        if (simble > 0) {
           /*
            output2 << "         " << iEvent <<"          " << simble << "         0         0" << endl;
            for(int i=0; i<pythia.event.size();i++) {
                if (pythia.event[i].isFinal() ){
                    c_id = pythia.event[i].id();
                    bool sss = findarray(pdgid,lenght, c_id);
                    if (sss) {
                        cbar_meson_px = pythia.event[i].px();
                        cbar_meson_py = pythia.event[i].py();
                        cbar_meson_pz = pythia.event[i].pz();
                        cbar_meson_energy = pythia.event[i].e();
                        // get the posistion of the final hadrons
                        hmt=sqrt(pythia.event[i].m()*pythia.event[i].m()+cbar_meson_px*cbar_meson_px+cbar_meson_py*cbar_meson_py);
                        x_hadron=x_str/m_str+hbarc*cbar_meson_px/hmt;
                        y_hadron=y_str/m_str+hbarc*cbar_meson_py/hmt;
                        z_hadron=z_str/m_str+hbarc*cbar_meson_pz/hmt;
                        t_hadron=t_str/m_str+hbarc*cbar_meson_energy/hmt;
			output2 << "         " << i <<"        "<<c_id<<"    "<<cbar_meson_px<<"    "<<cbar_meson_py<<"    "<<cbar_meson_pz
			        << "    " << cbar_meson_energy <<"    "<< pythia.event[i].m() 
			        << "    " << x_hadron <<"    "<< y_hadron <<"    "<< z_hadron <<"    "<< t_hadron <<endl;
                    }
                }
            }
            */
            // *************** fastjet ******************** //
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
    
            // Loop over events
            vector<fastjet::PseudoJet> input_particles;
            // First use the fastjet to pre-select the event with jet pT > 500 GeV at parton level
            for (int i = 0; i < pythia.event.size(); i++) {
                fastjet::PseudoJet particlefastjet = PseudoJet(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(),
                                                               pythia.event[i].e());
                particlefastjet.set_user_index(pythia.event[i].id());
                input_particles.push_back(particlefastjet);
            }
            input_particles = particle_selector(input_particles);
        
            // Then do the jet finding
            fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R_jet);
            // select jet
            //input_particles = particle_selector(input_particles);
            fastjet::ClusterSequence clust_seq(input_particles, jet_def);
            // get the resulting jets ordered in pt
            vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt( jet_selector(clust_seq.inclusive_jets()) );
        
            if (inclusive_jets.size() == 0 || inclusive_jets[0].pt() < 550.) continue;
            
        
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
                        if (constituents[jj].pt() > 0.3 && abs(constituents[jj].eta()) < 2.4 && 
                            (abs(constituents[jj].user_index()) ==211 || abs(constituents[jj].user_index()) ==321 || 
                             abs(constituents[jj].user_index()) ==2212) ) Nparticle ++;
                    }
                    Nchevent.push_back(Nparticle);
                }
            }
        }
        /*
        if (Nchevent.size() > 0) {
           output2 << "         " << iEvent <<"          " << findMax(Nchevent) << endl;
        } else {
            output2 << "         " << iEvent << "      " << 0 << endl;
        }
        */
        if (findMax(Nchevent) > Nchcut) { //pre-select high Nch events
            output << "# event id " << iEvent  <<  ", Number_of_parton  "  << std::endl;
            output << Npart << std::endl;
            for (int ll=0; ll < Npart; ll++) {
                output << idpo[ll] << "  " 
                       << pxpo[ll] << "  " << pypo[ll] << "  " << pzpo[ll] << "  " << epo[ll] <<"  "
                       << xxpo[ll] << "  "<< yypo[ll] << "  " << zzpo[ll] << "  " << ttpo[ll] << "  "
                       << nncol[ll] << "  " << aacol[ll]
                       << std::endl;
            }
        }
    }
    output2.close();
    output.close();
    fclose(infile1);
    // End Pythia8
}

int searchmin(double*p,int*q,int*s,int len)
{
    double m = 100000000.0;
    int k;
    for (int i = 0; i < len; ++i)
    {
        if ((m > p[i])&&(s[i]==0))
        {
            m = p[i];
            k = i;
        }
    }
    return q[k];
} 


int searchmax(double*p,int*q,int len)
{
    double m = 0.0;
    int k;
    for (int i = 0; i < len; ++i)
    {
        if (m < p[i])
        {
            m = p[i];
            k = i;
        }
    }
    return q[k];
} 

int findMax(std::vector<int>& vec) {
    if (vec.empty()) {
        // Handle the case where the vector is empty.
        std::cerr << "Error: Empty vector\n";
        return -1;  // You can choose a different value or use an exception.
    }

    int max = vec[0];  // Assume the first element is the maximum.

    for (int i = 1; i < vec.size(); ++i) {
        if (vec[i] > max) {
            max = vec[i];  // Update the maximum value.
        }
    }

    return max;
}


// check particle ID
bool findarray(int*p, int len,int val)
{
	int i;
	bool ret = false;
	for (i = 0; i!= len; i++)
	{
		if (p[i] == val)
		{ret = true; break;}
	}
	return ret;
} 

    
int searchmin2(double*p,int len)
{
    double m = 100000000.0;
    int k;
    for (int i = 0; i < len; ++i)
    {
        if (m > p[i])
        {
            m = p[i];
            k = i;
        }
    }
    return k;
} 
    
