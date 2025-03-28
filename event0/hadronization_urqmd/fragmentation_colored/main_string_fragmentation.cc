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
int searchmin(double*p, int *q,int*s,int len);
int searchmax(double*p, int *q,int len);
int searchmin2(double*p,int len);

char infiles[128];

//===========================================================================================
int main(int argv, char* argc[])
{
    string random_str = string(argc[1]);
    int n_event=atoi(argc[1]);
    string output_filename2;
    output_filename2 = "hadrons_frag1.dat";// output files of final hadrons
    cout << output_filename2 << endl;
    string ramdomseed_str = "Random:seed = "+random_str;
    ofstream output2(output_filename2.c_str()); 
    if (!output2.is_open() ) {
        cout << "cannot open output file:"<< endl
         << output_filename2 << endl;
        return -1;
    }
    string ipput_filename2;
    ipput_filename2 = "parton_info.dat";// output files of final hadrons
    sprintf(infiles, ipput_filename2.c_str()); //input parton file 
    FILE* infile1;
    infile1 = fopen(infiles,"r");

    Pythia pythia;
    pythia.readString("Random:setSeed = on");
    //pythia.readString("Random:seed = 0");  
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
    //pythia.readString("StringFlav:BtoMratio=0.5");
    pythia.readString("StringFlav:probQQtoQ=0.34");
    pythia.readString("HadronLevel:Decay = off");
    pythia.readString("HadronLevel:Hadronize = on");
    pythia.init();

    double hbarc = 0.19732;
    double c_px, c_py ,c_pz, c_e, c_m,c_x,c_y,c_z,c_t;
    int c_id,tt;
    int c_col,c_acol;
    int acol_ip[5000]={0};
    int pp_collision=0;    //used for total cross section
    int  ccbar_num=0;
    double cmeson_px,cmeson_py,cmeson_pz, cmeson_P,tempp;
    double cbar_meson_px,cbar_meson_py,cbar_meson_pz, cbar_meson_P,pt_square,cbar_meson_energy;
    double pt,mass,amid,Qmid;
    int mid,ie,II,status,col,acol,Npart,NN,Ntotal,PAosi,simble,mmaxindex;
    int idp[10000]={0},idpo[10000]={0};
    double pxpo[10000]={0.0},pypo[10000]={0.0},pzpo[10000]={0.0},epo[10000]={0.0},ptpo[10000]={0.0},xxpo[10000]={0.0},yypo[10000]={0.0},zzpo[10000]={0.0},ttpo[10000]={0.0},phio[10000]={0.0},etao[10000]={0.0};//,mass[10000]={0.0};
    double pxp[10000]={0.0},pyp[10000]={0.0},pzp[10000]={0.0},ep[10000]={0.0},ptp[10000]={0.0},xxp[10000]={0.0},yyp[10000]={0.0},zzp[10000]={0.0},ttp[10000]={0.0},phi[10000]={0.0},distance[10000]={0.0},dsting[100][1000]={0.0},Qscale[10000]={0.0};//,mass[10000]={0.0};
    int nncol[10000]={0},aacol[10000]={0},index[10000]={0},qindex[10000]={0},aqindex[10000]={0},used[10000]={0},pair[1000]={0},apair[1000]={0},gindex[1000]={0},strings[100][1000]={0},Snum[1000]={0};//,nncol_mid[10000]={0},aacol_mid[10000]={0};
    // event loop
    for (int iEvent=0; iEvent<n_event; iEvent++) {     
        fscanf(infile1,"%d\n",&Npart);
        if (Npart==0 ) {output2 << iEvent+1<<" "<<0 << endl;continue;} 
        int Nquark=0;
        int Naquark =0;
        int Ngluon=0;
        int Npair=0;
        for (int ll=0;ll<Npart;ll++) {
                fscanf(infile1,"%d %lf %lf %lf %lf %lf %lf %lf %lf %d %d\n",
                       &c_id,&c_px,&c_py,&c_pz,&c_e,&c_x, &c_y, &c_z,&c_t, &c_col, &c_acol);// format of input partons
                Qmid = 0.0; // The scale for the parton shower. 
                idpo[ll]=c_id;
                if(c_id==21){gindex[Ngluon]=ll;Ngluon++;}
                pxpo[ll]=c_px;
                pypo[ll]=c_py;
                pzpo[ll]=c_pz;
                ptpo[ll]=c_px*c_px+c_py*c_py;
                double pmg=sqrt(c_px*c_px+c_py*c_py+c_pz*c_pz);
                epo[ll]=c_e;
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
                nncol[ll]=c_col;
                aacol[ll]=c_acol;
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
        
        //first, find unpaired color and anticolor tags. 
        // From JETSCPAE colored hadronization
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
        
        
//****** fragment the colored partons ***********
        //pythia.forceTimeShower(1,Npart,maxQ0);//Continue the FSR to the defaulted scale 
        pythia.forceHadronLevel();
        int simble = 0;
        for(int i=0; i<pythia.event.size();i++) {
            if (pythia.event[i].isFinal() ) {
                c_id = pythia.event[i].id();
                if( (abs(c_id) !=22)&&(abs(c_id) !=11)) {
                    simble=simble+1;
                }
            }
        }
        if(simble==0){output2 << iEvent+1<<" "<<simble << endl;}
        if(simble>0){
            output2 << iEvent+1<<" "<<simble << endl;
            for(int i=0; i<pythia.event.size();i++)
                {
                if (pythia.event[i].isFinal() ){
                    c_id = pythia.event[i].id();
                    if ( (abs(c_id) !=22)&&(abs(c_id) !=11)) {
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
                        output2 << c_id << " "<<cbar_meson_px<<"  "<<cbar_meson_py<<"  "<<cbar_meson_pz <<" "<<cbar_meson_energy << "  "<< pythia.event[i].m()<<" "<< x_hadron<<" "<<y_hadron<<" "<<z_hadron<<" "<<t_hadron<<" "<<endl;
                    }
                }
            }
        }
        pythia.next();

}
output2.close();


  return 0;
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
    
