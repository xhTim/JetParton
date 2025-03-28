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
    // Read the binary output
    stringstream strparticle_name;
    strparticle_name << "final_state_hard_hadrons.bin";
    std::ifstream InStream;
    InStream.precision(15);
    string particle_name = strparticle_name.str();
    InStream.open(particle_name.c_str(), std::ios::in | std::ios::binary);
    double mass, pT, phi, eta;
    float temp;
    int iev = 0;
    while( InStream.read(reinterpret_cast<char*>(&total_number_of_particles), sizeof(int))) {
        iev = iev +1;
        for (auto i=0; i<total_number_of_particles; i++) {
            InStream.read(reinterpret_cast<char*>(&pid), sizeof(int));
            InStream.read(reinterpret_cast<char*>(&temp), sizeof(float));
            mass = temp;
            InStream.read(reinterpret_cast<char*>(&temp), sizeof(float));
            pT = temp;
            InStream.read(reinterpret_cast<char*>(&temp), sizeof(float));
            phi = temp;
            InStream.read(reinterpret_cast<char*>(&temp), sizeof(float));
            eta = temp;
//            cout <<total_number_of_particles << " " <<  iev << "  " << pid << "  " << mass << "  " << pT << "  " << phi << "  " << eta << endl;
        }
    }
    InStream.close();


}

