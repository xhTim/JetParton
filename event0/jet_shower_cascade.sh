#!/bin/bash

#This is the script to run the jet shower and parton cascade framework
# Copy right Wenbin Zhao @ 2023

# First compile all the code
cd pythia_parton
make
cd ../ZPC
make 
cd ../hadronization_urqmd
cd fragmentation
make 
cd ../urqmd_code
FC=gfortran make
cd ../../
cd fastjet_hadron
make 
cd ../

# Then run the framework

# Generate the pythia parton
cd pythia_parton
./mymain06 10000 0
cd ../

# ZPC for parton cascade
cd  ZPC
ln -s ../pythia_parton/parton_info.dat ./
./exec
cd ../

# fragmentation and urqmd
cd hadronization_urqmd
cd fragmentation
ln -s ../../ZPC/ana/zpc.dat ./
./main_string_fragmentation 10000

cd ../urqmd_code
    # script to run urqmd
    cd osc2u
    ln -s ../../fragmentation/hadrons_frag1.dat ./
    ./osc2u.e < hadrons_frag1.dat > run.log
    mv fort.14 ../urqmd/OSCAR.input
    cd ../urqmd
    ./runqmd.sh > run.log
    rm -fr OSCAR.input
    cd ..
cd ../
cd ../
# jet finding of final hadrons
# cd fastjet_hadron
# ln -s ../hadronization_urqmd/urqmd_code/urqmd/particle_list.dat ./
# ./fastjet_hadron 100
# cd ../

# Save the final results into folder
# mkdir results
# cp -r fastjet_hadron/jet_hardons ./results

