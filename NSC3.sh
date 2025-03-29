#!/bin/bash

cd /eos/cms/store/group/phys_heavyions/huangxi/Wiese
cp -r event0 Playground/job-$1
cd Playground/job-$1

# Generate the pythia parton
#sleep 1s
cd pythia_parton
./mymain06 {nevent} 0 #$(({random_number} + $1 * 12345 + $1 * 11))
cd ../

# ZPC for parton cascade
cd  ZPC
mkdir -p ana
ln -sf ../pythia_parton/parton_info.dat ./
./exec
rm -r ana/parton-collisionsHistory.dat
rm -r ana/zpc.res
cd ../
rm -rf pythia_parton/parton_info.dat

# fragmentation and urqmd
cd hadronization_urqmd
cd fragmentation
ln -sf ../../ZPC/ana/zpc.dat ./
./main_string_fragmentation {nevent}
rm -r ../../ZPC/ana/*

cd ../urqmd_code
    # script to run urqmd
    cd osc2u
    ln -sf ../../fragmentation/hadrons_frag1.dat ./
    ./osc2u.e < hadrons_frag1.dat > run.log
    rm -r ../../fragmentation/hadrons_frag1.dat
    mv fort.14 ../urqmd/OSCAR.input
    cd ../urqmd
    ./runqmd.sh > run.log
    rm -fr OSCAR.input
    rm -rf run.log
    cd ..
cd ../
cd ../
# jet finding of final hadrons
cd fastjet_hadron
ln -sf ../hadronization_urqmd/urqmd_code/urqmd/particle_list.dat ./
./fastjet_hadron_trackTree {nevent} $1
rm -r ../hadronization_urqmd/urqmd_code/urqmd/particle_list.dat
cd ../

# Save the final results into folder
# mkdir -p results
# mv fastjet_hadron/final_state_hard_hadrons.bin ./results/$1.bin
rm -r fastjet_hadron
rm -r hadronization_urqmd
rm -r pythia_parton
rm -r ZPC
cd ../
#rm -rf job-$1
cd ../
