#!/bin/bash

rm -rf event0
# First link the file
mkdir event0
cd event0 

mkdir fastjet_hadron
mkdir hadronization_urqmd
mkdir pythia_parton
mkdir ZPC

cd fastjet_hadron
ln -s /wsu/home/he/he92/he9215/Berkeley_work/V2_in_jet/Jet_parton_cascade/fastjet_hadron/* ./

cd ../hadronization_urqmd
mkdir fragmentation
mkdir urqmd_code
cd fragmentation 
ln -s /wsu/home/he/he92/he9215/Berkeley_work/V2_in_jet/Jet_parton_cascade/hadronization_urqmd/fragmentation/* ./ 
cd ../urqmd_code
cp -r /wsu/home/he/he92/he9215/Berkeley_work/V2_in_jet/Jet_parton_cascade/hadronization_urqmd/urqmd_code/run_afterburner.sh ./
mkdir osc2u
mkdir urqmd
cd osc2u
ln -s /wsu/home/he/he92/he9215/Berkeley_work/V2_in_jet/Jet_parton_cascade/hadronization_urqmd/urqmd_code/osc2u/* ./
cd ../urqmd
ln -s /wsu/home/he/he92/he9215/Berkeley_work/V2_in_jet/Jet_parton_cascade/hadronization_urqmd/urqmd_code/urqmd/* ./
cd ../../../

cd pythia_parton
ln -s /wsu/home/he/he92/he9215/Berkeley_work/V2_in_jet/Jet_parton_cascade/pythia_parton/* ./

cd ../ZPC
ln -s /wsu/home/he/he92/he9215/Berkeley_work/V2_in_jet/Jet_parton_cascade/ZPC/* ./
mkdir ana
cd ../


