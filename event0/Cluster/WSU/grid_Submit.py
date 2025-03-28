#!/usr/bin/env python
#author: Wenbin Zhao
#email: wenbinzhao237@gmail.com

from subprocess import call
import sys
import random
import time
def submit(fold_id_start=0, fold_id_end = 1, nevent = 10, random_number = 0):
    jobs = '''#!/bin/bash
#SBATCH --job-name SoftZPC{fold_id_start}
#SBATCH -t 16:59:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3gb
#SBATCH --constraint="intel"
#SBATCH --mail-type=all
#SBATCH --mail-user=XXXX@wayne.edu
#SBATCH -o log/%j.out
#SBATCH -e log/%j.err
#SBATCH -q primary



hostname
date
 
#module load cmake/3.20.0
#module load python/3.9.5
#module load grid-default
#module load eigen/3.3.9
#module load autotools
#module load gnu7/7.3.0
#module load gsl/2.4
#module load boost/1.74.0
#module load prun/1.2
#module load openmpi3/3.1.6
#module load hdf5/1.10.7

# First link the file
#mkdir event0
#cd event0 

#mkdir fastjet_hadron
#mkdir hadronization_urqmd
#mkdir pythia_parton
#mkdir ZPC

#cd fastjet_hadron
#ln -s ../../Jet_parton_cascade/fastjet_hadron/* ./

#cd ../hadronization_urqmd
#mkdir fragmentation
#mkdir urqmd_code
#cd fragmentation 
#ln -s ../../../Jet_parton_cascade/hadronization_urqmd/fragmentation/* ./ 
#cd ../urqmd_code
#cp -r ../../../Jet_parton_cascade/hadronization_urqmd/urqmd_code/run_afterburner.sh ./
#mkdir osc2u
#mkdir urqmd
#cd osc2u
#ln -s ../../../../Jet_parton_cascade/hadronization_urqmd/urqmd_code/osc2u/* ./
#cd ../urqmd
#ln -s ../../../../Jet_parton_cascade/hadronization_urqmd/urqmd_code/urqmd/* ./
#cd ../../../

#cd pythia_parton
#ln -s ../../Jet_parton_cascade/pythia_parton/* ./

#cd ../ZPC
#ln -s ../../Jet_parton_cascade/ZPC/* ./
#mkdir ana
#cd ../

# Then run the framework
for (( ii={fold_id_start}; ii<{fold_id_end}; ii++ ))
do

cp -r event0 PlaygroundZPC/job-$ii
cd PlaygroundZPC/job-$ii

# Generate the pythia parton
#sleep 1s
cd pythia_parton
./mymain06 {nevent} $(({random_number} + $ii * 12345))
cd ../

# ZPC for parton cascade
cd  ZPC
ln -s ../pythia_parton/parton_info.dat ./
./exec
rm -r ana/parton-collisionsHistory.dat
rm -r ana/zpc.res
cd ../
rm -rf pythia_parton/parton_info.dat

# fragmentation and urqmd
cd hadronization_urqmd
cd fragmentation
ln -s ../../ZPC/ana/zpc.dat ./
./main_string_fragmentation {nevent}
rm -r ../../ZPC/ana/*

cd ../urqmd_code
    # script to run urqmd
    cd osc2u
    ln -s ../../fragmentation/hadrons_frag1.dat ./
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
ln -s ../hadronization_urqmd/urqmd_code/urqmd/particle_list.dat ./
./fastjet_hadron {nevent}
rm -r ../hadronization_urqmd/urqmd_code/urqmd/particle_list.dat
cd ../

# Save the final results into folder
mkdir results
mv fastjet_hadron/final_state_hard_hadrons.bin ./results/$ii.bin
rm -r fastjet_hadron
rm -r hadronization_urqmd 
rm -r pythia_parton 
rm -r ZPC 
cd ../../
done
'''.format(fold_id_start=fold_id_start, fold_id_end = fold_id_end, nevent = nevent, random_number = random_number)

    job_name = "NSC3_%s.sh"%(fold_id_start)
    with open(job_name, 'w') as fout:
        fout.write(jobs)
    call(['sbatch', job_name])
    call(['mv', job_name, 'jobs/'])
# comparing the grid size dependence of the hypersf cube
# etaos_ymin = 0.08 for both ampt runs, not 0.2

if __name__=='__main__':
    import sys
    fold_id_start = int(int(sys.argv[1]))
    nfold = int(int(sys.argv[2]))
    nevent = int(int(sys.argv[3]))
    # Get the current system time as a seed
    seed = int(time.time())

    # Seed the random number generator
    random.seed(seed)

    # Generate a random number between 0 and 10,000,000
    random_number = random.randint(0, 10**8)
    #for n in range(0,nods):
    #mmid = fold_id_start * tot_ev + n 
    submit(fold_id_start * nfold, nfold + fold_id_start * nfold, nevent, random_number + fold_id_start)
    #start_num += jobs_per_cpu 

