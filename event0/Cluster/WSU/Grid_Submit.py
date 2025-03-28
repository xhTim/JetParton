#!/usr/bin/env python
#author: Wenbin Zhao
#email: wenbinzhao237@gmail.com

from subprocess import call
import sys
def submit(fold_id=0,nevent = 10):
    jobs = '''#!/bin/bash
#SBATCH --job-name SoftHard{fold_id}
#SBATCH -t 23:59:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2gb
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

cp -r event0 Playground/job-{fold_id}
cd Playground/job-{fold_id}

# Generate the pythia parton
cd pythia_parton
./mymain06 {nevent}
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
./main_string_fragmentation {nevent}

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
cd fastjet_hadron
ln -s ../hadronization_urqmd/urqmd_code/urqmd/particle_list.dat ./
./fastjet_hadron {nevent}
cd ../

# Save the final results into folder
mkdir results
mv fastjet_hadron/final_state_hard_hadrons.bin ./results/{fold_id}.bin

'''.format(fold_id=fold_id, nevent = nevent)

    job_name = "NSC3_%s.sh"%(fold_id)
    with open(job_name, 'w') as fout:
        fout.write(jobs)
    call(['sbatch', job_name])
    call(['mv', job_name, 'jobs/'])
# comparing the grid size dependence of the hypersf cube
# etaos_ymin = 0.08 for both ampt runs, not 0.2

if __name__=='__main__':
    import sys
    fold_id = int(int(sys.argv[1]))
    nevent = int(int(sys.argv[2]))
    #for n in range(0,nods):
    #mmid = fold_id * tot_ev + n 
    submit(fold_id,nevent)
    #start_num += jobs_per_cpu 

