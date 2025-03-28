#!/usr/bin/env python
#author: Wenbin Zhao
#email: wenbinzhao237@gmail.com

from subprocess import call
import sys
import random
import time
def submit(fold_id=0,nevent = 10, random_number = 0):
    jobs = '''#!/bin/bash
#SBATCH --job-name SoftHard{fold_id}
#SBATCH -t 23:59:00
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

cp -r event0/pythia_hadron Playground/job-{fold_id}
cd Playground/job-{fold_id}

# Generate the pythia parton
#cd pythia_parton
./mymain07 {nevent} {random_number}

mkdir results
mv fastjet_hadron/final_state_hard_hadrons.bin ./results/{fold_id}.bin
cd ../
'''.format(fold_id=fold_id, nevent = nevent, random_number = random_number)

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
    # Get the current system time as a seed
    seed = int(time.time())

    # Seed the random number generator
    random.seed(seed)

    # Generate a random number between 0 and 10,000,000
    random_number = random.randint(0, 10**4)
    
    submit(fold_id,nevent, random_number*fold_id+fold_id)
    print(random_number+fold_id)
    #start_num += jobs_per_cpu 
