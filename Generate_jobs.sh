#!/bin/bash

#sleep 0.1h
#rm -rf /pscratch/sd/w/wenbinz/V2inJet/event0
#rm -rf /pscratch/sd/w/wenbinz/V2inJet/grid_Submit.py
#cp -r event0 /pscratch/sd/w/wenbinz/V2inJet
#cp -r grid_Submit.py /pscratch/sd/w/wenbinz/V2inJet
#cd /pscratch/sd/w/wenbinz/V2inJet
#mkdir log
#mkdir jobs
mkdir /eos/cms/store/group/phys_heavyions/huangxi/Wiese/Playground
((Randum_number=$RANDOM))

((Randum_number=$RANDOM))
for (( ii=0; ii<1; ii++ ))
do
    ./grid_Submit.py $ii 60 100000 # !!! 60 is the number of tasks, 100000 is the number of events in each output. Total number of events is 60^2*100000 !!!
    #./testgrid_Submittest.py $ii 20 20000
    #sleep 10s
    #./grid_Submit_pythia_only.py $ii 600000

    #sbatch testCompile_inte.sh $ii $Randum_number
done



