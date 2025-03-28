This is the framework to couple the jet parton shower and parton cascade

# Install LHAPDF
1. **Download the LHAPDF Package**  
   Download the LHAPDF compressed file
2. mkdir /Tool/LHAPDF_Lib
3. cd /path/to/extracted-LHAPDF-directory
4. ./configure --prefix=/Tool/LHAPDF_Lib --disable-python
5. make -j8

# Install PDF Files
1. wget http://lhapdfsets.web.cern.ch/lhapdfsets/current/NNPDF31_nnlo_as_0118.tar.gz
2. tar -xzvf NNPDF31_nnlo_as_0018.tar.gz
3. mv NNPDF31_nnlo_as_0018 /Tool/LHAPDF_Lib/share/LHAPDF/

# Compile pythia8 with LHAPDF
1. download pythia8.310
2. After extracting Pythia, enter the resulting directory.
3. ./configure --prefix=/Tool/pythia8310_install --with-lhapdf6=/Tool/LHAPDF_install
4. make -j8
5. make install

# Modify and compile the code
1. Clone the repository to your lxplus home. Then mv the 'event0' fold to your eos disk (in my case, mv event0 /eos/cms/store/group/phys_heavyions/huangxi/Wiese/).
2. Modify grid_Submit.py line 20 'cd /eos/cms/store/group/phys_heavyions/huangxi/Wiese': cd to your event0 directory.
3. Modify Generate_job.sh line 11 'mkdir /eos/cms/store/group/phys_heavyions/huangxi/Wiese/Playground': create the folder 'Playground' in the same directory as event0.
4. Modify paths to the final output root files in event0/fastjet_hadron/fastjet_hadron_trackTree.cc line 222 'TFile * fout = TFile::Open( Form("/eos/cms/store/group/phys_heavyions/huangxi/PC/pp_parton_cascade_%d.root",jobnumber) ,"recreate");'
5. Modify paths to Pythia and fastjet in Makefile in pythia_parton/, hadronization_urqmd/fragmentation/ and fastjet_hadron/.
6. Go to event0/
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

# Run the code
You can set number of tasks and number of events in each task in Generate_job.sh.
Then submit jobs: ./Generate_jobs.sh
