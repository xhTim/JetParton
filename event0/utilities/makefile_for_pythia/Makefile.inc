# PYTHIA configuration file.
# Generated on 2023年 10月 19日 星期四 00:30:02 PDT with the user supplied options:
# --prefix=/home/wenbin/Downloads/Wenbin_working/Work/Berkeley_work/v2_in_jet/parton_shower_in_pythia8/pythia8_install
# --with-lhapdf6=/home/wenbin/Downloads/Wenbin_working/Work/Berkeley_work/v2_in_jet/parton_shower_in_pythia8/RESULTS/T2/pythia_hadron/LHAPDF_install

# Install directory prefixes.
PREFIX_BIN=/home/wenbin/Downloads/Wenbin_working/Work/Berkeley_work/v2_in_jet/parton_shower_in_pythia8/pythia8_install/bin
PREFIX_INCLUDE=/home/wenbin/Downloads/Wenbin_working/Work/Berkeley_work/v2_in_jet/parton_shower_in_pythia8/pythia8_install/include
PREFIX_LIB=/home/wenbin/Downloads/Wenbin_working/Work/Berkeley_work/v2_in_jet/parton_shower_in_pythia8/pythia8_install/lib
PREFIX_SHARE=/home/wenbin/Downloads/Wenbin_working/Work/Berkeley_work/v2_in_jet/parton_shower_in_pythia8/pythia8_install/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
CXX=g++
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so
OBJ_COMMON=

EVTGEN_USE=false
EVTGEN_CONFIG=
EVTGEN_BIN=
EVTGEN_INCLUDE=
EVTGEN_LIB=

FASTJET3_USE=false
FASTJET3_CONFIG=
FASTJET3_BIN=
FASTJET3_INCLUDE=
FASTJET3_LIB=

HEPMC2_USE=false
HEPMC2_CONFIG=
HEPMC2_BIN=
HEPMC2_INCLUDE=
HEPMC2_LIB=

HEPMC3_USE=false
HEPMC3_CONFIG=
HEPMC3_BIN=
HEPMC3_INCLUDE=
HEPMC3_LIB=

LHAPDF5_USE=false
LHAPDF5_CONFIG=
LHAPDF5_BIN=
LHAPDF5_INCLUDE=
LHAPDF5_LIB=

LHAPDF6_USE=true
LHAPDF6_CONFIG=lhapdf-config
LHAPDF6_BIN=/home/wenbin/Downloads/Wenbin_working/Work/Berkeley_work/v2_in_jet/parton_shower_in_pythia8/RESULTS/T2/pythia_hadron/LHAPDF_install/bin/
LHAPDF6_INCLUDE=-I/home/wenbin/Downloads/Wenbin_working/Work/Berkeley_work/v2_in_jet/parton_shower_in_pythia8/RESULTS/T2/pythia_hadron/LHAPDF_install/include
LHAPDF6_LIB=-L/home/wenbin/Downloads/Wenbin_working/Work/Berkeley_work/v2_in_jet/parton_shower_in_pythia8/RESULTS/T2/pythia_hadron/LHAPDF_install/lib -Wl,-rpath,/home/wenbin/Downloads/Wenbin_working/Work/Berkeley_work/v2_in_jet/parton_shower_in_pythia8/RESULTS/T2/pythia_hadron/LHAPDF_install/lib -lLHAPDF

POWHEG_USE=false
POWHEG_CONFIG=
POWHEG_BIN=
POWHEG_INCLUDE=
POWHEG_LIB=

RIVET_USE=false
RIVET_CONFIG=
RIVET_BIN=
RIVET_INCLUDE=
RIVET_LIB=

ROOT_USE=false
ROOT_CONFIG=
ROOT_BIN=
ROOT_INCLUDE=
ROOT_LIB=

GZIP_USE=false
GZIP_CONFIG=
GZIP_BIN=
GZIP_INCLUDE=
GZIP_LIB=

PYTHON_USE=false
PYTHON_CONFIG=
PYTHON_BIN=
PYTHON_INCLUDE=
PYTHON_LIB=

MG5MES_USE=false
MG5MES_CONFIG=
MG5MES_BIN=
MG5MES_INCLUDE=
MG5MES_LIB=

OPENMP_USE=false
OPENMP_CONFIG=
OPENMP_BIN=
OPENMP_INCLUDE=
OPENMP_LIB=

MPICH_USE=false
MPICH_CONFIG=
MPICH_BIN=
MPICH_INCLUDE=
MPICH_LIB=

HDF5_USE=false
HDF5_CONFIG=
HDF5_BIN=
HDF5_INCLUDE=
HDF5_LIB=

HIGHFIVE_USE=false
HIGHFIVE_CONFIG=
HIGHFIVE_BIN=
HIGHFIVE_INCLUDE=
HIGHFIVE_LIB=
