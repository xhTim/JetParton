EXE=fastjet_hadron_trackTree
GCC=g++

PYTHIAINC=${shell ~/Tools/pythia8310/bin/pythia8-config --includedir}
PYTHIALIBS=${shell ~/Tools/pythia8310/bin/pythia8-config --libdir}

FASTJETINC=${shell ~/Tools/fastjet-install/bin/fastjet-config --cxxflags}
FASTJETLIBS=${shell ~/Tools/fastjet-install/bin/fastjet-config --libs}

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)

all: $(EXE)

fastjet_hadron_trackTree.o: fastjet_hadron_trackTree.cc
	$(GCC) -I$(PYTHIAINC) $(FASTJETINC) $(ROOTCFLAGS) -c fastjet_hadron_trackTree.cc

$(EXE): fastjet_hadron_trackTree.o
	$(GCC) -O3 -I$(PYTHIAINC) $(FASTJETINC) $(ROOTCFLAGS) -o $(EXE) fastjet_hadron_trackTree.o $(PYTHIALIBS)/libpythia8.a $(ROOTLIBS) -ldl $(FASTJETLIBS) 

.PHONY:clean
clean:
	-rm *.o $(EXE)
