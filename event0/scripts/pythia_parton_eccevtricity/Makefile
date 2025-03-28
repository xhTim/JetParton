GCC=g++


PYTHIAINC=${shell /home/wenbin/Downloads/Wenbin_working/Work/Berkeley_work/v2_in_jet/parton_shower_in_pythia8/pythia8310/bin/pythia8-config --includedir}
PYTHIALIBS=${shell /home/wenbin/Downloads/Wenbin_working/Work/Berkeley_work/v2_in_jet/parton_shower_in_pythia8/pythia8310/bin/pythia8-config --libdir}

FASTJETINC=${shell /home/wenbin/Downloads/Wenbin_working/Work/code/fastjet_install/bin/fastjet-config --cxxflags}
FASTJETLIBS=${shell /home/wenbin/Downloads/Wenbin_working/Work/code/fastjet_install/bin/fastjet-config --libs}

EXE=main
all: $(EXE)

main.o: main.cc
	$(GCC) -I$(PYTHIAINC) $(ROOTINC) $(FASTJETINC) -c main.cc

$(EXE): main.o
	$(GCC) -O3 -I$(PYTHIAINC) $(ROOTINC) $(FASTJETINC) -o $(EXE) main.o $(PYTHIALIBS)/libpythia8.a $(ROOTLIBS) $(FASTJETLIBS) 

.PHONY:clean
clean:
	-rm *.o $(EXE)	
