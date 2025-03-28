EXE=main_frag_fastjet
GCC=g++

PYTHIAINC=${shell /home/wenbin/Downloads/Wenbin_working/Work/Berkeley_work/v2_in_jet/parton_shower_in_pythia8/pythia8310/bin/pythia8-config --includedir}
PYTHIALIBS=${shell /home/wenbin/Downloads/Wenbin_working/Work/Berkeley_work/v2_in_jet/parton_shower_in_pythia8/pythia8310/bin/pythia8-config --libdir}

FASTJETINC=${shell /home/wenbin/Downloads/Wenbin_working/Work/code/fastjet_install/bin/fastjet-config --cxxflags}
FASTJETLIBS=${shell /home/wenbin/Downloads/Wenbin_working/Work/code/fastjet_install/bin/fastjet-config --libs}


all: $(EXE)

main_frag_fastjet.o: main_frag_fastjet.cc
	$(GCC) -I$(PYTHIAINC) $(FASTJETINC) -c main_frag_fastjet.cc

$(EXE): main_frag_fastjet.o
	$(GCC) -O3 -I$(PYTHIAINC) $(FASTJETINC) -o $(EXE) main_frag_fastjet.o $(PYTHIALIBS)/libpythia8.a -ldl $(FASTJETLIBS) 

.PHONY:clean
clean:
	-rm *.o $(EXE)
