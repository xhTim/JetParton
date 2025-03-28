EXE=main_string_fragmentation
GCC=g++

PYTHIAINC=${shell /home/wenbin/Downloads/Wenbin_working/Work/Berkeley_work/v2_in_jet/parton_shower_in_pythia8/pythia8310/bin/pythia8-config --includedir}
PYTHIALIBS=${shell /home/wenbin/Downloads/Wenbin_working/Work/Berkeley_work/v2_in_jet/parton_shower_in_pythia8/pythia8310/bin/pythia8-config --libdir}



all: $(EXE)

main_string_fragmentation.o: main_string_fragmentation.cc
	$(GCC) -I$(PYTHIAINC) -c main_string_fragmentation.cc

$(EXE): main_string_fragmentation.o
	$(GCC) -O3 -I$(PYTHIAINC) -o $(EXE) main_string_fragmentation.o $(PYTHIALIBS)/libpythia8.a -ldl 

.PHONY:clean
clean:
	-rm *.o $(EXE)
