EXE=mymain06
GCC=g++

PYTHIAINC=${shell ~/Tools/pythia8310/bin/pythia8-config --includedir}
PYTHIALIBS=${shell ~/Tools/pythia8310/bin/pythia8-config --libdir}

FASTJETINC=${shell ~/Tools/fastjet-install/bin/fastjet-config --cxxflags}
FASTJETLIBS=${shell ~/Tools/fastjet-install/bin/fastjet-config --libs}

all: $(EXE)

mymain06.o: mymain06.cc
	$(GCC) -I$(PYTHIAINC) $(FASTJETINC) -c mymain06.cc

$(EXE): mymain06.o
	$(GCC) -O3 -I$(PYTHIAINC) $(FASTJETINC) -o $(EXE) mymain06.o $(PYTHIALIBS)/libpythia8.a -ldl $(FASTJETLIBS)

.PHONY:clean
clean:
	-rm *.o $(EXE)
