DEPS = fitter.h 
OBJS = VacuumScrubbingFitter.o 

ROOTCFLAGS   := $(shell root-config --cflags) -Wl,--no-as-needed
ROOTLIBS     := $(shell root-config --glibs) -lTreePlayer 

XX := g++

CXXFLAGS = -Wall -O2 $(ROOTCFLAGS)

LDLIBS = $(ROOTLIBS) 
#OTHERLIBS = ./Loader_C.so 

EXEC_FILES = VacuumScrubbingFitter.exe

%.o: %.cxx $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(LDLIBS) 


VacuumScrubbingFitter.exe: $(OBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDLIBS) 


.PHONY: clean

clean:
	rm -f *.o $(OBJS) $(EXEC_FILES)
