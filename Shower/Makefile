OBJECTS = libHwShower.so 

SOURCES := $(wildcard *.cc)

HEADERS := $(wildcard *.h)

LIBOBJECTS := $(SOURCES:.cc=.o)

INCCLHEP := -lCLHEP

include ../Config/Makefile.common

$(OBJECTS): $(LIBOBJECTS)
	@cd $(SRCDIR)/lib; rm -f $(LIBOBJECTS)
	$(CXX) $(CXXFLAGS) $(CXXDEBUGFLAG) $(CXXOPTFLAG) $(CXXDLFLAGS) $(LIBOBJECTS) $(INCCLHEP) -o $@
	@ln -f -s `pwd`/$@ $(SRCDIR)/lib

# DO NOT DELETE THIS LINE -- g++dep uses it.
