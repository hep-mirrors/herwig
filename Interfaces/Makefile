#LIBOBJECT := libHwInterfaces.so
LIBOBJECTS := libHwAmegic.so libHwKtJet.so

SOURCES := $(wildcard *.cc)

OBJECTS := $(SOURCES:.cc=.o)

#PROGS   := $(SOURCES:.cc=.exe)

HWLIBS = -lHw -lHwHadronization -lHwShower -lHwDecay

EXTRALIBS = -lAmegic -lAmplitude -lAmplTools -lISR -lModel -lPhasespace -lProcess -lString -lVector -lZfunctions

AMEGICPATH = /usera/stephens/projects/AMEGIC++-1.0

KTJETPATH = /usera/stephens/projects/ktjet

EXTRAINCLUDES = -I$(AMEGICPATH)/Include -I$(KTJETPATH)

.PHONY: Lib

%.o: %.cc
	$(CXX) -c -ansi $(CXXEXT) $(COMPILER_FLAG) $(CXXDEBUGFLAG) $(CXXOPTFLAG) $(DEFS) $(INCLUDES) $(EXTRAINCLUDES) $<

include ../Config/Makefile.common

CXXDEBUGFLAG = -g

lib: $(LIBOBJECTS)

%.exe:

libHwAmegic.so: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(CXXDEBUGFLAG) $(CXXOPTFLAG) $(CXXDLFLAGS) AmegicInterface.o -o $@

libHwKtJet.so: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(CXXDEBUGFLAG) $(CXXOPTFLAG) $(CXXDLFLAGS) KtJetInterface.o -o $@

# DO NOT DELETE THIS LINE -- g++dep uses it.
# DO NOT PUT ANYTHING AFTER THIS LINE, IT WILL GO AWAY.


# IF YOU PUT ANYTHING HERE IT WILL GO AWAY
