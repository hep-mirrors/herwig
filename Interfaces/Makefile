# Generated automatically from Makefile.in by configure.
LIBOBJECTS := libHwAmegic.so libHwKtJet.so

SOURCES := $(wildcard *.cc)

OBJECTS := $(SOURCES:.cc=.o)

HWLIBS = -lHwUtils -lHwHadronization -lHwShower -lHwDecay -lHwME

AMEGICLIBS = -lAmegic -lAmplitude -lAmplTools -lISR -lModel -lPhasespace -lProcess -lString -lVector -lZfunctions
AMEGICPATH = /usera/stephens/projects/AMEGIC++-1.0
AMEGICINCLUDES = -I$(AMEGICPATH)/Include

KTJETPATH = /var/nwork/pcck/gieseke/mc/ktjet
KTJETINCLUDES = -I$(KTJETPATH)

.PHONY: Lib

AmegicInterface.o: AmegicInterface.cc
	$(CXX) -c $(CXXEXT) $(COMPILER_FLAG) $(CXXDEBUGFLAG) $(CXXOPTFLAG) $(DEFS) $(INCLUDES) $(AMEGICINCLUDES) $<

KtJetInterface.o: KtJetInterface.cc
	$(CXX) -c $(CXXEXT) $(COMPILER_FLAG) $(CXXDEBUGFLAG) $(CXXOPTFLAG) $(DEFS) $(INCLUDES) $(KTJETINCLUDES) $<
include ../Config/Makefile.common

CXXDEBUGFLAG = -g 

lib: $(LIBOBJECTS)

%.exe:

libHwAmegic.so: AmegicInterface.o
	@cd $(SRCDIR)/lib; rm -f $@
	$(CXX) $(CXXFLAGS) $(CXXDEBUGFLAG) $(CXXOPTFLAG) $(CXXDLFLAGS) $< -o $@
	@ln -f -s `pwd`/$@ $(SRCDIR)/lib

libHwKtJet.so: KtJetInterface.o
	@cd $(SRCDIR)/lib; rm -f $@
	$(CXX) $(CXXFLAGS) $(CXXDEBUGFLAG) $(CXXOPTFLAG) $(CXXDLFLAGS) $< -o $@
	@ln -f -s `pwd`/$@ $(SRCDIR)/lib


