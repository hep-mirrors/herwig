#LIBOBJECT := libHwInterfaces.so
LIBOBJECTS := libHwKtJet.so libHwEvtShapes.so

SOURCES := $(wildcard *.cc)

OBJECTS := $(SOURCES:.cc=.o)

#PROGS   := $(SOURCES:.cc=.exe)

HWLIBS = -lHw -lHwHadronization -lHwShower -lHwDecay

EXTRALIBS = -lEventShapes

KTJETPATH = /usera/gieseke/work/mc/ktjet

JETPATH = /usera/gieseke/work/mc/Jet_EventShape_Finders

#EXTRAINCLUDES = -I$(AMEGICPATH)/Include -I$(KTJETPATH) -I$(JETPATH)  -I$(ROOTPATH)/include
EXTRAINCLUDES = -I$(KTJETPATH) -I$(JETPATH)  -I$(ROOTSYS)/include

.PHONY: Lib

%.o: %.cc
	$(CXX) -c -ansi $(CXXEXT) $(COMPILER_FLAG) $(CXXDEBUGFLAG) $(CXXOPTFLAG) $(DEFS) $(INCLUDES) $(EXTRAINCLUDES) $<

include ../Config/Makefile.common

CXXDEBUGFLAG = -g

lib: $(LIBOBJECTS)

%.exe:

#depend:
#	$(SRCDIR)/Config/g++dep $(INCLUDES) $(EXTRAINCLUDES) $(SOURCES) $(CHECKSOURCES)

#libHwAmegic.so: $(OBJECTS)
#	$(CXX) $(CXXFLAGS) $(CXXDEBUGFLAG) $(CXXOPTFLAG) $(CXXDLFLAGS) AmegicInterface.o -o $@

libHwKtJet.so: KtJetInterface.o 
	$(CXX) $(CXXFLAGS) $(CXXDEBUGFLAG) $(CXXOPTFLAG) $(CXXDLFLAGS) KtJetInterface.o -o $@

libHwEvtShapes.so: EvtShapes.o
	$(CXX) $(CXXFLAGS) $(CXXDEBUGFLAG) $(CXXOPTFLAG) $(CXXDLFLAGS) EvtShapes.o -o $@

# DO NOT DELETE THIS LINE -- g++dep uses it.
# DO NOT PUT ANYTHING AFTER THIS LINE, IT WILL GO AWAY.


# IF YOU PUT ANYTHING HERE IT WILL GO AWAY
