SOURCES := $(wildcard *.cc)

OBJECTS := $(SOURCES:.cc=.o)

PROGS   := $(SOURCES:.cc=.exe)

HWLIBS = -lHwHadronization -lHw

.PHONY: Lib

%.o: %.cc
	$(CXX) -c $(CXXFLAGS) $(COMPILER_FLAG) $(CXXDEBUGFLAG) $(CXXOPTFLAG) $(DEFS) $(INCLUDES) $(CXXEXT) $<

include ../Config/Makefile.common

CXXDEBUGFLAG = -g

Lib:
	cd .. ; $(MAKE) -$(MAKEFLAGS) lib

current: $(PROGS) 

$(PROGS): ../lib/libHwHadronization.so ../lib/libHw.so

%.exe: %.o 
	export LD_RUN_PATH=$(RPATH); $(CXX) $(CXXEXT) $(LINKFLAGS) $(CXXDEBUGFLAG) -o $@ $< -L$(SRCDIR)/lib $(HWLIBS) $(LIBS)

# DO NOT DELETE THIS LINE -- g++dep uses it.
