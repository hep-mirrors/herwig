OBJECTS = libHwDecay.so 

SOURCES := $(wildcard *.cc)

LIBOBJECTS := $(SOURCES:.cc=.o)

include ../Config/Makefile.common

$(OBJECTS): $(LIBOBJECTS)
	@cd $(SRCDIR)/lib; rm -f $(LIBOBJECTS)
	$(CXX) $(CXXFLAGS) $(CXXDEBUGFLAG) $(CXXOPTFLAG) $(CXXDLFLAGS) $(LIBOBJECTS) -o $@
	@ln -f -s `pwd`/$@ $(SRCDIR)/lib

# DO NOT DELETE THIS LINE -- g++dep uses it.
# DO NOT PUT ANYTHING AFTER THIS LINE, IT WILL GO AWAY.
