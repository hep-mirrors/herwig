LIBRARIES = libHw.so

include ../Config/Makefile.common

lib: $(LIBRARIES)

$(LIBRARIES): $(wildcard *.o)
	$(CXX) $(CXXFLAGS) $(CXXDEBUGFLAG) $(CXXOPTFLAG) $(CXXDLFLAGS) $(wildcard *.o) -o $@

# DO NOT DELETE THIS LINE -- g++dep uses it.
# DO NOT PUT ANYTHING AFTER THIS LINE, IT WILL GO AWAY.


# IF YOU PUT ANYTHING HERE IT WILL GO AWAY
