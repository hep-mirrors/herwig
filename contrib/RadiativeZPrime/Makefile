# -*- Makefile -*- (for emacs)

#
# This Makefile is intended for compiling Herwig++ plugins
# You can find plugins here: INSERT URL
#
# This Makefile received very little testing, 
# any bug reports are very welcome!
#

# location of include files
THEPEGINCLUDE = -I/home/richardn/montecarlo/c++/ThePEG/trunk
CLHEPINCLUDE  = -I/home/richardn/montecarlo/c++/utilities/include
HERWIGINCLUDE = -I/home/richardn/montecarlo/c++/Herwig++/trunk
INCLUDE = $(THEPEGINCLUDE) $(CLHEPINCLUDE) $(HERWIGINCLUDE)
#
# C++ flags
# 
CXX = g++
CXXFLAGS = -O2

ALLCCFILES=$(shell echo *.cc)

default : RadiativeZPrime.so

%.o : %.cc %.h
	$(CXX) -fPIC $(CPPFLAGS) $(INCLUDE) $(CXXFLAGS) -c -shared $< -o $@  

RadiativeZPrime.so : AnomalousVVVVertex.o  RadiativeZPrimeModel.o  FFZPrimeVertex.o\
		     GammaZPrimeZVertex.o MEqq2ZPrime2ff.o MEqq2ZPrime2ZGamma.o\
                     MEqq2ZPrime2ZGamma2ffGamma.o SimpleZPrimeAnalysis.o
	$(CXX) -fPIC $(CPPFLAGS) $(INCLUDE) $(CXXFLAGS) -shared AnomalousVVVVertex.o \
	RadiativeZPrimeModel.o  FFZPrimeVertex.o GammaZPrimeZVertex.o \
	MEqq2ZPrime2ZGamma2ffGamma.o MEqq2ZPrime2ZGamma.o MEqq2ZPrime2ff.o \
	 SimpleZPrimeAnalysis.o -o RadiativeZPrime.so

clean:
	rm -f $(ALLCCFILES:.cc=.o) RadiativeZPrime.so

Zprime.tar.gz: 
	tar czf ./Zprime.tar.gz *.cc *.icc *.fh *.h Makefile LHC-RadiativeZPrime.in \
	RadiativeZPrime.model