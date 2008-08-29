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

default : LeptonME.so

%.o : %.cc %.h
	$(CXX) -fPIC $(CPPFLAGS) $(INCLUDE) $(CXXFLAGS) -c -shared $< -o $@  

LeptonME.so: MEee2Higgs2SM.o  MEee2VectorMeson.o  MEee2Z.o
	$(CXX) -fPIC $(CPPFLAGS) $(INCLUDE) $(CXXFLAGS) \
	MEee2Higgs2SM.o  MEee2VectorMeson.o  MEee2Z.o -shared -o LeptonME.so

clean:
	rm -f $(ALLCCFILES:.cc=.o) LeptonME.so
