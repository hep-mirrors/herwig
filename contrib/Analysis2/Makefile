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

default : Analysis2.so

%.o : %.cc %.h
	$(CXX) -fPIC $(CPPFLAGS) $(INCLUDE) $(CXXFLAGS) -c -shared $< -o $@  

Analysis2.so: Analysis2Base.o  FourJetCorrelations.o  Histogram2Output.o  YMerge.o \
	      EventShapes2.o   GnuplotOutput.o        JetFinder.o \
              FastJetFinder.o  Histogram2.o           KtJetFinder.o
	$(CXX) -fPIC $(CPPFLAGS) $(INCLUDE) $(CXXFLAGS) \
	Analysis2Base.o  FourJetCorrelations.o  Histogram2Output.o  YMerge.o \
	EventShapes2.o   GnuplotOutput.o        JetFinder.o \
        FastJetFinder.o  Histogram2.o           KtJetFinder.o -shared -o Analysis2.so

clean:
	rm -f $(ALLCCFILES:.cc=.o) Analysis2.so
