# -*- Makefile -*- (for emacs)

#
# This Makefile is intended for compiling Herwig++ plugins
# You can find plugins here: INSERT URL
#
# This Makefile received very little testing, 
# any bug reports are very welcome!
#

# location of include files
HPPDIR = /usera/seyi/Herwig
THEPEGDIR = /usera/seyi/Herwig

#CLHEPINCLUDE  = 
#ROOTINCLUDE = -I/usera/andreas/root/include
#ROOTCFLAGS  := $(shell /usera/andreas/root/bin/root-config --ldflags)
#ROOTGLIBS   := $(shell /usera/andreas/root/bin/root-config --glibs)
#ROOT = $(ROOTCFLAGS) $(ROOTGLIBS)
INCLUDE = -I$(THEPEGDIR)/include -I$(HPPDIR)/include
#
# C++ flags
# 
CXX = g++
CXXFLAGS = -O2

ALLCCFILES=$(shell echo *.cc)
all : $(ALLCCFILES:.cc=.so)

%.so : %.cc %.h 
	$(CXX) -fPIC $(CPPFLAGS) $(INCLUDE) $(ROOT) $(CXXFLAGS) -shared $< -o $@  

install: 
#	cp LEPEventShapes.so $(HPPDIR)/lib/Herwig++
	cp SimpleILCAnalysis.so $(HPPDIR)/lib/Herwig++
clean:
	rm -f $(ALLCCFILES:.cc=.so)
