CC = g++
FLAGS = -ansi -pedantic -Wall -g 
SOFLAGS = -fpic -rdynamic

# Change these paths to point to your installations
HERWIGPATH = /usr/local
THEPEGPATH = /usr/local
KTJETPATH = /usr/local/ktjet
CLHEPPATH = /usr/local

INCLUDE = -I$(THEPEGPATH) -I$(HERWIGPATH) -I$(KTJETPATH) -I$(CLHEPPATH)/include
SOURCES := $(wildcard *.cc)
OBJECTS := $(SOURCES:.cc=.o)
PROGS := *** Your program here ***

HWLIBS = -lHwHadronization -lHwUtils -lHwDecay -lHwShower -lHwME
THEPEGLIBS = -lThePEGRepo -lThePEGHandlers -lThePEGPDF -lThePEGME -lThePEGEvent -lThePEGPDT -lThePEGSM -lThePEGInter -lThePEG -lThePEGMEQCD
DLLIB = -ldl
CLHEPLIB = -lCLHEP
LIBS = -L$(THEPEGPATH)/ThePEG/lib -L$(CLHEPPATH)/lib -L$(HERWIGPATH)/Herwig++/lib $(THEPEGLIBS) $(HWLIBS) $(DLLIB) $(CLHEPLIB)

RPATH = $(THEPEGPATH)/ThePEG/lib:$(HERWIGPATH)/Herwig++/lib

all: $(PROGS) 

%.o: %.cc
	$(CC) -c $(FLAGS) $(INCLUDE) $<

***Your program here ***: $(OBJECTS)
	export LD_RUN_PATH=$(RPATH); $(CC) $(SOFLAGS) -o $@ $(LIBS) $(OBJECTS)
