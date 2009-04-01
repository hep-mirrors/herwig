###############################################################################
#				DIRECTORIES

#Where LHAPDF is installed
LHAPDFDIR = /usera/seyi/lhapdfnew
#Where the /MCPWNLO/Herwig folder is 
HERWIGDIR = /usera/seyi/MCPWNLO/Herwig

###############################################################################

##                           ROOT DIRECTORIES

#ROOTSYS = /usera/seyi/root

#ROOTLIB = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --libs)

#ROOTINC = $(ROOTSYS)/include/

###############################################################################
#				LIBRARIES

#Now the libraries  -L(directory) -lname NOT preceeded by lib!!!!
#i.e if your lib is called libNAME you need -lNAME

##Path to CLHEP library should be $(CLHEPDIR)/lib/ and check lib name
#CLHEPLIB = -L$(HOME)/ktjet/lib/ -lKtEvent -L$(CLHEPDIR)/lib/ -lCLHEP
#CLHEPLIBNEW =  -L$(CLHEPDIRNEW)/lib/ -lCLHEP
##Path to LHAPDF library should be the same as the LHAPDF directory
LHAPDFLIB = -L$(LHAPDFDIR)/lib -lLHAPDF
LHAPDFINC = $(LHAPDFDIR)/include

#MAC OS users....
#If you had to get an external version of libg2c - i.e it doesn't come with you compiler......We stress, you've probably already got this.
LIBG2C    = /usr/local/lib/
# You would then have to add -L$(LIBG2C) just in front of -lg2c on line 58

################################################################################
#No need to edit anything below this line!

OPT = -funroll-loops -O2 -Wall


FSOURCES  = $(shell ls ./src/*.f)
CXXSOURCES = $(shell ls ./src/*.cxx)
FOBJS     = $(FSOURCES:.f=.o)
CXXOBJS   = $(CXXSOURCES:.cxx=.o)
CC	  = gcc	
FF	  = g77

all: MCHG MCHGPP run_eeZH run_PPVH

%.o:		%.cxx
		$(CC) -c $(OPT) \
		-I$(LHAPDFINC)\
		$^ -o $@

%.o:		%.f
		$(FF) -c $(OPT) $^ -o $@


MCHG:		MCHG.o
		g++ $(OPT) -I./src/ $^\
		$(CLHEPLIB) \
		$(LHAPDFLIB) \
		-L$(LIBG2C) -lg2c -o $@

MCHGPP:		MCHGPP.o
		g++ $(OPT) -I./src/ $^\
		$(CLHEPLIB) \
		$(LHAPDFLIB) \
		-L$(LIBG2C) -lg2c -o $@

run_eeZH:	run_eeZH.o
		g++ $(OPT) -I./src/ $^\
		$(LHAPDFLIB) \
		-L$(LIBG2C) -lg2c -o $@
		mv run_eeZH $(HERWIGDIR)

run_PPVH:	run_PPVH.o
		g++ $(OPT) -I./src/ $^\
		$(LHAPDFLIB) \
		-L$(LIBG2C) -lg2c -o $@
		mv run_PPVH $(HERWIGDIR)



clean:
		rm *.o
		rm MCHG
		rm MCHGPP




