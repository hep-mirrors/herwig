SUBDIRS = include \
	Utilities PDT Decay PDF Models \
	Shower Hadronization MatrixElement \
	UnderlyingEvent Interfaces Analysis Looptools \
	lib src Doc

noinst_HEADERS = Config/Herwig.h Config/HepMCHelper_CLHEP.h \
	Config/HepMCHelper_HepMC.h

EXTRA_DIST=GUIDELINES

DISTCHECK_CONFIGURE_FLAGS = --enable-debug --with-thepeg=$(THEPEGPATH)
