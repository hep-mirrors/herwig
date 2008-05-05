SUBDIRS = include \
	Utilities PDT Decay PDF Models \
	Shower Hadronization MatrixElement \
	UnderlyingEvent Interfaces Analysis Analysis2 Looptools \
	lib src Doc

noinst_HEADERS = Config/HepMCHelper_CLHEP.h Config/HepMCHelper_HepMC.h

EXTRA_DIST=GUIDELINES

DISTCHECK_CONFIGURE_FLAGS = --enable-debug --with-thepeg=$(THEPEGPATH)
