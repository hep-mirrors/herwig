SUBDIRS = include \
	Utilities PDT Decay PDF Models \
	Shower Hadronization MatrixElement \
	UnderlyingEvent Interfaces Analysis Looptools \
	lib src Doc Contrib

noinst_HEADERS = Config/HepMCHelper.h

EXTRA_DIST=GUIDELINES

DISTCHECK_CONFIGURE_FLAGS = --enable-debug --with-thepeg=$(THEPEGPATH)

ACLOCAL_AMFLAGS = -I m4

libclean:
	find . -name '*.la' -print0 | xargs -0 rm -rf
	(cd lib && make clean)
