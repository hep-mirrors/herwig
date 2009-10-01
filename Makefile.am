SUBDIRS = include \
	Utilities PDT Decay PDF Models \
	Shower Hadronization MatrixElement \
	UnderlyingEvent Analysis Looptools \
	lib src Doc Contrib Tests

EXTRA_DIST = GUIDELINES

DISTCHECK_CONFIGURE_FLAGS = --enable-debug --with-thepeg=$(THEPEGPATH)

ACLOCAL_AMFLAGS = -I m4

DISTCLEANFILES = config.herwig

libclean:
	find . -name '*.la' -print0 | xargs -0 rm -rf
	cd lib && $(MAKE) clean

tests:
	cd Tests && $(MAKE) tests
