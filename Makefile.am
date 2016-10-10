SUBDIRS = include \
	Utilities PDT Decay PDF Models \
	Shower DipoleShower Hadronization MatrixElement \
	UnderlyingEvent Analysis Looptools Sampling \
	API lib src Doc Contrib Tests

EXTRA_DIST = GUIDELINES

## DISTCHECK_CONFIGURE_FLAGS = --enable-debug --with-thepeg=$(THEPEGPATH) --with-gsl=$(GSLPATH)

ACLOCAL_AMFLAGS = -I m4

DISTCLEANFILES = config.herwig

libclean:
	find . -name '*.la' -print0 | xargs -0 rm -rf
	cd lib && $(MAKE) $(AM_MAKEFLAGS) clean
	cd src && $(MAKE) $(AM_MAKEFLAGS) clean

tests:
	cd Tests && $(MAKE) $(AM_MAKEFLAGS) tests


## ThePEG registration
unregister:
	cd src && $(MAKE) $(AM_MAKEFLAGS) unregister

register:
	cd src && $(MAKE) $(AM_MAKEFLAGS) register
