SUBDIRS = include \
	Utilities PDT Decay PDF Models \
	Shower Hadronization MatrixElement \
	UnderlyingEvent Analysis Looptools Sampling \
	API lib src Doc Contrib Tests

EXTRA_DIST = GUIDELINES cat_with_cpplines

## DISTCHECK_CONFIGURE_FLAGS = --enable-debug --with-thepeg=$(THEPEGPATH) --with-gsl=$(GSLPATH) --with-boost=/mt/data1/herwig

ACLOCAL_AMFLAGS = -I m4

DISTCLEANFILES = config.herwig

libclean:
	find . -name '*.la' -delete
	cd lib && $(MAKE) $(AM_MAKEFLAGS) clean
	cd src && $(MAKE) $(AM_MAKEFLAGS) clean

tests:
	cd Tests && $(MAKE) $(AM_MAKEFLAGS) tests


## ThePEG registration
unregister:
	cd src && $(MAKE) $(AM_MAKEFLAGS) unregister

register:
	cd src && $(MAKE) $(AM_MAKEFLAGS) register
