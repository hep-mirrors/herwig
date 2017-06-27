pkglib_LTLIBRARIES = HwUA5.la
HwUA5_la_SOURCES = UA5Handler.cc  UA5Handler.h UA5Handler.icc
HwUA5_la_LDFLAGS= $(AM_LDFLAGS) -module -version-info 6:0:0

pkglib_LTLIBRARIES += HwMPI.la
HwMPI_la_SOURCES = MPISampler.cc MPISampler.h MPISampler.icc \
	MPISampler.fh MPIHandler.cc \
	MPIHandler.fh MPIHandler.h \
	ProcessHandler.cc ProcessHandler.icc \
	ProcessHandler.fh ProcessHandler.h \
	MPIXSecReweighter.h MPIXSecReweighter.cc \
	stat.h
HwMPI_la_LDFLAGS= $(AM_LDFLAGS) -module -version-info 14:0:0
HwMPI_la_LIBADD= $(GSLLIBS)
HwMPI_la_CPPFLAGS= $(AM_CPPFLAGS) $(GSLINCLUDE)
