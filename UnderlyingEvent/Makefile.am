pkglib_LTLIBRARIES = HwMPI.la
HwMPI_la_SOURCES = MPISampler.cc MPISampler.h MPISampler.icc \
	MPISampler.fh MPIHandler.cc \
	MPIHandler.fh MPIHandler.h \
	ProcessHandler.cc ProcessHandler.icc \
	ProcessHandler.fh ProcessHandler.h \
	MPIXSecReweighter.h MPIXSecReweighter.cc \
	stat.h
HwMPI_la_LDFLAGS= $(AM_LDFLAGS) -module -version-info 16:0:0
HwMPI_la_LIBADD= $(GSLLIBS)
HwMPI_la_CPPFLAGS= $(AM_CPPFLAGS) $(GSLINCLUDE)
