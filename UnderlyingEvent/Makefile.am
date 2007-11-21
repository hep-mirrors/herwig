pkglib_LTLIBRARIES = HwUA5.la
HwUA5_la_SOURCES = UA5Handler.cc  UA5Handler.h UA5Handler.icc
HwUA5_la_LDFLAGS= -module -version-info 3:0:0

pkglib_LTLIBRARIES += HwMPI.la
HwMPI_la_SOURCES = MPISampler.cc MPISampler.h MPISampler.icc \
	MPIHandler.cc MPIHandler.icc MPIHandler.fh MPIHandler.h \
	stat.h
HwMPI_la_LDFLAGS= -module -version-info 3:0:0
HwMPI_la_LIBADD= $(GSLLIBS)
HwMPI_la_CPPFLAGS= $(AM_CPPFLAGS) $(GSLINCLUDE)
