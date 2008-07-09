SUBDIRS = HardGenerators Matching .
pkglib_LTLIBRARIES = HwNasonShower.la
libHwNasonShower_la_LIBADD = HardGenerators/libHwNasonHard.la
HwNasonShower_la_SOURCES = \
NasonEvolver.h NasonEvolver.fh NasonEvolver.icc NasonEvolver.cc\
HardestEmissionGenerator.h HardestEmissionGenerator.fh \
HardestEmissionGenerator.icc HardestEmissionGenerator.cc\
NasonTree.h NasonTree.icc NasonTree.fh NasonTree.cc

HwNasonShower_la_LDFLAGS = -module
HwNasonShower_la_CPPFLAGS=$(AM_CPPFLAGS) $(GSLINCLUDE)
