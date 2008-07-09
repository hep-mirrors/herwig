SUBDIRS = HardGenerators Matching .
pkglib_LTLIBRARIES = HwPowhegShower.la
libHwPowhegShower_la_LIBADD = HardGenerators/libHwPowhegHard.la
HwPowhegShower_la_SOURCES = \
PowhegEvolver.h PowhegEvolver.fh PowhegEvolver.icc PowhegEvolver.cc\
HardestEmissionGenerator.h HardestEmissionGenerator.fh \
HardestEmissionGenerator.icc HardestEmissionGenerator.cc\
HardTree.h HardTree.icc HardTree.fh HardTree.cc

HwPowhegShower_la_LDFLAGS = -module
HwPowhegShower_la_CPPFLAGS=$(AM_CPPFLAGS) $(GSLINCLUDE)
