SUBDIRS = HardGenerators Matching .
pkglib_LTLIBRARIES = HwPowhegShower.la
HwPowhegShower_la_LDFLAGS = -module -version-info 1:0:0
HwPowhegShower_la_LIBADD = HardGenerators/libHwPowhegHard.la \
                              Matching/libHwPowhegMatch.la
HwPowhegShower_la_SOURCES = \
PowhegEvolver.h PowhegEvolver.fh PowhegEvolver.cc\
HardestEmissionGenerator.h HardestEmissionGenerator.fh \
HardestEmissionGenerator.icc HardestEmissionGenerator.cc\
HardTree.h HardTree.icc HardTree.fh HardTree.cc

HwPowhegShower_la_LDFLAGS = -module
HwPowhegShower_la_CPPFLAGS=$(AM_CPPFLAGS) $(GSLINCLUDE)
