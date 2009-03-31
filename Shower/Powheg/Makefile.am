SUBDIRS = HardGenerators Matching .
pkglib_LTLIBRARIES = HwPowhegShower.la
HwPowhegShower_la_LDFLAGS = -module -version-info 1:0:0
HwPowhegShower_la_LIBADD = HardGenerators/libHwPowhegHard.la \
                              Matching/libHwPowhegMatch.la
HwPowhegShower_la_SOURCES = \
PowhegEvolver.h PowhegEvolver.fh PowhegEvolver.cc\
QEDEvolver.h QEDEvolver.cc\
HardestEmissionGenerator.h HardestEmissionGenerator.fh \
HardestEmissionGenerator.cc
