SUBDIRS = SplittingFunctions Couplings \
          Default Base .
pkglib_LTLIBRARIES = HwShower.la
HwShower_la_LDFLAGS = -module
HwShower_la_LIBADD = ./SplittingFunctions/libHwSplitting.la \
                        ./Couplings/libHwCouplings.la \
                        ./Default/libHwDefaultShower.la \
                        ./Base/libHwShowerBase.la
HwShower_la_SOURCES = ShowerConfig.h \
ShowerHandler.h ShowerHandler.fh ShowerHandler.icc ShowerHandler.cc
