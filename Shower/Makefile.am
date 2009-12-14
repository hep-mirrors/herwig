SUBDIRS = SplittingFunctions Couplings Powheg Fortran \
          Default Base .
pkglib_LTLIBRARIES = HwShower.la
HwShower_la_LDFLAGS = -module -version-info 11:1:0
HwShower_la_LIBADD = SplittingFunctions/libHwSplitting.la \
                     Couplings/libHwCouplings.la \
                     Default/libHwDefaultShower.la \
                     ../PDF/libHwRemDecayer.la \
                     ../PDF/libHwMPIPDF.la\
                     Base/libHwShowerBase.la
HwShower_la_SOURCES = ShowerConfig.h \
ShowerHandler.h ShowerHandler.fh ShowerHandler.cc \
UEBase.h UEBase.cc UEBase.fh
