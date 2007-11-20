SUBDIRS = SplittingFunctions Couplings \
          Default CKKW Base .
pkglib_LTLIBRARIES = HwShower.la
HwShower_la_LDFLAGS = -module -version-info 4:0:0
HwShower_la_LIBADD = SplittingFunctions/libHwSplitting.la \
                     Couplings/libHwCouplings.la \
                     Default/libHwDefaultShower.la \
		     CKKW/libHwCKKW.la \
                     Base/libHwShowerBase.la \
		     $(GSLLIBS)
HwShower_la_SOURCES = ShowerConfig.h \
ShowerHandler.h ShowerHandler.fh ShowerHandler.icc ShowerHandler.cc

HwShower_la_CPPFLAGS=$(AM_CPPFLAGS) $(GSLINCLUDE)
