SUBDIRS = Data .

pkglib_LTLIBRARIES = HwAnalysis2.la
HwAnalysis2_la_LDFLAGS = -module -version-info 3:1:0
HwAnalysis2_la_CPPFLAGS = $(AM_CPPFLAGS)
HwAnalysis2_la_LIBADD =
HwAnalysis2_la_SOURCES = \
JetFinder.h JetFinder.fh JetFinder.icc JetFinder.cc \
Analysis2Base.h Analysis2Base.fh Analysis2Base.icc Analysis2Base.cc \
GnuplotOutput.h GnuplotOutput.fh GnuplotOutput.icc GnuplotOutput.cc \
EventShapes2.h EventShapes2.fh EventShapes2.icc EventShapes2.cc \
YMerge.h YMerge.fh YMerge.icc YMerge.cc \
FourJetCorrelations.h FourJetCorrelations.fh \
FourJetCorrelations.icc FourJetCorrelations.cc

# the fastjet interface

if WANT_LIBFASTJET
pkglib_LTLIBRARIES += HwFastJetInterface.la
HwFastJetInterface_la_LDFLAGS = -module -version-info 3:0:0
HwFastJetInterface_la_CPPFLAGS = $(AM_CPPFLAGS) $(FASTJETINCLUDE)
HwFastJetInterface_la_LIBADD = $(FASTJETLIBS)
HwFastJetInterface_la_SOURCES = \
FastJetInterface.h FastJetInterface.fh FastJetInterface.icc FastJetInterface.cc
endif

# the ktjet interface

if WANT_LIBKTJET
pkglib_LTLIBRARIES += HwKtJetInterface.la
HwKtJetInterface_la_LDFLAGS = -module -version-info 3:0:0
HwKtJetInterface_la_CPPFLAGS = $(AM_CPPFLAGS) $(KTJETINCLUDE) -I$(KTJETPATH)
HwKtJetInterface_la_LIBADD = $(KTJETLIBS) -L$(CLHEPPATH)  $(CLHEPLIB)
HwKtJetInterface_la_SOURCES = \
KtJetInterface.h KtJetInterface.fh KtJetInterface.icc KtJetInterface.cc
endif
