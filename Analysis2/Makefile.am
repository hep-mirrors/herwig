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
pkglib_LTLIBRARIES += HwFastJetFinder.la
HwFastJetFinder_la_LDFLAGS = -module -version-info 3:0:0
HwFastJetFinder_la_CPPFLAGS = $(AM_CPPFLAGS) $(FASTJETINCLUDE)
HwFastJetFinder_la_LIBADD = $(FASTJETLIBS)
HwFastJetFinder_la_SOURCES = \
FastJetFinder.h FastJetFinder.fh FastJetFinder.icc FastJetFinder.cc
endif

# the ktjet interface

if WANT_LIBKTJET
pkglib_LTLIBRARIES += HwKtJetFinder.la
HwKtJetFinder_la_LDFLAGS = -module -version-info 3:0:0
HwKtJetFinder_la_CPPFLAGS = $(AM_CPPFLAGS) $(KTJETINCLUDE) -I$(KTJETPATH)
HwKtJetFinder_la_LIBADD = $(KTJETLIBS) -L$(CLHEPPATH)  $(CLHEPLIB)
HwKtJetFinder_la_SOURCES = \
KtJetFinder.h KtJetFinder.fh KtJetFinder.icc KtJetFinder.cc
endif
