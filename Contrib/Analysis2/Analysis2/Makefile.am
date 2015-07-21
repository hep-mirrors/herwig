SUBDIRS = Data

pkglib_LTLIBRARIES = Analysis2.la
Analysis2_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 0:0:0
Analysis2_la_CPPFLAGS = $(AM_CPPFLAGS)
Analysis2_la_LIBADD = $(top_builddir)/Histogram2/libHistogram2.la
Analysis2_la_SOURCES = \
EventExtractor.h EventExtractor.cc \
JetFinder.h JetFinder.fh JetFinder.icc JetFinder.cc \
Analysis2Base.h Analysis2Base.fh Analysis2Base.icc Analysis2Base.cc \
GnuplotOutput.h GnuplotOutput.fh GnuplotOutput.icc GnuplotOutput.cc \
EventShapes2.h EventShapes2.fh EventShapes2.icc EventShapes2.cc \
YMerge.h YMerge.fh YMerge.icc YMerge.cc \
YJetRates.h YJetRates.cc \
FourJetCorrelations.h FourJetCorrelations.fh \
FourJetCorrelations.icc FourJetCorrelations.cc

# the fastjet interface

# if WITH_LIBFASTJET
# pkglib_LTLIBRARIES += FastJetFinder.la
# FastJetFinder_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 3:0:0
# FastJetFinder_la_CPPFLAGS = $(AM_CPPFLAGS) $(FASTJETINCLUDE)
# FastJetFinder_la_LIBADD = $(FASTJETLIBS)
# FastJetFinder_la_SOURCES = \
# FastJetFinder.h FastJetFinder.fh FastJetFinder.icc FastJetFinder.cc
# endif

# the ktjet interface

if WITH_LIBKTJET
pkglib_LTLIBRARIES += KtJetFinder.la
KtJetFinder_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 3:0:0
KtJetFinder_la_CPPFLAGS = $(AM_CPPFLAGS) $(KTJETINCLUDE) -I$(KTJETPATH)
KtJetFinder_la_LIBADD = $(KTJETLIBS) -L$(CLHEPPATH)  $(CLHEPLIB)
KtJetFinder_la_SOURCES = \
KtJetFinder.h KtJetFinder.fh KtJetFinder.icc KtJetFinder.cc
endif
