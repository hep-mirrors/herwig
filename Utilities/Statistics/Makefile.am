noinst_LTLIBRARIES = libHwStatistics.la
libHwStatistics_la_SOURCES = \
Bin.h \
Counter.h \
CrossSections.h \
Distribution.h \
EventContribution.h \
Histogram.h \
Run.h \
Bin.cc \
Counter.cc \
CrossSections.cc \
Distribution.cc \
Histogram.cc \
Run.cc

bin_PROGRAMS = combineDistributions combineRuns makeDistributions

combineDistributions_SOURCES = \
combineDistributions.cc

combineDistributions_LDADD = \
$(top_builddir)/Utilities/XML/libHwXML.la \
$(top_builddir)/Utilities/Statistics/libHwStatistics.la

combineRuns_SOURCES = \
combineRuns.cc

combineRuns_LDADD = \
$(top_builddir)/Utilities/XML/libHwXML.la \
$(top_builddir)/Utilities/Statistics/libHwStatistics.la

makeDistributions_SOURCES = \
makeDistributions.cc

makeDistributions_LDADD = \
$(top_builddir)/Utilities/XML/libHwXML.la \
$(top_builddir)/Utilities/Statistics/libHwStatistics.la
