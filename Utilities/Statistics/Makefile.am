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

combineDistributions_LDADD = \
$(top_builddir)/Utilities/XML/libHwXML.la \
libHwStatistics.la

combineDistributions_SOURCES = \
combineDistributions.cc

combineRuns_LDADD = \
$(top_builddir)/Utilities/XML/libHwXML.la \
libHwStatistics.la

combineRuns_SOURCES = \
combineRuns.cc

makeDistributions_LDADD = \
$(top_builddir)/Utilities/XML/libHwXML.la \
libHwStatistics.la

makeDistributions_SOURCES = \
makeDistributions.cc

