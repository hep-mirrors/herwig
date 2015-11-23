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

bin_PROGRAMS = herwig-combinedistributions herwig-combineruns herwig-makedistributions

herwig_combinedistributions_LDADD = \
$(top_builddir)/Utilities/XML/libHwXML.la \
libHwStatistics.la

herwig_combinedistributions_SOURCES = \
herwig-combinedistributions.cc

herwig_combineruns_LDADD = \
$(top_builddir)/Utilities/XML/libHwXML.la \
libHwStatistics.la

herwig_combineruns_SOURCES = \
herwig-combineruns.cc

herwig_makedistributions_LDADD = \
$(top_builddir)/Utilities/XML/libHwXML.la \
libHwStatistics.la

herwig_makedistributions_SOURCES = \
herwig-makedistributions.cc

