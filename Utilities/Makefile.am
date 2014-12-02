SUBDIRS = XML Statistics
noinst_LTLIBRARIES = libHwUtils.la libHwRunDirectories.la
libHwUtils_la_SOURCES = \
EnumParticles.h \
Interpolator.tcc Interpolator.h \
Kinematics.cc Kinematics.h \
Smearing.cc Smearing.h \
Maths.h Maths.cc \
StandardSelectors.cc StandardSelectors.h\
Histogram.cc Histogram.fh Histogram.h \
GaussianIntegrator.cc GaussianIntegrator.h \
GaussianIntegrator.tcc \
Statistic.h HerwigStrategy.cc HerwigStrategy.h \
GSLIntegrator.h GSLIntegrator.tcc \
GSLBisection.h GSLBisection.tcc GSLHelper.h

libHwRunDirectories_la_SOURCES = \
RunDirectories.h RunDirectories.cc

libHwUtils_la_LIBADD = \
XML/libHwXML.la \
Statistics/libHwStatistics.la

CLEANFILES=
include $(srcdir)/Makefile.am.versionstring
