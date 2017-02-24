SUBDIRS = XML Statistics
noinst_LTLIBRARIES = libHwUtils.la

libHwUtils_la_SOURCES = \
EnumParticles.h \
Interpolator.tcc Interpolator.h \
Kinematics.cc Kinematics.h \
Maths.h Maths.cc \
StandardSelectors.cc StandardSelectors.h\
Histogram.cc Histogram.fh Histogram.h \
GaussianIntegrator.cc GaussianIntegrator.h \
GaussianIntegrator.tcc \
Statistic.h HerwigStrategy.cc HerwigStrategy.h \
GSLIntegrator.h GSLIntegrator.tcc \
GSLBisection.h GSLBisection.tcc GSLHelper.h \
expm-1.h

nodist_libHwUtils_la_SOURCES = hgstamp.inc
BUILT_SOURCES = hgstamp.inc
CLEANFILES = hgstamp.inc

HGVERSION := $(shell hg -R $(top_srcdir) parents --template '"Herwig {node|short} ({branch})"' 2> /dev/null || echo \"$(PACKAGE_STRING)\" || true )

.PHONY: update_hgstamp
hgstamp.inc: update_hgstamp
	@[ -f $@ ] || touch $@
	@echo '$(HGVERSION)' | cmp -s $@ - || echo '$(HGVERSION)' > $@

libHwUtils_la_LIBADD = \
XML/libHwXML.la \
Statistics/libHwStatistics.la


check_PROGRAMS = utilities_test
utilities_test_SOURCES = \
tests/utilitiesTestsMain.cc \
tests/utilitiesTestsGlobalFixture.h \
tests/utilitiesTestsKinematics.h \
tests/utilitiesTestMaths.h \
tests/utilitiesTestsStatistic.h
utilities_test_LDADD = $(BOOST_UNIT_TEST_FRAMEWORK_LIBS) $(BOOST_FILESYSTEM_LIBS) $(BOOST_SYSTEM_LIBS) $(THEPEGLIB) -ldl libHwUtils.la
utilities_test_LDFLAGS = $(AM_LDFLAGS) -export-dynamic $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS) $(BOOST_SYSTEM_LDFLAGS) $(BOOST_FILESYSTEM_LDFLAGS) $(THEPEGLDFLAGS)
utilities_test_CPPFLAGS = $(AM_CPPFLAGS) $(BOOST_CPPFLAGS) 
TESTS = utilities_test
