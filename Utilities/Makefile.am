SUBDIRS = Histogram2 .
noinst_LTLIBRARIES = libHwUtils.la
libHwUtils_la_LIBADD = Histogram2/libHwHistogram2.la
libHwUtils_la_SOURCES = \
EnumParticles.h \
Interpolator.tcc Interpolator.h \
Kinematics.cc Kinematics.h Kinematics.icc \
Smearing.cc Smearing.h Smearing.icc \
Maths.h Maths.icc Maths.cc \
StandardSelectors.cc StandardSelectors.h StandardSelectors.icc \
Histogram.cc Histogram.fh Histogram.h Histogram.icc \
GaussianIntegrator.cc GaussianIntegrator.h \
GaussianIntegrator.icc GaussianIntegrator.tcc \
Statistic.h HerwigStrategy.cc HerwigStrategy.h \
GSLIntegrator.h GSLIntegrator.tcc \
GSLBisection.h GSLBisection.tcc GSLHelper.h

########## SVN version string. Only use in cc files! ####################
EXTRA_DIST = versionstring.h.in
BUILT_SOURCES = versionstring.h
if USE_SVNVERSION # # # # # #
version:
	@svnversion -n $(top_srcdir) > version.tmp.new

version.tmp: version
	@if ! diff -q version.tmp{.new,}; then mv -f version.tmp{.new,}; fi

versionstring.h: version.tmp versionstring.h.in
	@echo "creating versionstring.h"
	@sed s%@HERWIGVERSIONSTRING@%"$(PACKAGE_STRING) `cd $(top_srcdir); svn info | grep URL | sed -e 's@.*herwig/@@' -e 's@branches/@@'` r`cat version.tmp`"% \
	  $(srcdir)/versionstring.h.in > versionstring.h
else # # # # # # # # # # # # #
versionstring.h:
	@echo "creating versionstring.h"
	@sed s%@HERWIGVERSIONSTRING@%"$(PACKAGE_STRING)"% \
	  $(srcdir)/versionstring.h.in > versionstring.h
	@touch versionstring.h
endif # # # # # # # # # # # # # 
########################################################################

CLEANFILES = versionstring.h version.tmp version.tmp.new
