## ugly hack to get linking done by CXX
LINK = $(LIBTOOL) --tag=CXX --mode=link $(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) \
        $(AM_LDFLAGS) $(LDFLAGS) -o $@

pkglib_LTLIBRARIES = libHerwig.la
libHerwig_la_SOURCES =
libHerwig_la_LDFLAGS = -version-info 4:2:0 -export-dynamic 
libHerwig_la_LDFLAGS += $(THEPEGLDFLAGS)
libHerwig_la_LIBADD = \
$(top_builddir)/Hadronization/libHwHadronization.la \
$(top_builddir)/Models/StandardModel/libHwStandardModel.la \
$(top_builddir)/Decay/libHwDecay.la \
$(top_builddir)/Decay/FormFactors/libHwFormFactor.la \
$(top_builddir)/Decay/Radiation/libHwDecRad.la \
$(top_builddir)/Utilities/libHwUtils.la \
$(top_builddir)/Models/General/libHwModelGenerator.la \
$(top_builddir)/Decay/General/libHwGeneralDecay.la \
$(top_builddir)/MatrixElement/General/libHwGeneralME.la \
$(top_builddir)/Decay/WeakCurrents/libHwWeakCurrent.la \
$(THEPEGLIB) -ldl

if WANT_LOOPTOOLS
libHerwig_la_LIBADD += \
$(top_builddir)/Looptools/libHwLooptools.la \
$(top_builddir)/Looptools/libHwLooptoolsHelper.la
endif

all-local: done-all-links

done-all-links:
	if test ! -L Herwig++ ; then $(LN_S) . Herwig++ ; fi
	find $(top_builddir) -name '*.so' -not -name 'lib*' -not -path '$(top_builddir)/lib/*' -exec $(LN_S) \{\} \;
	echo "stamp" > done-all-links

clean-local:
	rm -f *.so done-all-links
