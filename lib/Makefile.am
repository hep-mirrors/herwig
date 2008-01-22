## ugly hack to get linking done by CXX
LINK = $(LIBTOOL) --tag=CXX --mode=link $(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) \
        $(AM_LDFLAGS) $(LDFLAGS) -o $@

pkglib_LTLIBRARIES = Herwig.la
Herwig_la_SOURCES =
Herwig_la_LDFLAGS = -module -version-info 4:2:0
Herwig_la_LDFLAGS += $(THEPEGLDFLAGS)
Herwig_la_LIBADD = \
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
Herwig_la_LIBADD += \
$(top_builddir)/Looptools/libHwLooptools.la \
$(top_builddir)/Looptools/libHwLooptoolsHelper.la
Herwig_la_LDFLAGS += $(FLIBS)
endif

all-local: done-all-links

done-all-links: | Herwig.la
	if test ! -L Herwig++ ; then $(LN_S) . Herwig++ ; fi
	find $(top_builddir) \( -name '*.so.*' -or -name '*.so' \) \
-not -name 'lib*' -not -path '$(top_builddir)/lib/*' -exec $(LN_S) \{\} \;
	$(LN_S) .libs/Herwig.so .libs/Herwig.so.* .
	echo "stamp" > done-all-links

clean-local:
	rm -f *.so *.so.* done-all-links
