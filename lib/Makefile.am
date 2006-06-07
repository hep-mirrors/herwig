## ugly hack to get linking done by CXX
LINK = $(LIBTOOL) --tag=CXX --mode=link $(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) \
        $(AM_LDFLAGS) $(LDFLAGS) -o $@

pkglib_LTLIBRARIES = libHerwig.la
libHerwig_la_SOURCES =
libHerwig_la_LIBADD = \
$(top_builddir)/Hadronization/libHwHadronization.la \
$(top_builddir)/Models/StandardModel/libHwStandardModel.la \
$(top_builddir)/Decay/libHwDecay.la \
$(top_builddir)/Helicity/Correlations/libHwCorrelations.la \
$(top_builddir)/Helicity/WaveFunction/libHwWaveFunction.la \
$(top_builddir)/Utilities/libHwUtils.la \
$(top_builddir)/Helicity/Vertex/libHwVertex.la 

if WANT_LOOPTOOLS
libHerwig_la_LIBADD += \
$(top_builddir)/Looptools/libHwLooptools.la
endif

if WANT_NEWDECAYERS
libHerwig_la_LIBADD += \
$(top_builddir)/Decay/FormFactors/libHwFormFactor.la \
$(top_builddir)/Decay/WeakCurrents/libHwWeakCurrent.la \
$(top_builddir)/Decay/Radiation/libHwDecRad.la 
endif

all-local: done-all-links

done-all-links:
	find $(top_builddir) -name '*.so' -not -name 'lib*' -not -path '$(top_builddir)/lib/*' -exec $(LN_S) \{\} \;
	echo "stamp" > done-all-links

clean-local:
	rm -f *.so done-all-links
