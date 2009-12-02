pkglib_LTLIBRARIES = Herwig.la
Herwig_la_SOURCES =
Herwig_la_LIBTOOLFLAGS = --tag=CXX
Herwig_la_LDFLAGS = -module -version-info 6:0:0
Herwig_la_LDFLAGS += $(THEPEGLDFLAGS) $(FCLIBS)
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
$(top_builddir)/MatrixElement/libHwME.la \
$(top_builddir)/Decay/WeakCurrents/libHwWeakCurrent.la \
$(top_builddir)/Looptools/libHwLooptoolsXFC.la \
$(THEPEGLIB) -ldl

all-local: done-all-links

done-all-links: Herwig.la
	if test ! -L Herwig++ ; then $(LN_S) -f . Herwig++ ; fi
	find $(top_builddir) \( -name '*.so.*' -or -name '*.so' \) \
-not -name 'lib*' -not -path '$(top_builddir)/lib/*' -exec $(LN_S) -f \{\} \;
	$(LN_S) -f .libs/Herwig*so* .
	echo "stamp" > done-all-links

clean-local:
	rm -f *.so *.so.* done-all-links
