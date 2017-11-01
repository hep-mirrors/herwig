pkglib_LTLIBRARIES = Herwig.la
Herwig_la_SOURCES =
Herwig_la_LIBTOOLFLAGS = --tag=CXX
Herwig_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 23:0:0
Herwig_la_LDFLAGS += $(THEPEGLDFLAGS) $(BOOST_SYSTEM_LDFLAGS) $(BOOST_FILESYSTEM_LDFLAGS) $(FCLIBS)
Herwig_la_LIBADD = \
$(top_builddir)/Hadronization/libHwHadronization.la \
$(top_builddir)/Models/libHwStandardModel.la \
$(top_builddir)/Decay/libHwDecay.la \
$(top_builddir)/Decay/libHwFormFactor.la \
$(top_builddir)/Decay/libHwDecRad.la \
$(top_builddir)/Utilities/libHwUtils.la \
$(top_builddir)/Models/libHwModelGenerator.la \
$(top_builddir)/Decay/General/libHwGeneralDecay.la \
$(top_builddir)/MatrixElement/General/libHwGeneralME.la \
$(top_builddir)/MatrixElement/libHwME.la \
$(top_builddir)/MatrixElement/Reweighters/libHwReweighters.la \
$(top_builddir)/MatrixElement/Matchbox/libHwMatchbox.la \
$(top_builddir)/Decay/libHwWeakCurrent.la \
$(top_builddir)/Looptools/libHwLooptools.la \
$(top_builddir)/Shower/libHwShower.la \
$(THEPEGLIB) $(BOOST_SYSTEM_LIBS) $(BOOST_FILESYSTEM_LIBS) -ldl

dist_noinst_SCRIPTS = fix-osx-path


POSTPROCESSING = done-all-links

if NEED_APPLE_FIXES
POSTPROCESSING += apple-fixes
endif



all-local: $(POSTPROCESSING)

done-all-links: Herwig.la
	find $(top_builddir) \( -name '*.so.*' -or -name '*.so' \) \
  -not -name 'lib*' -not -path '$(top_builddir)/lib/*' \
  -not -path '$(top_builddir)/.hg/*' -exec $(LN_S) -f \{\} \;
	$(LN_S) -f .libs/Herwig*so* .
	echo "stamp" > $@

apple-fixes: fix-osx-path done-all-links
	./$<
	echo "stamp" > $@

clean-local:
	rm -f *.so *.so.* done-all-links apple-fixes
