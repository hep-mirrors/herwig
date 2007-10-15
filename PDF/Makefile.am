pkglib_LTLIBRARIES =
pkglib_LTLIBRARIES += HwMRST.la
HwMRST_la_SOURCES = \
MRST.cc MRST.h MRST.icc \
PartonExtractor.h PartonExtractor.icc PartonExtractor.cc PartonExtractor.fh
# add this to produce tests of the PDFs
#HwMRST_la_CPPFLAGS=$(AM_CPPFLAGS) -DMRST_TESTING
HwMRST_la_LDFLAGS = -module -version-info 2:0:0

EXTRA_DIST = mrst

pkglib_LTLIBRARIES += HwRemDecayer.la
HwRemDecayer_la_SOURCES = \
HwRemDecayer.h HwRemDecayer.icc HwRemDecayer.cc HwRemDecayer.fh \
ForcedSplitting.h ForcedSplitting.icc ForcedSplitting.cc ForcedSplitting.fh

HwRemDecayer_la_LDFLAGS = -module -version-info 2:0:0

pkglib_LTLIBRARIES += HwMPIPDF.la
HwMPIPDF_la_SOURCES = \
MPIPDF.h MPIPDF.icc MPIPDF.cc MPIPDF.fh

HwMPIPDF_la_LDFLAGS = -module -version-info 2:0:0

if WANT_LOCAL_PDF
install-data-local:
	for i in `find $(srcdir)/mrst -name '*.dat'`; \
	do \
	$(install_sh_DATA) $$i $(DESTDIR)$(pkgdatadir)/PDF/$${i#$(srcdir)/}; \
	done

uninstall-local:
	rm -rf $(DESTDIR)$(pkgdatadir)/PDF
endif

dist-hook:
	rm -rf `find $(distdir)/mrst -name '.svn' -or -name '.cvsignore' -or -name 'CVS'`
