pkglib_LTLIBRARIES = HwMRST.la
HwMRST_la_SOURCES = \
MRST.cc MRST.h MRST.icc 
## add this to produce tests of the PDFs
## HwMRST_la_CPPFLAGS=$(AM_CPPFLAGS) -DMRST_TESTING
HwMRST_la_LDFLAGS = -module -version-info 4:0:0

EXTRA_DIST = mrst

## bound into HwShower.la
noinst_LTLIBRARIES = libHwRemDecayer.la
libHwRemDecayer_la_SOURCES = \
HwRemDecayer.h HwRemDecayer.icc HwRemDecayer.cc HwRemDecayer.fh

## bound into HwShower.la
noinst_LTLIBRARIES += libHwMPIPDF.la
libHwMPIPDF_la_SOURCES = \
MPIPDF.h MPIPDF.cc MPIPDF.fh \
MinBiasPDF.h MinBiasPDF.cc MinBiasPDF.fh

pkglib_LTLIBRARIES += HwSatPDF.la
HwSatPDF_la_SOURCES = \
SatPDF.h SatPDF.cc

HwSatPDF_la_LDFLAGS = -module -version-info 1:0:0

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
