pkglib_LTLIBRARIES = 

noinst_LTLIBRARIES = libHwPDF.la 
libHwPDF_la_SOURCES = \
PartonExtractor.cc PartonExtractor.h \
PartonExtractor.icc PartonExtractor.fh

pkglib_LTLIBRARIES += HwMRST.la
HwMRST_la_SOURCES = \
MRST.cc MRST.h MRST.icc 
## add this to produce tests of the PDFs
## HwMRST_la_CPPFLAGS=$(AM_CPPFLAGS) -DMRST_TESTING
HwMRST_la_LDFLAGS = -module -version-info 4:0:0

pkglib_LTLIBRARIES += HwLeptonPDF.la
HwLeptonPDF_la_SOURCES = \
WeiszackerWilliamsPDF.cc WeiszackerWilliamsPDF.h \
LeptonRemnant.cc LeptonRemnant.h LeptonRemnant.icc LeptonRemnant.fh
HwLeptonPDF_la_LDFLAGS = -module -version-info 3:0:0


EXTRA_DIST = mrst

pkglib_LTLIBRARIES += HwRemDecayer.la
HwRemDecayer_la_SOURCES = \
HwRemDecayer.h HwRemDecayer.icc HwRemDecayer.cc HwRemDecayer.fh

HwRemDecayer_la_LDFLAGS = -module -version-info 4:0:0

pkglib_LTLIBRARIES += HwMPIPDF.la
HwMPIPDF_la_SOURCES = \
MPIPDF.h MPIPDF.cc MPIPDF.fh MPIPDF.icc \
MinBiasPDF.h MinBiasPDF.cc MinBiasPDF.fh

HwMPIPDF_la_LDFLAGS = -module -version-info 3:0:0

pkglib_LTLIBRARIES += HwSatPDF.la
HwSatPDF_la_SOURCES = \
SatPDF.h SatPDF.cc SatPDF.fh

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
