pkglib_LTLIBRARIES = HwMRST.la
HwMRST_la_SOURCES = \
MRST.cc MRST.h MRST.icc
## add this to produce tests of the PDFs
## HwMRST_la_CPPFLAGS=$(AM_CPPFLAGS) -DMRST_TESTING
HwMRST_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 8:0:0

EXTRA_DIST = mrst diffraction

pkglib_LTLIBRARIES += HwPomeronPDF.la
HwPomeronPDF_la_SOURCES = \
PomeronPDF.cc PomeronPDF.h
## add this to produce tests of the PDFs
## HwDIFFRACTIVEPDF_la_CPPFLAGS=$(AM_CPPFLAGS) -DDIFFRACTIVEPDF_TESTING
HwPomeronPDF_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 2:0:0

pkglib_LTLIBRARIES += HwReggeonPDF.la
HwReggeonPDF_la_SOURCES = \
ReggeonPDF.cc ReggeonPDF.h
## add this to produce tests of the PDFs
## HwDIFFRACTIVEPDF_la_CPPFLAGS=$(AM_CPPFLAGS) -DDIFFRACTIVEPDF_TESTING
HwReggeonPDF_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 2:0:0

pkglib_LTLIBRARIES += HwPomeronFlux.la
HwPomeronFlux_la_SOURCES = \
PomeronFlux.h PomeronFlux.cc
HwPomeronFlux_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 3:0:0

## bound into HwShower.la
noinst_LTLIBRARIES = libHwRemDecayer.la
libHwRemDecayer_la_SOURCES = \
HwRemDecayer.h HwRemDecayer.cc HwRemDecayer.fh

## bound into HwShower.la
noinst_LTLIBRARIES += libHwMPIPDF.la
libHwMPIPDF_la_SOURCES = \
MPIPDF.h MPIPDF.cc MPIPDF.fh \
MinBiasPDF.h MinBiasPDF.cc MinBiasPDF.fh

pkglib_LTLIBRARIES += HwSatPDF.la
HwSatPDF_la_SOURCES = \
SatPDF.h SatPDF.cc

HwSatPDF_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 5:0:0

pkglib_LTLIBRARIES += HwIncomingPhotonEvolver.la
HwIncomingPhotonEvolver_la_SOURCES = \
IncomingPhotonEvolver.h IncomingPhotonEvolver.cc

HwIncomingPhotonEvolver_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 5:0:0

if WANT_LOCAL_PDF
install-data-local:
	for i in `find $(srcdir)/mrst -name '*.dat'`; \
	do \
	$(install_sh_DATA) $$i $(DESTDIR)$(pkgdatadir)/PDF/$${i#$(srcdir)/}; \
	done

	for i in `find $(srcdir)/diffraction -name '*.data'`; \
	do \
	$(install_sh_DATA) $$i $(DESTDIR)$(pkgdatadir)/PDF/$${i#$(srcdir)/}; \
	done

uninstall-local:
	rm -rf $(DESTDIR)$(pkgdatadir)/PDF
endif

dist-hook:
	rm -rf `find $(distdir)/mrst -name '.svn' -or -name '.cvsignore' -or -name 'CVS'`
	rm -rf `find $(distdir)/diffraction -name '.svn' -or -name '.cvsignore' -or -name 'CVS'`
