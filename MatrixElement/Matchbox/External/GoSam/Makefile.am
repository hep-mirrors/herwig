pkglib_LTLIBRARIES = HwMatchboxGoSam.la
HwMatchboxGoSam_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 11:1:0

HwMatchboxGoSam_la_SOURCES = \
GoSamAmplitude.h GoSamAmplitude.cc

HwMatchboxGoSam_la_CPPFLAGS = $(AM_CPPFLAGS) \
-DHERWIG_BINDIR="\"$(bindir)\"" \
-DHERWIG_PKGDATADIR="\"$(pkgdatadir)\"" \
-DGOSAM_PREFIX="\"$(GOSAMPREFIX)\""


install-exec-local:
	$(install_sh_SCRIPT) $(srcdir)/GoSamHelper.py $(DESTDIR)$(bindir)/GoSamHelper.py

uninstall-local:
	rm -f  $(DESTDIR)$(bindir)/GoSamHelper.py

