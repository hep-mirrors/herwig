pkglib_LTLIBRARIES = HwMatchboxGoSam.la
HwMatchboxGoSam_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 13:0:0

HwMatchboxGoSam_la_SOURCES = \
GoSamAmplitude.h GoSamAmplitude.cc

HwMatchboxGoSam_la_CPPFLAGS = $(AM_CPPFLAGS) \
-DHERWIG_BINDIR="\"$(bindir)\"" \
-DHERWIG_PKGDATADIR="\"$(pkgdatadir)\"" \
-DGOSAM_PREFIX="\"$(GOSAMPREFIX)\""

EXTRA_DIST = gosam2herwig.py

install-exec-local:
	$(install_sh_SCRIPT) $(srcdir)/gosam2herwig.py $(DESTDIR)$(bindir)/gosam2herwig

uninstall-local:
	rm -f  $(DESTDIR)$(bindir)/gosam2herwig
