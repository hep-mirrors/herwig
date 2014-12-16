pkglib_LTLIBRARIES = HwMatchboxGoSam.la
HwMatchboxGoSam_la_LDFLAGS = -module -version-info 11:1:0

HwMatchboxGoSam_la_SOURCES = \
GoSamAmplitude.h GoSamAmplitude.cc


install-exec-local:
	$(install_sh_SCRIPT) $(builddir)/GoSamHelper.py $(DESTDIR)$(bindir)/GoSamHelper.py

uninstall-local:
	rm -f  $(DESTDIR)$(bindir)/GoSamHelper.py