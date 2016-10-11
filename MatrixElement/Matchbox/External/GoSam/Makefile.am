EXTRA_DIST = \
GoSamAmplitude.h \
GoSamAmplitude.cc \
gosam2herwig.py

install-exec-local:
	$(install_sh_SCRIPT) $(srcdir)/gosam2herwig.py $(DESTDIR)$(bindir)/gosam2herwig

uninstall-local:
	rm -f  $(DESTDIR)$(bindir)/gosam2herwig
