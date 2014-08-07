EXTRA_DIST = \
exampleInputs \
python

if HAVE_PYTHON

install-data-local:
	for i in `find $(srcdir)/python -type f`; \
	do \
	$(install_sh_DATA) $$i $(DESTDIR)$(pkglibdir)/$${i#$(srcdir)/}; \
	done

install-exec-local:
	$(install_sh_SCRIPT) $(srcdir)/python/ufo2herwig $(DESTDIR)$(bindir)/ufo2herwig
	sed -e's!@PKGLIBDIR@!$(pkglibdir)!' $(srcdir)/python/ufo2herwig > $(DESTDIR)$(bindir)/ufo2herwig
	$(install_sh_SCRIPT) $(srcdir)/python/slha2herwig $(DESTDIR)$(bindir)/slha2herwig

uninstall-local:
	rm -rf $(DESTDIR)$(pkglibdir)/python
	rm -f  $(DESTDIR)$(bindir)/ufo2herwig
	rm -f  $(DESTDIR)$(bindir)/slha2herwig

endif

## pkgdata_DATA = Makefile-FR
