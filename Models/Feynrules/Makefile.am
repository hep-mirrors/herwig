EXTRA_DIST = \
exampleInputs \
python

install-data-local:
	for i in `find $(srcdir)/python -type f`; \
	do \
	$(install_sh_DATA) $$i $(DESTDIR)$(pkglibdir)/$${i#$(srcdir)/}; \
	done

install-exec-local:
	$(install_sh_SCRIPT) $(srcdir)/python/ufo2herwig $(DESTDIR)$(bindir)/ufo2herwig
	sed -e's!@PKGLIBDIR@!$(pkglibdir)!' $(srcdir)/python/ufo2herwig > $(DESTDIR)$(bindir)/ufo2herwig

uninstall-local:
	rm -rf $(DESTDIR)$(pkglibdir)/python
	rm -f  $(DESTDIR)$(bindir)/ufo2herwig

## pkgdata_DATA = Makefile-FR
