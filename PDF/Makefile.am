pkglib_LTLIBRARIES = MRST.la
MRST_la_SOURCES = BaryonRemnants.cc BaryonRemnants.h BaryonRemnants.icc MRST99_1.cc MRST99_1.h MRST.cc MRSTData.cc MRSTData.h MRST.h
MRST_la_LDFLAGS = -module

EXTRA_DIST = mrst

install-data-local:
	for i in `find $(srcdir)/mrst -name '*.dat'`; \
	do \
	$(install_sh_DATA) $$i $(DESTDIR)$(pkgdatadir)/PDF/$${i#$(srcdir)/}; \
	done

uninstall-local:
	rm -rf $(DESTDIR)$(pkgdatadir)/PDF

dist-hook:
	rm -rf `find $(distdir)/mrst -name '.svn' -or -name '.cvsignore' -or -name 'CVS'`
