BUILT_SOURCES = done-all-links

snippetsdir = ${pkgdatadir}/snippets

INPUTFILES = \
DipoleShowerParameters.in

dist_snippets_DATA = $(INPUTFILES)

CLEANFILES = done-all-links

done-all-links: $(INPUTFILES)
	@echo "Linking input files"
	@for i in $(INPUTFILES); do \
	if test -f $(srcdir)/$$i -a ! -e $$i; then \
	$(LN_S) -f $(srcdir)/$$i; fi; done
	@touch done-all-links
