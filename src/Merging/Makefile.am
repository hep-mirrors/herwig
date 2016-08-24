BUILT_SOURCES = done-all-links

Matchboxdir = ${pkgdatadir}/Merging

INPUTFILES = \
3LO-LEP.in \
2LO-LEP.in \
MergingLEP.in 


dist_Matchbox_DATA = $(INPUTFILES)

CLEANFILES = done-all-links

done-all-links: $(INPUTFILES)
	@echo "Linking input files"
	@for i in $(INPUTFILES); do \
	if test -f $(srcdir)/$$i -a ! -e $$i; then \
	$(LN_S) -f $(srcdir)/$$i; fi; done
	@touch done-all-links
