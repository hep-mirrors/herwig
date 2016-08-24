BUILT_SOURCES = done-all-links

Matchboxdir = ${pkgdatadir}/Merging

INPUTFILES = \
2LO-LEP.in \
2LO-LHC-Z.in \
3LO-LHC-Z.in \
4LO-LHC-Z.in \
2LO-LHC-W.in \
3LO-LHC-W.in \
3LO-LEP.in \
MergingLEP.in \
MergingLHC.in



dist_Matchbox_DATA = $(INPUTFILES)

CLEANFILES = done-all-links

done-all-links: $(INPUTFILES)
	@echo "Linking input files"
	@for i in $(INPUTFILES); do \
	if test -f $(srcdir)/$$i -a ! -e $$i; then \
	$(LN_S) -f $(srcdir)/$$i; fi; done
	@touch done-all-links
