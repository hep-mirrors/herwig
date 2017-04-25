BUILT_SOURCES = done-all-links

Matchboxdir = ${pkgdatadir}/Merging

INPUTFILES = \
\
LEP91-Analysis.in \
LHC7-Z-Analysis.in \
LHC7-W-Analysis.in \
LHC7-J-Analysis.in \
LHC7-T-Analysis.in \
LHC8-H-Analysis.in \
\
FactorCMWScheme.in \
LinearCMWScheme.in \
Merging-Dipole-FactorCMWSchemeTune.in \
\
LEP-Merging.in \
LHC-W-Merging.in \
LHC-Z-Merging.in \
LHC-H-Merging.in \
LHC-T-Merging.in \
LHC-J-Merging.in
 
dist_Matchbox_DATA = $(INPUTFILES)

CLEANFILES = done-all-links

done-all-links: $(INPUTFILES)
	@echo "Linking input files"
	@for i in $(INPUTFILES); do \
	if test -f $(srcdir)/$$i -a ! -e $$i; then \
	$(LN_S) -f $(srcdir)/$$i; fi; done
	@touch done-all-links
