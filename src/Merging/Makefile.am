BUILT_SOURCES = done-all-links

Matchboxdir = ${pkgdatadir}/Merging

INPUTFILES = \
1LO-0NLO-LEP.in \
1LO-0NLO-LEP_Qt.in \
1LO-0NLO-LEP_Matchbox.in \
1LO-0NLO-LEP_FixedOrder.in \
1LO-0NLO-LHC-W.in \
1LO-0NLO-LHC-Z.in \
2LO-0NLO-LEP.in \
2LO-0NLO-LHC-W.in \
2LO-0NLO-LHC-Z.in \
2LO-1NLO-LEP.in \
2LO-1NLO-LEP_Matchbox.in \
2LO-1NLO-LHC-Z.in \
2LO-1NLO-LHC-W.in \
3LO-0NLO-LEP.in \
3LO-0NLO-LEP_only0.in \
3LO-0NLO-LEP_only1.in \
3LO-0NLO-LEP_only2.in \
3LO-0NLO-LHC-W.in \
3LO-0NLO-LHC-Z.in \
3LO-1NLO-LEP.in \
3LO-1NLO-LHC-Z.in \
3LO-1NLO-LHC-W.in \
3LO-2NLO-LEP.in \
3LO-2NLO-LEP_light.in \
3LO-2NLO-LHC-Z.in \
3LO-2NLO-LHC-W.in \
4LO-0NLO-LHC-Z.in \
1LO-0NLO-LHC-Z_OISZ.in \
2LO-0NLO-LHC-Z_OISZ.in \
3LO-0NLO-LHC-Z_OISZ.in \
1LO-0NLO-LHC-J.in \
1LO-0NLO-LHC-T.in \
1LO-0NLO-LHC-T_QtDecay.in \
1LO-0NLO-LHC-J_OISZ.in \
2LO-0NLO-LHC-J.in \
2LO-0NLO-LHC-J_only0.in \
2LO-0NLO-LHC-J_only1.in \
2LO-0NLO-LHC-T.in \
2LO-0NLO-LHC-T_QtDecay.in \
3LO-0NLO-LHC-T.in \
2LO-0NLO-LHC-J_OISZ.in 
#2LO-0NLO-LHC-J_OISZIC.in

 
dist_Matchbox_DATA = $(INPUTFILES)

CLEANFILES = done-all-links

done-all-links: $(INPUTFILES)
	@echo "Linking input files"
	@for i in $(INPUTFILES); do \
	if test -f $(srcdir)/$$i -a ! -e $$i; then \
	$(LN_S) -f $(srcdir)/$$i; fi; done
	@touch done-all-links
