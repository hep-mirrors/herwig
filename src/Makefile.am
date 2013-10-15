SUBDIRS = defaults

AUTOMAKE_OPTIONS = -Wno-portability

defaultsdir = ${pkgdatadir}/defaults

bin_PROGRAMS = Herwig++
Herwig___SOURCES = Herwig++.cc herwigopts.c herwigopts.h
BUILT_SOURCES = herwigopts.c herwigopts.h
Herwig___LDFLAGS = $(AM_LDFLAGS) -export-dynamic $(THEPEGLDFLAGS)
Herwig___LDADD = $(THEPEGLIB) -ldl
Herwig___CPPFLAGS = $(AM_CPPFLAGS) \
-DHERWIG_PKGDATADIR="\"$(pkgdatadir)\"" \
-DHERWIG_PKGLIBDIR="\"$(pkglibdir)\"" \
-DTHEPEG_PKGLIBDIR="\"$(THEPEGPATH)/lib/ThePEG\""
bin_SCRIPTS = herwig-config

HELPERFILES = CMSSM40.1.1.slha RPV3.1.slha NMSSM.spc \
ADD.model \
Leptoquark.model \
LH.model \
LHTP.model \
MSSM.model \
MUED.model \
NMSSM.model \
RPV-Bi.model \
RPV-Tri.model \
RS.model \
Sextet.model \
TTBA.model \
Zprime.model \
RPV-BI.slha \
RPV-TRI.slha \
RPV-UDD.slha

INPUTFILES = \
DIS.in \
DIS-Matchbox.in \
GammaGamma.in \
ILC.in \
ILC-MSSM.in \
ILC-MUED.in \
ILC-RS.in \
LEP.in \
LEP-Matchbox.in \
LHC-ADD.in \
LHC-diffractive.in \
LHC-GammaGamma.in \
LHC.in \
LHC-Matchbox.in \
LHC-LQ.in \
LHC-MSSM.in \
LHC-MUED.in \
LHC-NMSSM.in \
LHC-Powheg.in \
LHC-RPV.in \
LHC-RS.in \
LHC-Sextet.in \
LHC-TRP.in \
LHC-TTBA.in \
LHC-MB.in \
LHC-ZP.in \
TVT.in \
TVT-Matchbox.in \
TVT-Powheg.in \
TVT-TTBA.in \
LHC-LH.in \
LHC-LHTP.in

dist_pkgdata_DATA = $(INPUTFILES) $(HELPERFILES)
pkgdata_DATA = Makefile-UserModules

CLEANFILES = HerwigDefaults.rpo \
  *.run *.log *.out *.tex \
  multi.test *.output probs.test chisq.value \
  LHC-RS-BR.spc LHC-MSSM-BR.spc LHC-RPV-BR.spc

## checking targets ##

HerwigDefaults.rpo: Herwig++ $(srcdir)/defaults/*.in defaults/PDF.in defaults/Analysis.in $(top_builddir)/lib/*.so
	./Herwig++ init -L$(top_builddir)/lib defaults/HerwigDefaults.in -D

check_BSM_Full=
check_BSM=

if WANT_BSM

check_BSM += check-LHC-RPV check-LHC-RS

check_BSM_Full += \
check-LHC-RPV check-LHC-MSSM check-ILC-MSSM \
check-LHC-NMSSM \
check-LHC-MUED check-ILC-MUED \
check-LHC-RS check-ILC-RS check-LHC-ADD \
check-LHC-LH check-LHC-LHTP \
check-LHC-TRP \
check-LHC-TTBA check-TVT-TTBA \
check-LHC-ZP \
check-LHC-LQ \
check-LHC-Sextet

endif

check-local: check-LEP check-LHC $(check_BSM) check-DIS check-ILC check-GammaGamma check-LHC-Powheg

check-Powheg: check-LHC-Powheg check-TVT-Powheg

check-BSM: $(check_BSM_Full)

check_Matchbox= 
if WANT_DIPOLE
check_Matchbox += \
check-LEP-Matchbox check-DIS-Matchbox check-LHC-Matchbox check-TVT-Matchbox
endif

check-Matchbox: $(check_Matchbox)

check-extra: check-LHC-diffractive check-LHC-GammaGamma check-LHC-MB check-TVT

check-all: check-local check-Powheg check-BSM check-Matchbox check-extra

link-helper-files:
	@for i in $(HELPERFILES); do \
	if test -f $(srcdir)/$$i -a ! -e $$i; then \
	$(LN_S) -f $(srcdir)/$$i; fi; done

CHECKCOMMAND = ./Herwig++ run  $(notdir $(subst .in,.run,$<)) -N500 -d1 -D

check-%: $(srcdir)/%.in HerwigDefaults.rpo link-helper-files
	./Herwig++ read $< -D
	@echo $(CHECKCOMMAND)
	@$(CHECKCOMMAND) && echo "# $@ OK #" \
			 || (echo "###### $@ BAD ######"; false)

## valgrind targets ##

VALGRIND=valgrind --leak-check=full --num-callers=25 --freelist-vol=100000000 --leak-resolution=med --trace-children=yes

valgrind: valgrind-init valgrind-read valgrind-run

valgrind-init:
	$(VALGRIND) ./Herwig++ init -d1 -D -L$(top_builddir)/lib defaults/HerwigDefaults.in \
&> /tmp/valgrind-init.log
valgrind-read:
	$(VALGRIND) ./Herwig++ read -d1 -D LHC.in &> /tmp/valgrind-read.log
valgrind-run:
	$(VALGRIND) ./Herwig++ run -d1 -D -N5 LHC.run &> /tmp/valgrind-run.log

SETUPTHEPEG=$(THEPEGPATH)/bin/setupThePEG
THEPEGREPO=$(THEPEGPATH)/lib/ThePEG/ThePEGDefaults.rpo

install-data-hook:
	@echo Creating repository
	@./Herwig++ init -L$(DESTDIR)$(pkglibdir) -i$(DESTDIR)$(pkgdatadir) $(DESTDIR)$(defaultsdir)/HerwigDefaults.in --repo=$(DESTDIR)$(pkgdatadir)/HerwigDefaults.rpo
	@if test -n "$(DESTDIR)"; then sed -i -e "s@$(DESTDIR)@@g" $(DESTDIR)$(pkgdatadir)/HerwigDefaults.rpo; fi

uninstall-hook:
	rm -f $(DESTDIR)$(pkgdatadir)/HerwigDefaults.rpo


register: register-with-thepeg-repo
register-with-thepeg-repo:
	@if test -x "$(SETUPTHEPEG)" -a -w "$(THEPEGREPO)"; \
		then echo Registering with ThePEG; \
		"$(SETUPTHEPEG)" --init \
		$(DESTDIR)$(defaultsdir)/HerwigDefaults.in \
		-r "$(THEPEGREPO)" -o "$(THEPEGREPO)" \
		-i $(DESTDIR)$(pkgdatadir) \
		-l$(DESTDIR)$(pkglibdir) ; \
		if test -n "$(DESTDIR)"; \
		then sed -i -e "s@$(DESTDIR)@@g" "$(THEPEGREPO)" ; fi ; \
	fi


unregister : unregister-from-thepeg-repo
unregister-from-thepeg-repo:
	@if test -x "$(SETUPTHEPEG)" -a -w "$(THEPEGREPO)"; \
		then echo Unregistering with ThePEG; \
		"$(SETUPTHEPEG)" --init defaults/HerwigCleanup.in \
		-r "$(THEPEGREPO)" -o "$(THEPEGREPO)" \
		-l$(DESTDIR)$(pkglibdir) ; \
	fi

EXTRA_DIST = herwigopts.ggo

GENGETOPT = gengetopt

%opts.h %opts.c : %opts.ggo
	$(GENGETOPT) < $<

