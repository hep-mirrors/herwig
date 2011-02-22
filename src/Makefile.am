SUBDIRS = defaults

AUTOMAKE_OPTIONS = -Wno-portability

defaultsdir = ${pkgdatadir}/defaults

bin_PROGRAMS = Herwig++
Herwig___SOURCES = Herwig++.cc herwigopts.c herwigopts.h
Herwig___LDFLAGS = $(AM_LDFLAGS) -export-dynamic $(THEPEGLDFLAGS)
Herwig___LDADD = $(THEPEGLIB) -ldl
Herwig___CPPFLAGS = $(AM_CPPFLAGS) \
-DHERWIG_PKGDATADIR="\"$(pkgdatadir)\"" \
-DHERWIG_PKGLIBDIR="\"$(pkglibdir)\"" \
-DTHEPEG_PKGLIBDIR="\"$(THEPEGPATH)/lib/ThePEG\""
bin_SCRIPTS = herwig-config

HELPERFILES = SPhenoSPS1a.spc NMSSM.spc MUED.model MSSM.model NMSSM.model\
              RS.model ADD.model Leptoquark.model ColourRec.LEPtune

INPUTFILES = ILC.in LEP.in LHC.in TVT.in LHC-MSSM.in \
             LHC-MUED.in LHC-NMSSM.in LHC-RS.in LHC-ADD.in\
             ILC-MSSM.in ILC-MUED.in ILC-RS.in DIS.in LHC-Powheg.in \
	     TVT-Powheg.in GammaGamma.in LHC-TRP.in LHC-diffractive.in \
             LHC-GammaGamma.in LHC-LQ.in LHC-UE7-1.in

## provide -no3body for faster 'make check'
dist_noinst_DATA = LHC-MUED-no3body.in 
dist_pkgdata_DATA = $(INPUTFILES) $(HELPERFILES)
pkgdata_DATA = Makefile-UserModules

CLEANFILES = HerwigDefaults.rpo \
  *.run *.log *.out *.tex \
  multi.test *.output probs.test chisq.value \
  LHC-RS-BR.spc LHC-MUED-no3body-BR.spc LHC-MSSM-BR.spc

## checking targets ##

HerwigDefaults.rpo: Herwig++ $(srcdir)/defaults/*.in defaults/PDF.in defaults/Analysis.in $(top_builddir)/lib/*.so
	./Herwig++ init -L$(top_builddir)/lib -i defaults/HerwigDefaults.in --exitonerror

check_BSM_Full=
check_BSM=
if WANT_MSSM
check_BSM += check-LHC-MSSM 
check_BSM_Full += check-LHC-MSSM
endif
if WANT_NMSSM
check_BSM_Full += check-LHC-NMSSM
endif
if WANT_UED
check_BSM += check-LHC-MUED-no3body
check_BSM_Full += check-LHC-MUED
endif
if WANT_RS
check_BSM += check-LHC-RS
check_BSM_Full += check-LHC-RS
endif
if WANT_TRP
check_BSM_Full += check-LHC-TRP
endif
if WANT_ADD
check_BSM_Full += check-LHC-ADD
endif

check-local: check-LHC check-LEP check-DIS check-ILC check-GammaGamma $(check_BSM) check-LHC-Powheg

check-Powheg: check-LHC-Powheg check-TVT-Powheg

check-BSM: $(check_BSM_Full)

link-helper-files:
	@for i in $(HELPERFILES); do \
	if test -f $(srcdir)/$$i -a ! -e $$i; then \
	$(LN_S) -f $(srcdir)/$$i; fi; done

check-%: $(srcdir)/%.in HerwigDefaults.rpo link-helper-files
	./Herwig++ read $< --exitonerror
	./Herwig++ run  $(notdir $(subst .in,.run,$<)) -N500 -d1 --exitonerror

## valgrind targets ##

VALGRIND=valgrind --leak-check=full --num-callers=25 --freelist-vol=100000000 --leak-resolution=med --trace-children=yes

valgrind: valgrind-init valgrind-read valgrind-run

valgrind-init:
	$(VALGRIND) ./Herwig++ init -d 1 -L$(top_builddir)/lib -i defaults/HerwigDefaults.in \
&> /tmp/valgrind-init.log
valgrind-read:
	$(VALGRIND) ./Herwig++ read -d 1 LHC.in &> /tmp/valgrind-read.log
valgrind-run:
	$(VALGRIND) ./Herwig++ run -d 1 -N5 LHC.run &> /tmp/valgrind-run.log

SETUPTHEPEG=$(THEPEGPATH)/bin/setupThePEG
THEPEGREPO=$(THEPEGPATH)/lib/ThePEG/ThePEGDefaults.rpo

install-data-hook:
	@echo Creating repository
	@./Herwig++ init -L$(DESTDIR)$(pkglibdir) -i $(DESTDIR)$(defaultsdir)/HerwigDefaults.in -r $(DESTDIR)$(pkgdatadir)/HerwigDefaults.rpo --exitonerror

uninstall-hook:
	rm -f $(DESTDIR)$(pkgdatadir)/HerwigDefaults.rpo


register: register-with-thepeg-repo
register-with-thepeg-repo:
	@if test -x "$(SETUPTHEPEG)" -a -w "$(THEPEGREPO)"; \
		then echo Registering with ThePEG; \
		"$(SETUPTHEPEG)" --init \
		$(DESTDIR)$(defaultsdir)/HerwigDefaults.in \
		-r "$(THEPEGREPO)" -o "$(THEPEGREPO)" \
		-l$(DESTDIR)$(pkglibdir) --exitonerror ; \
	fi


unregister : unregister-from-thepeg-repo
unregister-from-thepeg-repo:
	@if test -x "$(SETUPTHEPEG)" -a -w "$(THEPEGREPO)"; \
		then echo Unregistering with ThePEG; \
		"$(SETUPTHEPEG)" --init defaults/HerwigCleanup.in \
		-r "$(THEPEGREPO)" -o "$(THEPEGREPO)" \
		-l$(DESTDIR)$(pkglibdir) --exitonerror ; \
	fi

EXTRA_DIST = herwigopts.ggo

herwigopts.h herwigopts.c : herwigopts.ggo
	gengetopt < $<
