SUBDIRS = defaults

AUTOMAKE_OPTIONS = -Wno-portability

defaultsdir = ${pkgdatadir}/defaults

bin_PROGRAMS = Herwig++
Herwig___SOURCES = Herwig++.cc HerwigRun.cc HerwigRun.h
Herwig___LDFLAGS = -export-dynamic $(THEPEGLDFLAGS)
Herwig___LDADD = $(THEPEGLIB) -ldl
Herwig___CPPFLAGS = $(AM_CPPFLAGS) -DHERWIG_PKGDATADIR="\"$(pkgdatadir)\"" -DHERWIG_PKGLIBDIR="\"$(pkglibdir)\"" \
-DTHEPEG_PKGLIBDIR="\"$(THEPEGPATH)/lib/ThePEG\""
bin_SCRIPTS = herwig-config

HELPERFILES = SPhenoSPS1a.spc MUED.model MSSM.model RS.model

INPUTFILES = ILC.in LEP.in LHC.in TVT.in LHC-MSSM.in LHC-MUED.in LHC-RS.in \
             ILC-MSSM.in ILC-MUED.in ILC-RS.in DIS.in LHC-Powheg.in TVT-Powheg.in

## provide -no3body for faster 'make check'
dist_noinst_DATA = LHC-MUED-no3body.in 
dist_pkgdata_DATA = $(INPUTFILES) $(HELPERFILES)
pkgdata_DATA = Makefile-UserModules

CLEANFILES = HerwigDefaults.rpo \
  *.run *.log *.out *.tex \
  multi.test *.output probs.test chisq.value

## checking targets ##

HerwigDefaults.rpo: Herwig++ $(srcdir)/defaults/*.in defaults/PDF.in defaults/Analysis.in $(top_builddir)/lib/*.so
	./Herwig++ init -L$(top_builddir)/lib -i defaults/HerwigDefaults.in --exitonerror

check_BSM=
if WANT_MSSM
check_BSM += check-LHC-MSSM
endif
if WANT_UED
check_BSM += check-LHC-MUED-no3body
endif
if WANT_RS
check_BSM += check-LHC-RS
endif

check-local: check-LHC check-LEP check-DIS check-ILC $(check_BSM)

check-Powheg: check-LHC-Powheg check-TVT-Powheg

link-helper-files:
	@for i in $(HELPERFILES); do \
	if test -f $(srcdir)/$$i -a ! -e $$i; then \
	$(LN_S) -f $(srcdir)/$$i; fi; done

check-%: $(srcdir)/%.in HerwigDefaults.rpo link-helper-files
	./Herwig++ read $< --exitonerror
	./Herwig++ run  $(notdir $(subst .in,.run,$<)) -N100 --exitonerror

## Chisquare report ##
EXTRA_SCRIPTS = sum_chisq

chisq.value : $(srcdir)/CHISQ.in HerwigDefaults.rpo sum_chisq
	./Herwig++ read $< --exitonerror
	./Herwig++ run  $(notdir $(subst .in,.run,$<)) -N10000 --exitonerror
	@./sum_chisq $(notdir $(subst .in,.log,$<)) > $@

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
##	@if test -x "$(SETUPTHEPEG)" -a -w "$(THEPEGREPO)"; \
##		then echo Registering with ThePEG; \
##		"$(SETUPTHEPEG)" --init \
##		$(DESTDIR)$(defaultsdir)/HerwigDefaults.in \
##		-r "$(THEPEGREPO)" -o "$(THEPEGREPO)" \
##		-l$(DESTDIR)$(pkglibdir) --exitonerror ; \
##	fi

uninstall-hook:
	rm -f $(DESTDIR)$(pkgdatadir)/HerwigDefaults.rpo
##	@if test -x "$(SETUPTHEPEG)" -a -w "$(THEPEGREPO)"; \
##		then echo Unregistering with ThePEG; \
##		"$(SETUPTHEPEG)" --init defaults/HerwigCleanup.in \
##		-r "$(THEPEGREPO)" -o "$(THEPEGREPO)" \
##		-l$(DESTDIR)$(pkglibdir) --exitonerror ; \
##	fi

include $(top_srcdir)/Utilities/Makefile.am.versionstring
