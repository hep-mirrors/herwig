SUBDIRS = defaults
defaultsdir = ${pkgdatadir}/defaults

bin_PROGRAMS = Herwig++
Herwig___SOURCES = Herwig++.cc
Herwig___LDFLAGS = -export-dynamic $(THEPEGLDFLAGS)
Herwig___LDADD = $(top_builddir)/lib/libHerwig.la $(THEPEGLIB) -ldl

bin_SCRIPTS = herwig-config

HELPERFILES = SPhenoSPS1a.spc MUED.model MSSM.model RS.model

INPUTFILES = ILC.in LEP.in LHC.in TVT.in LHC-MSSM.in LHC-MUED.in LHC-RS.in \
             ILC-MSSM.in ILC-MUED.in ILC-RS.in DIS.in

dist_pkgdata_DATA = $(INPUTFILES) $(HELPERFILES)
pkgdata_DATA = Makefile-UserModules

########## SVN version string. Only use in cc files! ####################
EXTRA_DIST = versionstring.h.in
BUILT_SOURCES = versionstring.h
if USE_SVNVERSION # # # # # #
version:
	@svnversion -n $(top_srcdir) > version.tmp.new

version.tmp: version
	@if ! diff -q version.tmp{.new,}; then mv -f version.tmp{.new,}; fi

versionstring.h: version.tmp versionstring.h.in
	@echo "creating versionstring.h"
	@sed s%@HERWIGVERSIONSTRING@%"$(PACKAGE_STRING) `cd $(top_srcdir); svn info | grep URL | sed s@.*herwig/@@` r`cat version.tmp`"% \
	  $(srcdir)/versionstring.h.in > versionstring.h
else # # # # # # # # # # # # #
versionstring.h:
	@echo "creating versionstring.h"
	@sed s%@HERWIGVERSIONSTRING@%"$(PACKAGE_STRING)"% \
	  $(srcdir)/versionstring.h.in > versionstring.h
	@touch versionstring.h
endif # # # # # # # # # # # # # 
########################################################################

CLEANFILES = HerwigDefaults.rpo \
  *.run *.log *.out *.tex \
  versionstring.h version.tmp \
  version.tmp.new multi.test \
  *.output

## checking targets ##

HerwigDefaults.rpo: Herwig++ $(srcdir)/defaults/*.in defaults/PDF.in defaults/Analysis.in $(top_builddir)/lib/*.so
	./Herwig++ init -L$(top_builddir)/lib -i defaults/HerwigDefaults.in --exitonerror

check-local: check-LHC check-LEP check-TVT check-ILC check-LHC-MSSM \
	check-LHC-MUED check-LHC-RS

link-helper-files:
	@for i in $(HELPERFILES); do \
	if test -f $(srcdir)/$$i -a ! -e $$i; then \
	$(LN_S) -f $(srcdir)/$$i; fi; done

check-%: $(srcdir)/%.in HerwigDefaults.rpo link-helper-files
	./Herwig++ read $< --exitonerror
	./Herwig++ run  $(notdir $(subst .in,.run,$<)) -N100 --exitonerror

## valgrind targets ##

VALGRIND=valgrind --leak-check=full --num-callers=25 --freelist-vol=100000000 --leak-resolution=med --trace-children=yes

valgrind: valgrind-init valgrind-read valgrind-run

valgrind-init:
	$(VALGRIND) ./Herwig++ init -L$(top_builddir)/lib -i defaults/HerwigDefaults.in \
&> /tmp/valgrind-init.log
valgrind-read:
	$(VALGRIND) ./Herwig++ read LHC.in &> /tmp/valgrind-read.log
valgrind-run:
	$(VALGRIND) ./Herwig++ run -N1 LHC.run &> /tmp/valgrind-run.log

install-data-hook:
	./Herwig++ init -L$(DESTDIR)$(pkglibdir) -i $(DESTDIR)$(defaultsdir)/HerwigDefaults.in -r $(DESTDIR)$(pkgdatadir)/HerwigDefaults.rpo --exitonerror

uninstall-hook:
	rm -f $(DESTDIR)$(pkgdatadir)/HerwigDefaults.rpo
