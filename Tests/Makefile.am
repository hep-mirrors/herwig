AUTOMAKE_OPTIONS = -Wno-portability
AM_LDFLAGS += -module -avoid-version -rpath /dummy/path/not/used



EXTRA_LTLIBRARIES = LeptonTest.la

LeptonTest_la_SOURCES = \
Lepton/VVTest.h Lepton/VVTest.cc \
Lepton/VBFTest.h Lepton/VBFTest.cc \
Lepton/VHTest.h Lepton/VHTest.cc





REPO = $(top_builddir)/src/HerwigDefaults.rpo
HERWIG = $(top_builddir)/src/Herwig++
HWREAD = $(HERWIG) read -r $(REPO) -L $(builddir)/.libs
HWRUN = $(HERWIG) run

tests : tests-LEP

tests-LEP : test-LEP-VV test-LEP-VH test-LEP-VBF test-LEP-default test-LEP-BB

test-LEP-% : Inputs/LEP-%.in LeptonTest.la
	$(HWREAD) $<
	$(HWRUN) $(notdir $(subst .in,.run,$<)) -N 100000

clean-local:
	rm -f *.out *.log *.tex *.top *.run *.dump *.mult *.Bmult
