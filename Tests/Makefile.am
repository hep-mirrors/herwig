AUTOMAKE_OPTIONS = -Wno-portability

dummylibdir=/tmp

dummylib_LTLIBRARIES = LeptonTest.la

LeptonTest_la_LDFLAGS = -module -version-info 1:0:0
LeptonTest_la_SOURCES = \
Lepton/VVTest.h Lepton/VVTest.cc \
Lepton/VHTest.h Lepton/VHTest.cc

REPO = $(top_builddir)/src/HerwigDefaults.rpo
HERWIG = $(top_builddir)/src/Herwig++
HWREAD = $(HERWIG) read -r $(REPO) -L $(builddir)/.libs
HWRUN = $(HERWIG) run

tests : test-LEP-VV test-LEP-VH

test-% : Inputs/%.in
	$(HWREAD) $<
	$(HWRUN) $(notdir $(subst .in,.run,$<)) -N1000

test-LEP-% : LeptonTest.la

clean-local:
	rm -f *.out *.log *.tex *.top *.run
