
dummylibdir=/tmp

dummylib_LTLIBRARIES = LeptonTest.la

LeptonTest_la_SOURCES = Lepton/VVTest.h Lepton/VVTest.cc
LeptonTest_la_LDFLAGS = -module -version-info 1:0:0

REPO = $(top_builddir)/src/HerwigDefaults.rpo
HERWIG = $(top_builddir)/src/Herwig++
HWREAD = $(HERWIG) read -r $(REPO) -L $(builddir)/.libs
HWRUN = $(HERWIG) run

longtests : LeptonTest.la Inputs/LEP-VV.in
	$(HWREAD) Inputs/LEP-VV.in
	$(HWRUN) LEP-VV.run -N1000

clean-local:
	rm -f *.out *.log *.tex *.top *.run
