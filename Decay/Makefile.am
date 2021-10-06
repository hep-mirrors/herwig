SUBDIRS = FormFactors Tau Baryon VectorMeson HeavyMeson Perturbative \
	  WeakCurrents ScalarMeson TensorMeson Partonic General Radiation

if HAVE_EVTGEN
SUBDIRS += EvtGen
endif

noinst_LTLIBRARIES = libHwDecay.la

libHwDecay_la_LIBADD = \
$(top_builddir)/PDT/libHwPDT.la

nodist_libHwDecay_la_SOURCES =  \
hwdecay__all.cc

BUILT_SOURCES  = hwdecay__all.cc
CLEANFILES = hwdecay__all.cc

hwdecay__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
DecayIntegrator.fh DecayIntegrator.h \
PhaseSpaceMode.fh PhaseSpaceMode.h \
HwDecayerBase.fh HwDecayerBase.h \
HwDecayHandler.h \
DecayVertex.fh DecayVertex.h \
DecayMatrixElement.fh DecayMatrixElement.h \
TwoBodyDecayMatrixElement.h \
GeneralDecayMatrixElement.fh GeneralDecayMatrixElement.h \
BranchingRatioReweighter.h\
PerturbativeDecayer.h \
ResonanceHelpers.h PhaseSpaceChannel.h IsoSpin.h

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
DecayIntegrator.cc \
PhaseSpaceChannel.cc \
PhaseSpaceMode.cc \
HwDecayerBase.cc \
HwDecayHandler.cc \
DecayVertex.cc \
DecayMatrixElement.cc \
TwoBodyDecayMatrixElement.cc \
GeneralDecayMatrixElement.cc \
BranchingRatioReweighter.cc\
PerturbativeDecayer.cc

##################

pkglib_LTLIBRARIES = Hw64Decay.la

Hw64Decay_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 12:0:0

Hw64Decay_la_SOURCES = \
Hw64Decayer.h Hw64Decayer.cc

##################

pkglib_LTLIBRARIES += HwMamboDecay.la

HwMamboDecay_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 13:0:0

HwMamboDecay_la_SOURCES = \
MamboDecayer.h MamboDecayer.cc

##################

noinst_LTLIBRARIES += libHwFormFactor.la

libHwFormFactor_la_SOURCES = \
FormFactors/BaryonFormFactor.cc \
FormFactors/ScalarFormFactor.cc \
FormFactors/BtoSGammaHadronicMass.cc \
FormFactors/KMatrix.cc \
FormFactors/ScalarAmplitude.cc \
FormFactors/BaryonFormFactor.fh \
FormFactors/BaryonFormFactor.h \
FormFactors/ScalarFormFactor.fh \
FormFactors/ScalarFormFactor.h \
FormFactors/BtoSGammaHadronicMass.h \
FormFactors/BtoSGammaHadronicMass.fh \
FormFactors/KMatrix.h \
FormFactors/KMatrix.fh \
FormFactors/ScalarAmplitude.h \
FormFactors/ScalarAmplitude.fh

pkglib_LTLIBRARIES += HwFormFactors.la

HwFormFactors_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 13:0:0

HwFormFactors_la_CPPFLAGS = \
$(AM_CPPFLAGS) -I$(srcdir)/FormFactors

nodist_HwFormFactors_la_SOURCES = \
FormFactors/Formfactor__all.cc

##################

pkglib_LTLIBRARIES += HwTauDecay.la

HwTauDecay_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 13:0:0

HwTauDecay_la_SOURCES = \
Tau/TauDecayer.cc

##################

pkglib_LTLIBRARIES += HwBaryonDecay.la

HwBaryonDecay_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 12:0:0

HwBaryonDecay_la_LIBADD = \
$(top_builddir)/PDT/libHwBaryonWidth.la

HwBaryonDecay_la_CPPFLAGS = \
$(AM_CPPFLAGS) -I$(srcdir)/Baryon

nodist_HwBaryonDecay_la_SOURCES = \
Baryon/BaryonDecayer__all.cc

##################

pkglib_LTLIBRARIES += HwVMDecay.la

HwVMDecay_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 12:0:0

HwVMDecay_la_CPPFLAGS = \
$(AM_CPPFLAGS) -I$(srcdir)/VectorMeson

nodist_HwVMDecay_la_SOURCES = \
VectorMeson/VMDecayer__all.cc

##################

pkglib_LTLIBRARIES += HwPerturbativeDecay.la 

HwPerturbativeDecay_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 13:0:0

HwPerturbativeDecay_la_CPPFLAGS = \
$(AM_CPPFLAGS) -I$(srcdir)/Perturbative

nodist_HwPerturbativeDecay_la_SOURCES = \
Perturbative/Perturbative__all.cc

##################

noinst_LTLIBRARIES += libHwWeakCurrent.la
libHwWeakCurrent_la_SOURCES = \
WeakCurrents/WeakCurrent.cc \
WeakCurrents/LeptonNeutrinoCurrent.cc \
WeakCurrents/WeakCurrent.fh \
WeakCurrents/WeakCurrent.h \
WeakCurrents/LeptonNeutrinoCurrent.fh \
WeakCurrents/LeptonNeutrinoCurrent.h \
WeakCurrents/WeakBaryonCurrent.cc \
WeakCurrents/WeakBaryonCurrent.h

pkglib_LTLIBRARIES += HwWeakCurrents.la

HwWeakCurrents_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 13:0:0

HwWeakCurrents_la_CPPFLAGS = \
$(AM_CPPFLAGS) -I$(srcdir)/WeakCurrents

nodist_HwWeakCurrents_la_SOURCES = \
WeakCurrents/WeakCurrents__all.cc

##################

pkglib_LTLIBRARIES += HwHMDecay.la

HwHMDecay_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 11:0:0

HwHMDecay_la_CPPFLAGS = \
$(AM_CPPFLAGS) -I$(srcdir)/HeavyMeson

nodist_HwHMDecay_la_SOURCES = \
HeavyMeson/HMDecayer__all.cc

##################

pkglib_LTLIBRARIES += HwSMDecay.la

HwSMDecay_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 14:0:0

HwSMDecay_la_CPPFLAGS = \
$(AM_CPPFLAGS) -I$(srcdir)/ScalarMeson

nodist_HwSMDecay_la_SOURCES = \
ScalarMeson/SMDecayer__all.cc

##################

pkglib_LTLIBRARIES += HwTMDecay.la

HwTMDecay_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 12:0:0

HwTMDecay_la_CPPFLAGS = \
$(AM_CPPFLAGS) -I$(srcdir)/TensorMeson

nodist_HwTMDecay_la_SOURCES = \
TensorMeson/TMDecayer__all.cc

##################

pkglib_LTLIBRARIES += HwPartonicDecay.la

HwPartonicDecay_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 13:0:0

HwPartonicDecay_la_CPPFLAGS = \
$(AM_CPPFLAGS) -I$(srcdir)/Partonic

nodist_HwPartonicDecay_la_SOURCES = \
Partonic/Partonic__all.cc

##################

noinst_LTLIBRARIES += libHwDecRad.la

libHwDecRad_la_SOURCES = \
Radiation/DecayRadiationGenerator.cc \
Radiation/QEDRadiationHandler.cc \
Radiation/DecayRadiationGenerator.h \
Radiation/DecayRadiationGenerator.fh \
Radiation/QEDRadiationHandler.fh \
Radiation/QEDRadiationHandler.h


pkglib_LTLIBRARIES += HwSOPHTY.la

HwSOPHTY_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 7:0:0

HwSOPHTY_la_CPPFLAGS = \
$(AM_CPPFLAGS) -I$(srcdir)/Radiation

nodist_HwSOPHTY_la_SOURCES = \
Radiation/Sophty__all.cc
##################
