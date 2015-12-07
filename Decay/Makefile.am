SUBDIRS = FormFactors Tau Baryon VectorMeson Perturbative \
	  WeakCurrents ScalarMeson TensorMeson Partonic General Radiation

noinst_LTLIBRARIES = libHwDecay.la
pkglib_LTLIBRARIES = Hw64Decay.la HwMamboDecay.la

Hw64Decay_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 8:0:0
Hw64Decay_la_SOURCES = Hw64Decayer.h Hw64Decayer.cc

HwMamboDecay_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 9:0:0
HwMamboDecay_la_SOURCES = MamboDecayer.h MamboDecayer.cc

libHwDecay_la_LIBADD = $(top_builddir)/PDT/libHwPDT.la
libHwDecay_la_SOURCES =  \
DecayIntegrator.cc DecayIntegrator.fh DecayIntegrator.h \
DecayPhaseSpaceChannel.cc DecayPhaseSpaceChannel.fh \
DecayPhaseSpaceChannel.h \
DecayPhaseSpaceMode.cc DecayPhaseSpaceMode.fh \
DecayPhaseSpaceMode.h \
HwDecayerBase.cc HwDecayerBase.fh HwDecayerBase.h \
HwDecayHandler.cc HwDecayHandler.h\
DecayVertex.fh DecayVertex.h DecayVertex.cc \
DecayMatrixElement.h DecayMatrixElement.cc DecayMatrixElement.fh \
TwoBodyDecayMatrixElement.h TwoBodyDecayMatrixElement.cc \
GeneralDecayMatrixElement.h GeneralDecayMatrixElement.cc GeneralDecayMatrixElement.fh \
BranchingRatioReweighter.h BranchingRatioReweighter.cc
