SUBDIRS = FormFactors Tau Baryon VectorMeson Perturbative \
	  WeakCurrents ScalarMeson TensorMeson Partonic

pkglib_LTLIBRARIES = libHwDecay.la Hw64Decay.la

Hw64Decay_la_LDFLAGS = -module
Hw64Decay_la_SOURCES = Hw64Decayer.h Hw64Decayer.fh Hw64Decayer.icc Hw64Decayer.cc 

libHwDecay_la_LIBADD = $(top_builddir)/PDT/libHwPDT.la
libHwDecay_la_SOURCES =  \
DecayIntegrator.cc DecayIntegrator.fh DecayIntegrator.h DecayIntegrator.icc \
DecayPhaseSpaceChannel.cc DecayPhaseSpaceChannel.fh \
DecayPhaseSpaceChannel.h DecayPhaseSpaceChannel.icc \
DecayPhaseSpaceMode.cc DecayPhaseSpaceMode.fh \
DecayPhaseSpaceMode.h DecayPhaseSpaceMode.icc \
HwDecayerBase.cc HwDecayerBase.fh HwDecayerBase.h HwDecayerBase.icc \
HwDecayHandler.cc HwDecayHandler.h HwDecayHandler.icc \
Kinematics.cc Kinematics.h Kinematics.icc \
ThreeBodyIntegrator.cc ThreeBodyIntegrator.fh \
ThreeBodyIntegrator.h ThreeBodyIntegrator.icc
