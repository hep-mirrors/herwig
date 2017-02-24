SUBDIRS = Dipole QTilde

noinst_LTLIBRARIES = libHwShower.la
libHwShower_la_LIBADD = \
$(top_builddir)/PDF/libHwMPIPDF.la\
$(top_builddir)/PDF/libHwRemDecayer.la

libHwShower_la_SOURCES =  \
ShowerHandler.h ShowerHandler.fh ShowerHandler.cc \
ShowerVariation.h ShowerVariation.cc \
ShowerEventRecord.h ShowerEventRecord.cc \
PerturbativeProcess.h PerturbativeProcess.fh \
RealEmissionProcess.h RealEmissionProcess.fh \
UEBase.h UEBase.cc UEBase.fh \
Core/ShowerInteraction.h Core/ShowerConfig.h Core/ShowerConfig.cc \
Core/Base/Branching.h \
Core/Base/ShowerParticle.cc  Core/Base/ShowerParticle.fh  Core/Base/ShowerParticle.h \
Core/Base/ShowerKinematics.fh  Core/Base/ShowerKinematics.h Core/Base/ShowerKinematics.cc \
Core/Base/ShowerBasis.fh  Core/Base/ShowerBasis.h Core/Base/ShowerBasis.cc \
Core/Base/ShowerTree.h Core/Base/ShowerTree.fh Core/Base/ShowerTree.cc \
Core/Base/ShowerProgenitor.fh Core/Base/ShowerProgenitor.h \
Core/Base/HardTree.h Core/Base/HardTree.fh Core/Base/HardTree.cc\
Core/Base/SudakovFormFactor.cc Core/Base/SudakovFormFactor.h Core/Base/SudakovFormFactor.fh \
Core/Base/HardBranching.h Core/Base/HardBranching.fh Core/Base/HardBranching.cc\
Core/Couplings/ShowerAlpha.h Core/Couplings/ShowerAlpha.cc Core/Couplings/ShowerAlpha.fh\
Core/SplittingFunctions/SplittingGenerator.cc Core/SplittingFunctions/SplittingGenerator.h\
Core/SplittingFunctions/SplittingGenerator.fh \
Core/SplittingFunctions/SplittingFunction.h Core/SplittingFunctions/SplittingFunction.fh \
Core/SplittingFunctions/SplittingFunction.cc  \
Core/Base/ShowerVertex.cc Core/Base/ShowerVertex.fh Core/Base/ShowerVertex.h

