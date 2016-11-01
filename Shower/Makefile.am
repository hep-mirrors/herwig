SUBDIRS = Matching .
pkglib_LTLIBRARIES = HwShower.la
HwShower_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 23:0:0
HwShower_la_LIBADD = \
$(top_builddir)/PDF/libHwRemDecayer.la \
$(top_builddir)/PDF/libHwMPIPDF.la
HwShower_la_SOURCES =  \
UEBase.h UEBase.cc UEBase.fh \
Couplings/ShowerAlphaQCD.h Couplings/ShowerAlphaQCD.cc \
Couplings/ShowerAlphaQED.h Couplings/ShowerAlphaQED.cc\
ShowerHandler.h ShowerHandler.fh ShowerHandler.cc \
SplittingFunctions/HalfHalfOneSplitFn.h SplittingFunctions/HalfHalfOneSplitFn.cc\
SplittingFunctions/OneOneOneSplitFn.h SplittingFunctions/OneOneOneSplitFn.cc\
SplittingFunctions/ZeroZeroOneSplitFn.h SplittingFunctions/ZeroZeroOneSplitFn.cc\
SplittingFunctions/OneHalfHalfSplitFn.h SplittingFunctions/OneHalfHalfSplitFn.cc\
SplittingFunctions/HalfOneHalfSplitFn.h SplittingFunctions/HalfOneHalfSplitFn.cc\
Default/QTildeSudakov.cc Default/QTildeSudakov.h\
Default/QTildeModel.cc Default/QTildeModel.h\
Default/Decay_QTildeShowerKinematics1to2.cc \
Default/Decay_QTildeShowerKinematics1to2.h  \
Default/IS_QTildeShowerKinematics1to2.cc    Default/IS_QTildeShowerKinematics1to2.h  \
Default/FS_QTildeShowerKinematics1to2.cc    Default/FS_QTildeShowerKinematics1to2.h  \
Default/QTildeShowerKinematics1to2.cc Default/QTildeShowerKinematics1to2.fh\
Default/QTildeShowerKinematics1to2.h  \
Default/QTildeFinder.cc Default/QTildeFinder.h\
Default/QTildeReconstructor.cc Default/QTildeReconstructor.h Default/QTildeReconstructor.tcc \
Base/KinematicsReconstructor.cc \
Base/KinematicsReconstructor.h \
Base/KinematicsReconstructor.fh \
Base/ShowerModel.cc Base/ShowerModel.h Base/ShowerModel.fh \
Base/PartnerFinder.h Base/PartnerFinder.fh Base/PartnerFinder.cc \
Base/Evolver.h Base/Evolver.fh Base/Evolver.cc \
Base/ShowerVeto.h Base/ShowerVeto.fh Base/ShowerVeto.cc \
Base/FullShowerVeto.h Base/FullShowerVeto.fh Base/FullShowerVeto.cc

noinst_LTLIBRARIES = libHwShower.la
libHwShower_la_SOURCES = ShowerConfig.h \
Base/Branching.h \
Base/ShowerParticle.cc  Base/ShowerParticle.fh  Base/ShowerParticle.h \
Base/ShowerKinematics.fh  Base/ShowerKinematics.h Base/ShowerKinematics.cc \
Base/ShowerTree.h Base/ShowerTree.fh Base/ShowerTree.cc \
Base/ShowerProgenitor.fh Base/ShowerProgenitor.h \
Base/HardTree.h Base/HardTree.fh Base/HardTree.cc\
Base/SudakovFormFactor.cc Base/SudakovFormFactor.h Base/SudakovFormFactor.fh \
Base/HardBranching.h Base/HardBranching.fh Base/HardBranching.cc\
Couplings/ShowerAlpha.h Couplings/ShowerAlpha.cc Couplings/ShowerAlpha.fh\
SplittingFunctions/SplittingGenerator.cc SplittingFunctions/SplittingGenerator.h\
SplittingFunctions/SplittingGenerator.fh \
SplittingFunctions/SplittingFunction.h SplittingFunctions/SplittingFunction.fh \
SplittingFunctions/SplittingFunction.cc  \
Base/ShowerVertex.cc Base/ShowerVertex.fh Base/ShowerVertex.h
