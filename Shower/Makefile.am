SUBDIRS = SplittingFunctions Default Base .
pkglib_LTLIBRARIES = HwShower.la
HwShower_la_LDFLAGS = -module -version-info 11:1:0
HwShower_la_LIBADD = Default/libHwDefaultShower.la\
$(top_builddir)/Shower/SplittingFunctions/libHwSplitting.la \
$(top_builddir)/Shower/Base/libHwShowerBase.la \
$(top_builddir)/PDF/libHwRemDecayer.la \
$(top_builddir)/PDF/libHwMPIPDF.la
HwShower_la_SOURCES =  \
UEBase.h UEBase.cc UEBase.fh\
Couplings/ShowerAlphaQCD.h Couplings/ShowerAlphaQCD.cc\
Couplings/ShowerAlphaQED.h Couplings/ShowerAlphaQED.cc\
ShowerHandler.h ShowerHandler.fh ShowerHandler.cc


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
Couplings/ShowerAlpha.h Couplings/ShowerAlpha.cc\
SplittingFunctions/SplittingGenerator.cc SplittingFunctions/SplittingGenerator.h\
SplittingFunctions/SplittingGenerator.fh \
SplittingFunctions/SplittingFunction.h SplittingFunctions/SplittingFunction.fh \
SplittingFunctions/SplittingFunction.cc 


