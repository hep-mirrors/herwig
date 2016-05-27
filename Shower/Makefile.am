SUBDIRS = Dipole QTilde .

noinst_LTLIBRARIES = libHwShower.la
libHwShower_la_LIBADD = \
$(top_builddir)/PDF/libHwMPIPDF.la\
$(top_builddir)/PDF/libHwRemDecayer.la

libHwShower_la_SOURCES =  \
ShowerHandler.h ShowerHandler.fh ShowerHandler.cc \
ShowerVariation.h ShowerVariation.cc \
UEBase.h UEBase.cc UEBase.fh \
ShowerConfig.h ShowerConfig.cc \
QTilde/Base/Branching.h \
QTilde/Base/ShowerParticle.cc  QTilde/Base/ShowerParticle.fh  QTilde/Base/ShowerParticle.h \
QTilde/Base/ShowerKinematics.fh  QTilde/Base/ShowerKinematics.h QTilde/Base/ShowerKinematics.cc \
QTilde/Base/ShowerBasis.fh  QTilde/Base/ShowerBasis.h QTilde/Base/ShowerBasis.cc \
QTilde/Base/ShowerTree.h QTilde/Base/ShowerTree.fh QTilde/Base/ShowerTree.cc \
QTilde/Base/ShowerProgenitor.fh QTilde/Base/ShowerProgenitor.h \
QTilde/Base/HardTree.h QTilde/Base/HardTree.fh QTilde/Base/HardTree.cc\
QTilde/Base/SudakovFormFactor.cc QTilde/Base/SudakovFormFactor.h QTilde/Base/SudakovFormFactor.fh \
QTilde/Base/HardBranching.h QTilde/Base/HardBranching.fh QTilde/Base/HardBranching.cc\
QTilde/Couplings/ShowerAlpha.h QTilde/Couplings/ShowerAlpha.cc QTilde/Couplings/ShowerAlpha.fh\
QTilde/SplittingFunctions/SplittingGenerator.cc QTilde/SplittingFunctions/SplittingGenerator.h\
QTilde/SplittingFunctions/SplittingGenerator.fh \
QTilde/SplittingFunctions/SplittingFunction.h QTilde/SplittingFunctions/SplittingFunction.fh \
QTilde/SplittingFunctions/SplittingFunction.cc  \
QTilde/Base/ShowerVertex.cc QTilde/Base/ShowerVertex.fh QTilde/Base/ShowerVertex.h

