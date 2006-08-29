SUBDIRS = SplittingFunctions Couplings Kinematics MECorrections Nason .
pkglib_LTLIBRARIES = HwShower.la
HwShower_la_LDFLAGS = -module
HwShower_la_LIBADD = ./SplittingFunctions/libHwSplitting.la \
                        ./Couplings/libHwCouplings.la \
                        ./Kinematics/libHwShowerKinematics.la \
                        ./MECorrections/libHwMECorrections.la
HwShower_la_SOURCES = ShowerConfig.h \
ShowerVariables.h ShowerVariables.fh ShowerVariables.icc ShowerVariables.cc \
ShowerHandler.h ShowerHandler.fh ShowerHandler.icc ShowerHandler.cc \
Evolver.h Evolver.fh Evolver.icc Evolver.cc \
PartnerFinder.h PartnerFinder.fh PartnerFinder.cc PartnerFinder.icc \
ShowerTree.h ShowerTree.fh ShowerTree.icc ShowerTree.cc \
ShowerProgenitor.fh ShowerProgenitor.h ShowerProgenitor.icc
