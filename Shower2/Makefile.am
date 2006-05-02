SUBDIRS = SplittingFunctions Couplings Kinematics MECorrections .
pkglib_LTLIBRARIES = HwNewShower.la
HwNewShower_la_LDFLAGS = -module
HwNewShower_la_LIBADD = ./SplittingFunctions/libHwNewSplitting.la \
                        ./Couplings/libHwNewCouplings.la \
                        ./Kinematics/libHwNewShowerKinematics.la \
                        ./MECorrections/libHwNewMECorrections.la
HwNewShower_la_SOURCES = ShowerConfig.h \
ShowerVariables.h ShowerVariables.fh ShowerVariables.icc ShowerVariables.cc \
ShowerHandler.h ShowerHandler.fh ShowerHandler.icc ShowerHandler.cc \
Evolver.h Evolver.fh Evolver.icc Evolver.cc \
PartnerFinder.h PartnerFinder.fh PartnerFinder.cc PartnerFinder.icc \
ShowerTree.h ShowerTree.fh ShowerTree.icc ShowerTree.cc
