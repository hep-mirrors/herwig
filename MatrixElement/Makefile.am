SUBDIRS = General Lepton Hadron DIS Powheg Gamma Matchbox Reweighters

if WANT_LIBFASTJET
SUBDIRS += FxFx
endif

noinst_LTLIBRARIES = libHwME.la 

libHwME_la_SOURCES =  \
HwMEBase.h HwMEBase.fh HwMEBase.cc \
MEfftoVH.h MEfftoVH.cc \
MEfftoffH.h MEfftoffH.cc \
HardVertex.fh HardVertex.h HardVertex.cc \
ProductionMatrixElement.h ProductionMatrixElement.cc \
DrellYanBase.h DrellYanBase.cc \
BlobME.h BlobME.cc \
MEMinBias.h MEMinBias.cc
