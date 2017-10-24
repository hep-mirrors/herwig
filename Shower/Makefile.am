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
ShowerInteraction.h \
ShowerAlpha.h ShowerAlpha.cc ShowerAlpha.fh
