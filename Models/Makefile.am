SUBDIRS = RSModel  StandardModel General Susy UED Zprime \
          Transplanckian ADD Leptoquarks Sextet TTbAsymm \
          LH LHTP Feynrules

noinst_LTLIBRARIES = libHwStandardModel.la
nodist_libHwStandardModel_la_SOURCES = \
StandardModel/SM__all.cc

libHwStandardModel_la_CPPFLAGS = \
$(AM_CPPFLAGS) -I$(srcdir)/StandardModel

noinst_LTLIBRARIES += libHwModelGenerator.la
nodist_libHwModelGenerator_la_SOURCES = \
General/ModelGenerator__all.cc

libHwModelGenerator_la_CPPFLAGS = \
$(AM_CPPFLAGS) -I$(srcdir)/General




if WANT_BSM

pkglib_LTLIBRARIES = \
HwRSModel.la \
HwUED.la \
HwSusy.la \
HwNMSSM.la \
HwRPV.la \
HwZprimeModel.la \
HwTransplanck.la \
HwADDModel.la \
HwLeptoquarkModel.la \
HwSextetModel.la \
HwTTbAModel.la \
HwLHModel.la \
HwLHTPModel.la

endif

#############

HwRSModel_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 11:0:0

HwRSModel_la_CPPFLAGS = \
$(AM_CPPFLAGS) -I$(srcdir)/RSModel

nodist_HwRSModel_la_SOURCES = \
RSModel/RS__all.cc

#############

HwUED_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 11:0:0

HwUED_la_CPPFLAGS = \
$(AM_CPPFLAGS) -I$(srcdir)/UED

nodist_HwUED_la_SOURCES = \
UED/UED__all.cc

#############

HwSusy_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 14:0:0

HwSusy_la_CPPFLAGS = \
$(AM_CPPFLAGS) -I$(srcdir)/Susy

nodist_HwSusy_la_SOURCES = \
Susy/Susy__all.cc

#############

HwNMSSM_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 6:0:0

HwNMSSM_la_CPPFLAGS = \
$(AM_CPPFLAGS) -I$(srcdir)/Susy/NMSSM

nodist_HwNMSSM_la_SOURCES = \
Susy/NMSSM/NMSSM__all.cc

#############

HwRPV_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 4:0:0

HwRPV_la_CPPFLAGS = \
$(AM_CPPFLAGS) -I$(srcdir)/Susy/RPV

nodist_HwRPV_la_SOURCES = \
Susy/RPV/RPV__all.cc

#############

HwZprimeModel_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 4:0:0

HwZprimeModel_la_SOURCES = \
Zprime/ZprimeModel.cc  Zprime/ZprimeModelZPQQVertex.cc

#############

HwTransplanck_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 5:0:0

HwTransplanck_la_SOURCES = \
Transplanckian/METRP2to2.cc

#############

HwADDModel_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 5:0:0

HwADDModel_la_CPPFLAGS = \
$(AM_CPPFLAGS) -I$(srcdir)/ADD

nodist_HwADDModel_la_SOURCES = \
ADD/ADD__all.cc

#############

HwLeptoquarkModel_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 6:0:0

HwLeptoquarkModel_la_CPPFLAGS = \
$(AM_CPPFLAGS) -I$(srcdir)/Leptoquarks

nodist_HwLeptoquarkModel_la_SOURCES = \
Leptoquarks/Leptoquark__all.cc

#############

HwSextetModel_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 4:0:0

HwSextetModel_la_CPPFLAGS = \
$(AM_CPPFLAGS) -I$(srcdir)/Sextet

nodist_HwSextetModel_la_SOURCES = \
Sextet/Sextet__all.cc

#############

HwTTbAModel_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 4:0:0

HwTTbAModel_la_CPPFLAGS = \
$(AM_CPPFLAGS) -I$(srcdir)/TTbAsymm

nodist_HwTTbAModel_la_SOURCES = \
TTbAsymm/TTbA__all.cc

#############

HwLHModel_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 6:0:0

HwLHModel_la_CPPFLAGS = \
$(AM_CPPFLAGS) -I$(srcdir)/LH

nodist_HwLHModel_la_SOURCES = \
LH/LH__all.cc

#############

HwLHTPModel_la_LDFLAGS = \
$(AM_LDFLAGS) -module -version-info 6:0:0

HwLHTPModel_la_LIBADD = \
$(GSLLIBS)

HwLHTPModel_la_CPPFLAGS = \
$(AM_CPPFLAGS) $(GSLINCLUDE) -I$(srcdir)/LHTP

nodist_HwLHTPModel_la_SOURCES = \
LHTP/LHTP__all.cc

#############
