SUBDIRS =

pkglib_LTLIBRARIES = HwHEJ.la
HwHEJ_la_LDFLAGS = -module -version-info 11:1:0

HwHEJ_la_LIBADD = \
-L$(HEJLIBPATH) -lJets

HwHEJ_la_CPPFLAGS = $(AM_CPPFLAGS)
HwHEJ_la_CPPFLAGS += $(HEJCPPFLAGS)

HwHEJ_la_SOURCES = \
HEJPhasespace.h HEJPhasespace.cc \
HEJMEBase.h HEJMEBase.cc \
HEJMEPP2Jets.h HEJMEPP2Jets.cc \
HEJFactory.h HEJFactory.cc
