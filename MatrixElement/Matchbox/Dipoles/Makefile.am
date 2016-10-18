noinst_LTLIBRARIES = libHwMatchboxDipoles.la

nodist_libHwMatchboxDipoles_la_SOURCES = \
Dipole__all.cc

BUILT_SOURCES  = Dipole__all.cc
CLEANFILES = Dipole__all.cc

Dipole__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
FFggxDipole.h \
FFqgxDipole.h \
FFqqxDipole.h \
FIggxDipole.h \
FIqgxDipole.h \
FIqqxDipole.h \
IFggxDipole.h \
IFgqxDipole.h \
IFqgxDipole.h \
IFqqxDipole.h \
IIggxDipole.h \
IIgqxDipole.h \
IIqgxDipole.h \
IIqqxDipole.h \
FFMggxDipole.h \
FFMqgxDipole.h \
FFMsqgxDipole.h \
FFMqqxDipole.h \
FIMqgxDipole.h \
FIMsqgxDipole.h \
FIMqqxDipole.h \
IFMggxDipole.h \
IFMgqxDipole.h \
IFMqgxDipole.h \
IFMqqxDipole.h \
SubtractionDipole.fh \
SubtractionDipole.h


DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
FFggxDipole.cc \
FFqgxDipole.cc \
FFqqxDipole.cc \
FIggxDipole.cc \
FIqgxDipole.cc \
FIqqxDipole.cc \
IFggxDipole.cc \
IFgqxDipole.cc \
IFqgxDipole.cc \
IFqqxDipole.cc \
IIggxDipole.cc \
IIgqxDipole.cc \
IIqgxDipole.cc \
IIqqxDipole.cc \
FFMggxDipole.cc \
FFMqgxDipole.cc \
FFMsqgxDipole.cc \
FFMqqxDipole.cc \
FIMqgxDipole.cc \
FIMsqgxDipole.cc \
FIMqqxDipole.cc \
IFMggxDipole.cc \
IFMgqxDipole.cc \
IFMqgxDipole.cc \
IFMqqxDipole.cc \
SubtractionDipole.cc
