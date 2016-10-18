noinst_LTLIBRARIES = libHwGeneralDecay.la
libHwGeneralDecay_la_SOURCES = \
GeneralTwoBodyDecayer.h GeneralTwoBodyDecayer.fh GeneralTwoBodyDecayer.cc \
GeneralThreeBodyDecayer.h GeneralThreeBodyDecayer.fh \
GeneralThreeBodyDecayer.cc \
GeneralCurrentDecayer.h GeneralCurrentDecayer.fh GeneralCurrentDecayer.cc \
GeneralFourBodyDecayer.h GeneralFourBodyDecayer.fh \
GeneralFourBodyDecayer.cc \
StoFFFFDecayer.h StoFFFFDecayer.cc

nodist_libHwGeneralDecay_la_SOURCES = \
GeneralDecayer__all.cc 

BUILT_SOURCES  = GeneralDecayer__all.cc
CLEANFILES = GeneralDecayer__all.cc

GeneralDecayer__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
FFSDecayer.h \
FFVDecayer.h \
SFFDecayer.h \
SSSDecayer.h \
SSVDecayer.h \
SVVDecayer.h \
TFFDecayer.h \
TSSDecayer.h \
TVVDecayer.h \
VFFDecayer.h \
VSSDecayer.h \
VVSDecayer.h \
VVVDecayer.h \
SRFDecayer.h \
FRSDecayer.h \
FRVDecayer.h \
FtoFFFDecayer.h \
StoSFFDecayer.h \
StoFFVDecayer.h \
VtoFFVDecayer.h \
FtoFVVDecayer.h \
FFVCurrentDecayer.h

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
FFSDecayer.cc \
FFVDecayer.cc \
SFFDecayer.cc \
SSSDecayer.cc \
SSVDecayer.cc \
SVVDecayer.cc \
TFFDecayer.cc \
TSSDecayer.cc \
TVVDecayer.cc \
VFFDecayer.cc \
VSSDecayer.cc \
VVSDecayer.cc \
VVVDecayer.cc \
SRFDecayer.cc \
FRSDecayer.cc \
FRVDecayer.cc \
FtoFFFDecayer.cc \
StoSFFDecayer.cc \
StoFFVDecayer.cc \
VtoFFVDecayer.cc \
FtoFVVDecayer.cc \
FFVCurrentDecayer.cc
