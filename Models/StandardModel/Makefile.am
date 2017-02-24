BUILT_SOURCES  = SM__all.cc
CLEANFILES = SM__all.cc

SM__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
StandardModel.fh StandardModel.h \
RunningMassBase.h \
RunningMass.h \
StandardCKM.h \
O2AlphaS.h \
AlphaEM.h \
SMFFGVertex.h \
SMFFHVertex.h \
SMFFPVertex.h \
SMFFWVertex.h \
SMFFZVertex.h \
SMGGGGVertex.h \
SMGGGVertex.h \
SMWWHVertex.h \
SMWWWVertex.h \
SMWWWWVertex.h \
SMHGGVertex.h \
SMHPPVertex.h\
SMHHHVertex.h \
SMWWHHVertex.h \
GenericSVVVertex.h \
GenericVVVVertex.h

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
StandardModel.cc\
RunningMassBase.cc \
RunningMass.cc\
StandardCKM.cc \
O2AlphaS.cc \
AlphaEM.cc  \
SMFFGVertex.cc \
SMFFHVertex.cc \
SMFFPVertex.cc \
SMFFWVertex.cc \
SMFFZVertex.cc \
SMGGGGVertex.cc \
SMGGGVertex.cc \
SMWWHVertex.cc \
SMWWWVertex.cc \
SMWWWWVertex.cc \
SMHGGVertex.cc \
SMHPPVertex.cc \
SMHHHVertex.cc \
SMWWHHVertex.cc \
GenericSVVVertex.cc \
GenericVVVVertex.cc 
