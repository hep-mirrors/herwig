BUILT_SOURCES  = DalitzDecay__all.cc
CLEANFILES = DalitzDecay__all.cc

DalitzDecay__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
DalitzResonance.h\
DalitzResonance.fh\
FlatteResonance.h\
MIPWA.h\
PiPiI2.h\
DalitzKMatrix.h\
DalitzLASS.h\
DalitzGS.h\
DalitzSigma.h\
DalitzBase.h\
ScalarTo3ScalarDalitz.h\
VectorTo3PseudoScalarDalitz.h

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
DalitzResonance.cc\
FlatteResonance.cc\
MIPWA.cc\
PiPiI2.cc\
DalitzKMatrix.cc\
DalitzLASS.cc\
DalitzGS.cc\
DalitzSigma.cc\
DalitzBase.cc\
ScalarTo3ScalarDalitz.cc\
VectorTo3PseudoScalarDalitz.cc
