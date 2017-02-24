BUILT_SOURCES  = ModelGenerator__all.cc
CLEANFILES = ModelGenerator__all.cc

ModelGenerator__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
ModelGenerator.h ModelGenerator.fh \
DecayConstructor.h DecayConstructor.fh  \
NBodyDecayConstructorBase.h \
NBodyDecayConstructorBase.fh \
TwoBodyDecayConstructor.h \
TwoToTwoProcessConstructor.h TwoToTwoProcessConstructor.fh \
HardProcessConstructor.h HardProcessConstructor.fh \
HiggsVectorBosonProcessConstructor.h \
HiggsVBFProcessConstructor.h \
QQHiggsProcessConstructor.h \
HPDiagram.h \
ThreeBodyDecayConstructor.h \
FourBodyDecayConstructor.h  \
WeakCurrentDecayConstructor.h \
ResonantProcessConstructor.h \
ResonantProcessConstructor.fh \
VVSLoopVertex.fh VVSLoopVertex.h \
GenericHGGVertex.h \
GenericHPPVertex.h \
BSMWidthGenerator.h BSMWidthGenerator.fh \
PrototypeVertex.h \
TBDiagram.h TwoBodyDecay.h\
BSMModel.h 

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
ModelGenerator.cc\
DecayConstructor.cc \
NBodyDecayConstructorBase.cc \
TwoBodyDecayConstructor.cc \
TwoToTwoProcessConstructor.cc \
HardProcessConstructor.cc \
HiggsVectorBosonProcessConstructor.cc \
HiggsVBFProcessConstructor.cc \
QQHiggsProcessConstructor.cc \
ThreeBodyDecayConstructor.cc \
FourBodyDecayConstructor.cc \
WeakCurrentDecayConstructor.cc \
ResonantProcessConstructor.cc \
VVSLoopVertex.cc \
GenericHGGVertex.cc \
GenericHPPVertex.cc \
BSMWidthGenerator.cc \
PrototypeVertex.cc \
BSMModel.cc
