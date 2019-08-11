noinst_LTLIBRARIES = libHwDipoleShowerKernels.la

nodist_libHwDipoleShowerKernels_la_SOURCES = \
DipoleKernels__all.cc

BUILT_SOURCES  = DipoleKernels__all.cc
CLEANFILES = DipoleKernels__all.cc

DipoleKernels__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
	DipoleSplittingKernel.h \
	ColourMatrixElementCorrection.h \
	FFMqx2qgxDipoleKernel.h \
	FFMgx2ggxDipoleKernel.h \
	FFMgx2qqxDipoleKernel.h \
	FFqx2qgxDipoleKernel.h  \
	FFgx2ggxDipoleKernel.h  \
	FFgx2qqxDipoleKernel.h  \
	FIqx2qgxDipoleKernel.h  \
	FIgx2ggxDipoleKernel.h  \
	FIgx2qqxDipoleKernel.h  \
	IFqx2qgxDipoleKernel.h  \
	IFqx2gqxDipoleKernel.h  \
	IFgx2ggxDipoleKernel.h  \
	IFgx2qqxDipoleKernel.h  \
	IIqx2qgxDipoleKernel.h  \
	IIqx2gqxDipoleKernel.h  \
	IIgx2ggxDipoleKernel.h  \
	IIgx2qqxDipoleKernel.h  \
	FIMqx2qgxDipoleKernel.h \
	FIMgx2qqxDipoleKernel.h \
	IFMqx2qgxDipoleKernel.h \
	IFMqx2gqxDipoleKernel.h \
	IFMgx2ggxDipoleKernel.h \
	IFMgx2qqxDipoleKernel.h \
	FIMDecaygx2qqxDipoleKernel.h \
	FIMDecayqx2qgxDipoleKernel.h \
	FIMDecaygx2ggxDipoleKernel.h

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
	DipoleSplittingKernel.cc \
	ColourMatrixElementCorrection.cc \
	FFMqx2qgxDipoleKernel.cc \
	FFMgx2ggxDipoleKernel.cc \
	FFMgx2qqxDipoleKernel.cc \
	FFqx2qgxDipoleKernel.cc  \
	FFgx2ggxDipoleKernel.cc  \
	FFgx2qqxDipoleKernel.cc  \
	FIqx2qgxDipoleKernel.cc  \
	FIgx2ggxDipoleKernel.cc  \
	FIgx2qqxDipoleKernel.cc  \
	IFqx2qgxDipoleKernel.cc  \
	IFqx2gqxDipoleKernel.cc  \
	IFgx2ggxDipoleKernel.cc  \
	IFgx2qqxDipoleKernel.cc  \
	IIqx2qgxDipoleKernel.cc  \
	IIqx2gqxDipoleKernel.cc  \
	IIgx2ggxDipoleKernel.cc  \
	IIgx2qqxDipoleKernel.cc  \
	FIMqx2qgxDipoleKernel.cc \
	FIMgx2qqxDipoleKernel.cc \
	IFMqx2qgxDipoleKernel.cc \
	IFMqx2gqxDipoleKernel.cc \
	IFMgx2ggxDipoleKernel.cc \
	IFMgx2qqxDipoleKernel.cc \
	FIMDecaygx2qqxDipoleKernel.cc \
	FIMDecayqx2qgxDipoleKernel.cc \
	FIMDecaygx2ggxDipoleKernel.cc

