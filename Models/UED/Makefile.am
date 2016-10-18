BUILT_SOURCES  = UED__all.cc
CLEANFILES = UED__all.cc

UED__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
UEDBase.h UEDBase.fh \
UEDF1F1Z0Vertex.h \
UEDF1F1G0Vertex.h \
UEDF1F0G1Vertex.h \
UEDG1G1G0Vertex.h \
UEDG0G0G1G1Vertex.h \
UEDF1F1P0Vertex.h \
UEDF1F1W0Vertex.h \
UEDF1F0W1Vertex.h \
UEDP0H1H1Vertex.h \
UEDZ0H1H1Vertex.h \
UEDW0A1H1Vertex.h \
UEDZ0A1h1Vertex.h \
UEDW0W1W1Vertex.h \
UEDF1F0H1Vertex.h

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
UEDBase.cc \
UEDF1F1Z0Vertex.cc \
UEDF1F1G0Vertex.cc \
UEDF1F0G1Vertex.cc \
UEDG1G1G0Vertex.cc \
UEDG0G0G1G1Vertex.cc \
UEDF1F1P0Vertex.cc \
UEDF1F1W0Vertex.cc \
UEDF1F0W1Vertex.cc \
UEDP0H1H1Vertex.cc \
UEDZ0H1H1Vertex.cc \
UEDW0A1H1Vertex.cc \
UEDZ0A1h1Vertex.cc \
UEDW0W1W1Vertex.cc \
UEDF1F0H1Vertex.cc 
