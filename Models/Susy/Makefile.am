SUBDIRS = NMSSM RPV

BUILT_SOURCES  = Susy__all.cc
CLEANFILES = Susy__all.cc

Susy__all.cc : $(DIR_H_FILES) $(DIR_CC_FILES) Makefile
	@echo "Concatenating .cc files into $@"
	@$(top_srcdir)/cat_with_cpplines $(DIR_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

DIR_H_FILES = $(addprefix $(srcdir)/,$(ALL_H_FILES))
ALL_H_FILES = \
SusyBase.h SusyBase.fh \
MSSM.h MSSM.fh\
MixingMatrix.h MixingMatrix.fh \
SSCFSVertex.h    \
SSGFSVertex.h    \
SSHSFSFVertex.h  \
SSNFSVertex.h    \
SSWSSVertex.h    \
SSWWSSVertex.h    \
SSWGSSVertex.h   \
SSGSSVertex.h    \
SSGGSQSQVertex.h \
SSGSGSGVertex.h  \
SSNNZVertex.h    \
SSCCZVertex.h    \
SSCNWVertex.h    \
SSFFHVertex.h    \
SSGOGOHVertex.h  \
SSWWHVertex.h    \
SSWWHHVertex.h   \
SSWHHVertex.h    \
SSHHHVertex.h    \
SSHGGVertex.h    \
SSHPPVertex.h    \
SSNNPVertex.h \
SSGNGVertex.h \
SSNCTVertex.h \
SSGVNHVertex.h \
SSGVNVVertex.h \
SSGVFSVertex.h 

DIR_CC_FILES = $(addprefix $(srcdir)/,$(ALL_CC_FILES))
ALL_CC_FILES = \
SusyBase.cc  \
MSSM.cc  \
MixingMatrix.cc\
SSCFSVertex.cc    \
SSGFSVertex.cc    \
SSHSFSFVertex.cc  \
SSNFSVertex.cc    \
SSWSSVertex.cc    \
SSWWSSVertex.cc    \
SSWGSSVertex.cc   \
SSGSSVertex.cc    \
SSGGSQSQVertex.cc \
SSGSGSGVertex.cc  \
SSNNZVertex.cc    \
SSCCZVertex.cc    \
SSCNWVertex.cc    \
SSFFHVertex.cc    \
SSGOGOHVertex.cc  \
SSWWHVertex.cc    \
SSWWHHVertex.cc   \
SSWHHVertex.cc    \
SSHHHVertex.cc    \
SSHGGVertex.cc    \
SSHPPVertex.cc    \
SSNNPVertex.cc \
SSGNGVertex.cc \
SSNCTVertex.cc \
SSGVNHVertex.cc\
SSGVNVVertex.cc\
SSGVFSVertex.cc
