BUILT_SOURCES  = Perturbative__all.cc
CLEANFILES = Perturbative__all.cc

Perturbative__all.cc : $(ALL_CC_FILES) Makefile
	$(top_srcdir)/cat_with_cpplines $(ALL_CC_FILES) > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES) Makefile

ALL_H_FILES = \
 SMWDecayer.h\
 SMZDecayer.h\
 SMTopDecayer.h\
 SMTopPOWHEGDecayer.h\
 SMZFermionsPOWHEGDecayer.h\
 SMWFermionsPOWHEGDecayer.h


ALL_CC_FILES = \
 SMWDecayer.cc  \
 SMZDecayer.cc  \
 SMTopDecayer.cc \
 SMTopPOWHEGDecayer.cc \
 SMZFermionsPOWHEGDecayer.cc  \
 SMWFermionsPOWHEGDecayer.cc  





BUILT_SOURCES  += PerturbativeHiggs__all.cc
CLEANFILES += PerturbativeHiggs__all.cc

PerturbativeHiggs__all.cc : $(ALL_CC_FILES_2)
	$(top_srcdir)/cat_with_cpplines $(ALL_CC_FILES) > $@

EXTRA_DIST += $(ALL_H_FILES_2) $(ALL_CC_FILES_2)

ALL_H_FILES_2 = \
 SMHiggsFermionsDecayer.h \
 SMHiggsFermionsPOWHEGDecayer.h \
 SMHiggsGGHiggsPPDecayer.h\
 SMHiggsWWDecayer.h

ALL_CC_FILES_2 = \
 SMHiggsFermionsDecayer.cc \
 SMHiggsFermionsPOWHEGDecayer.cc \
 SMHiggsGGHiggsPPDecayer.cc \
 SMHiggsWWDecayer.cc 
