BUILT_SOURCES  = Perturbative__all.cc
DISTCLEANFILES = Perturbative__all.cc

Perturbative__all.cc : $(ALL_CC_FILES)
	$(top_srcdir)/cat_with_cpplines $^ > $@

EXTRA_DIST = $(ALL_H_FILES) $(ALL_CC_FILES)

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
DISTCLEANFILES += PerturbativeHiggs__all.cc

PerturbativeHiggs__all.cc : $(ALL_CC_FILES_2)
	$(top_srcdir)/cat_with_cpplines $^ > $@

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
