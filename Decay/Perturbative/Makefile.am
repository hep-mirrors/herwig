pkglib_LTLIBRARIES = HwPerturbativeDecay.la

if WANT_NEWDECAYERS
pkglib_LTLIBRARIES += HwPerturbativeHiggsDecay.la
endif

HwPerturbativeDecay_la_SOURCES = \
 SMWZDecayer.cc  SMWZDecayer.fh  SMWZDecayer.h  SMWZDecayer.icc\
 SMTopDecayer.cc SMTopDecayer.fh SMTopDecayer.h SMTopDecayer.icc
HwPerturbativeDecay_la_LDFLAGS = -module

HwPerturbativeHiggsDecay_la_SOURCES = \
 SMHiggsFermionsDecayer.cc  SMHiggsFermionsDecayer.fh  \
 SMHiggsFermionsDecayer.h  SMHiggsFermionsDecayer.icc  \
 SMHiggsGGHiggsPPDecayer.cc SMHiggsGGHiggsPPDecayer.fh \
 SMHiggsGGHiggsPPDecayer.h SMHiggsGGHiggsPPDecayer.icc
HwPerturbativeHiggsDecay_la_LDFLAGS = -module

#SMHiggsWWDecayer.cc  SMHiggsWWDecayer.fh  \
# SMHiggsWWDecayer.h  SMHiggsWWDecayer.icc  
