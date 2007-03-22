pkglib_LTLIBRARIES = HwSMVertex.la
HwSMVertex_la_SOURCES = \
SMFFGVertex.cc SMFFGVertex.h SMFFGVertex.icc \
SMFFHVertex.cc SMFFHVertex.h SMFFHVertex.icc \
SMFFPVertex.cc SMFFPVertex.h SMFFPVertex.icc \
SMFFWVertex.cc SMFFWVertex.h SMFFWVertex.icc \
SMFFZVertex.cc SMFFZVertex.h SMFFZVertex.icc \
SMGGGGVertex.cc SMGGGGVertex.h SMGGGGVertex.icc \
SMGGGVertex.cc SMGGGVertex.h SMGGGVertex.icc \
SMWWHVertex.cc SMWWHVertex.h SMWWHVertex.icc \
SMWWWVertex.cc SMWWWVertex.h SMWWWVertex.icc \
SMWWWWVertex.cc SMWWWWVertex.h SMWWWWVertex.icc 

EXTRA_HwSMVertex_la_SOURCES = \
SMHGGVertex.cc SMHGGVertex.h SMHGGVertex.icc SMHGGVertex.fh \
SMHPPVertex.cc SMHPPVertex.h SMHPPVertex.icc SMHPPVertex.fh

if WANT_LOOPTOOLS
HwSMVertex_la_SOURCES += \
SMHGGVertex.cc SMHGGVertex.h SMHGGVertex.icc SMHGGVertex.fh \
SMHPPVertex.cc SMHPPVertex.h SMHPPVertex.icc SMHPPVertex.fh
endif

HwSMVertex_la_LDFLAGS = -module -version-info 2:0:0
