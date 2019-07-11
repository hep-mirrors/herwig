pkglib_LTLIBRARIES = HwFxFx.la HwFxFxHandler.la
HwFxFx_la_SOURCES = \
FxFxFileReader.cc FxFxFileReader.fh FxFxFileReader.h\
FxFxEventHandler.cc FxFxEventHandler.h\
FxFxReader.cc FxFxReader.fh FxFxReader.h\
FxFx.h
HwFxFx_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 1:0:0

HwFxFxHandler_la_SOURCES = \
FxFxHandler.cc FxFxHandler.h
HwFxFxHandler_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 1:0:0

if HAVE_RIVET 
  pkglib_LTLIBRARIES += HwFxFxAnalysis.la
  HwFxFxAnalysis_la_SOURCES = \
  FxFxAnalysis.cc FxFxAnalysis.h
  HwFxFxAnalysis_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 1:0:0
endif
