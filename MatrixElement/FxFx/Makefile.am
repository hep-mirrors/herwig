pkglib_LTLIBRARIES = HwFxFx.la HwFxFxHandler.la
HwFxFx_la_SOURCES = \
FxFxFileReader.cc FxFxFileReader.fh FxFxFileReader.h\
FxFxEventHandler.cc FxFxEventHandler.h FxFxEventHandler.fh\
FxFxReader.cc FxFxReader.fh FxFxReader.h\
FxFx.h
HwFxFx_la_CPPFLAGS = $(AM_CPPFLAGS) $(FASTJETINCLUDE) -I$(FASTJETPATH)
HwFxFx_la_LIBADD = $(FASTJETLIBS) 
HwFxFx_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 2:0:0

HwFxFxHandler_la_SOURCES = \
FxFxHandler.cc FxFxHandler.h
HwFxFxHandler_la_CPPFLAGS = $(AM_CPPFLAGS) $(FASTJETINCLUDE) -I$(FASTJETPATH)
HwFxFxHandler_la_LIBADD = $(FASTJETLIBS) 
HwFxFxHandler_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 2:0:0
