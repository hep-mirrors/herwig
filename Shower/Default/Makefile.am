SUBDIRS = MECorrections .
## libHwDefaultShower.la will go into HwShower.so
noinst_LTLIBRARIES = libHwDefaultShower.la
libHwDefaultShower_la_LIBADD = MECorrections/libHwDefaultMECorrections.la
libHwDefaultShower_la_SOURCES = \
QTildeSudakov.cc QTildeSudakov.h\
QTildeModel.cc QTildeModel.h\
Decay_QtildaShowerKinematics1to2.cc Decay_QtildaShowerKinematics1to2.h  \
IS_QtildaShowerKinematics1to2.cc    IS_QtildaShowerKinematics1to2.h  \
FS_QtildaShowerKinematics1to2.cc    FS_QtildaShowerKinematics1to2.h  \
QtildaShowerKinematics1to2.cc QtildaShowerKinematics1to2.fh\
QtildaShowerKinematics1to2.h  \
QTildeFinder.cc QTildeFinder.h\
QTildeReconstructor.cc QTildeReconstructor.h

libHwDefaultShower_la_CPPFLAGS=$(AM_CPPFLAGS) $(GSLINCLUDE)
