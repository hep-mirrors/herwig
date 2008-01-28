noinst_LTLIBRARIES =
EXTRA_LTLIBRARIES =

if WANT_LOOPTOOLS
EXTRA_LTLIBRARIES += libHwLooptoolsCFC.la 
noinst_LTLIBRARIES += libHwLooptoolsXFC.la ## libHwLooptoolsCFC.la 
## EXTRA_LTLIBRARIES += libHwLooptoolsCheck2.la
endif

theXFCSOURCES = A/ffca0.F A/ffxa0.F A/A0.F A/A00.F \
B/ffxdb0.F B/Bcoeffa.F B/ffcb0.F B/Bcoeff.F B/ffxdb1.F \
B/ffdel2.F B/ffcel2.F B/ffcb1.F B/ffcdb0.F \
B/Bcoeffb.F B/ffxb1.F B/Bget.F B/ffxb0.F B/ffcb2p.F B/ffxb2p.F \
C/ffxc0p.F C/ffcc0p.F C/ffdel3.F C/Cget.F C/ffcel3.F C/ffdcc0.F \
C/ffxc0.F C/ffdxc0.F C/ffxc0p0.F C/C0.F C/ffcc0.F C/ffxc0i.F \
D/ffxdbd.F D/ffxd0i.F D/ffcdbd.F D/ffcd0.F D/D0.F D/ffxd0p.F D/ffcel4.F \
D/ffxd0.F D/Dget.F D/ffdel4.F D/ffxd0m0.F D/ffxd0h.F \
E/Ecoeffb.F E/E0.F E/ffxe0.F E/ffdel5.F E/Ecoeffa.F E/Eget.F \
util/ini.F util/ffcxs3.F util/GaussPivot.F util/Dump.F util/ffbndc.F \
util/ffcxr.F util/ffinit.F util/ffabcd.F util/Li2.F util/ff2dl2.F \
util/cache.c util/ffcxyz.F util/ffxxyz.F util/ffxli2.F util/ffdcxs.F \
util/ffcxs4.F util/fftran.F util/CGaussPivot.F util/auxCD.F util/ffcli2.F \
util/ffcrr.F clooptools.cc
##util/qcomplex.f90 

## the following need -DCOMPLEXPARA
theCFCSOURCES = A/A0.F A/A00.F B/Bget.F B/Bcoeffa.F C/Cget.F D/Dget.F \
E/E0.F E/Eget.F E/Ecoeffa.F E/Ecoeffb.F \
util/GaussPivot.F util/Dump.F util/Li2.F

commonCPPFLAGS = $(AM_CPPFLAGS) -I$(srcdir)/include
libHwLooptoolsXFC_la_SOURCES = $(theXFCSOURCES)
libHwLooptoolsXFC_la_CPPFLAGS = $(commonCPPFLAGS) 
libHwLooptoolsXFC_la_CFLAGS = $(AM_CFLAGS) -std=c99
libHwLooptoolsXFC_la_FFLAGS = $(AM_FFLAGS) ## -fdefault-integer-8

libHwLooptoolsCFC_la_SOURCES = $(theCFCSOURCES)
libHwLooptoolsCFC_la_CPPFLAGS = $(commonCPPFLAGS) -DCOMPLEXPARA

## installed already from "include/Makefile.am"
## noinst_HEADERS = include/clooptools.h

##if NEED_APPLE_FIXES
##libHwLooptoolsXFC_la_LDFLAGS = -Wl,-single_module
##libHwLooptoolsCFC_la_LDFLAGS = -Wl,-single_module
##endif

## libHwLooptoolsHelper2_la_SOURCES = clooptools.cc
EXTRA_PROGRAMS = tester
tester_SOURCES = tester.cc
tester_LDADD = libHwLooptoolsXFC.la $(FLIBS)
