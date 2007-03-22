pkglib_LTLIBRARIES =
EXTRA_LTLIBRARIES =

if WANT_LOOPTOOLS
pkglib_LTLIBRARIES += libHwLooptools.la
EXTRA_LTLIBRARIES += libHwLooptoolsCheck.la
endif

theFFSOURCES = ff/ffinit.f ff/ffbndc.f ff/ffxli2.f ff/ffcli2.f \
  ff/fftran.f ff/ffabcd.f ff/ff2dl2.f \
  ff/ffcrr.f ff/ffcxr.f ff/ffcxs3.f ff/ffcxs4.f \
  ff/ffxxyz.f ff/ffcxyz.f ff/ffdxc0.f ff/ffdcc0.f ff/ffdcxs.f \
  ff/ffdel2.f ff/ffdel3.f ff/ffdel4.f \
  ff/ffcel2.f ff/ffcel3.f ff/ffcel4.f \
  ff/ffxa0.f ff/ffxb0.f ff/ffxb1.f ff/ffxb2p.f ff/ffxdb0.f ff/ffxdb1.f \
  ff/ffxc0.f ff/ffxc0i.f ff/ffxc0p0.f ff/ffxc0p.f \
  ff/ffxd0.f ff/ffxd0h.f ff/ffxd0i.f ff/ffxd0m0.f ff/ffxd0p.f ff/ffxdbd.f \
  ff/ffca0.f ff/ffcb0.f ff/ffcb1.f  ff/ffcb2p.f ff/ffcdb0.f \
  ff/ffcc0.f ff/ffcc0p.f \
  ff/ffcd0.f ff/ffcdbd.f ff/fferr.h ff/ff.h ff/ffwarn.h

## ff/ffcb2.f ## breaks with undefined symbols

theLTCOMMON = lt/C0i.F lt/Cget.F lt/D0i.F lt/Dget.F \
  lt/CA0.F lt/CB0.F lt/CB1.F lt/CB00.F lt/CB11.F \
  lt/CDB0.F lt/CDB00.F lt/CDB1.F lt/CDB11.F \
  lt/CC0.F lt/CD0.F lt/CC0i.F lt/CCget.F lt/CD0i.F lt/CDget.F \
  lt/ini.F lt/cache.c lt/lt.h

theLTSOURCES = lt/A0.F lt/B0.F lt/B1.F lt/B00.F lt/B11.F \
  lt/DB0.F lt/DB00.F lt/DB1.F lt/DB11.F \
  lt/C0.F lt/D0.F $(theLTCOMMON)

theLTCHECKSOURCES = ltcheck/auxB.F ltcheck/auxCD.F ltcheck/check.F \
  ltcheck/A0check.F ltcheck/B0check.F ltcheck/B1check.F \
  ltcheck/B00check.F ltcheck/B11check.F ltcheck/DB0check.F \
  ltcheck/DB1check.F ltcheck/DB00check.F ltcheck/DB11check.F \
  ltcheck/C0check.F ltcheck/D0check.F $(LTCOMMON)

## installed already from "include/Makefile.am"
noinst_HEADERS = include/clooptools.h

## no optimization allowed, see LoopTools homepage
## this must be set as FFLAGS to avoid user switching on opt accidentally
## ignore the automake warning!
FFLAGS=-g -O0 

## -DQEXT=DBLE should really be checked in configure script
AM_FFLAGS = -I$(srcdir)/ff -I$(srcdir)/lt -DQEXT=DBLE

libHwLooptools_la_SOURCES = $(theFFSOURCES) $(theLTSOURCES)

libHwLooptools_la_LDFLAGS = -version-info 2:0:0
if NEED_APPLE_FIXES
libHwLooptools_la_LDFLAGS += -Wl,-single_module
endif

libHwLooptoolsCheck_la_SOURCES = $(theFFSOURCES) $(theLTCHECKSOURCES)

EXTRA_PROGRAMS = tester
tester_SOURCES = tester.cc
tester_LDADD = libHwLooptools.la 
