SUBDIRS = \
Base Utility Phasespace \
Dipoles InsertionOperators Matching \
Cuts Scales Builtin External Tests \
ColorFull

if WANT_DIPOLE
noinst_LTLIBRARIES = libHwMatchbox.la
endif

libHwMatchbox_la_LIBADD = \
Base/libHwMatchboxBase.la \
Utility/libHwMatchboxUtility.la \
Phasespace/libHwMatchboxPhasespace.la \
Dipoles/libHwMatchboxDipoles.la \
InsertionOperators/libHwMatchboxInsertionOperators.la \
Matching/libHwMatchboxMatching.la \
External/BLHAGeneric/libHwMatchboxExternalBLHAGeneric.la \
Tests/libHwMatchboxTests.la

libHwMatchbox_la_SOURCES =  \
MatchboxFactory.h MatchboxFactory.fh MatchboxFactory.cc
