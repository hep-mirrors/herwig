SUBDIRS = \
Base Utility Phasespace \
Dipoles InsertionOperators Matching \
Builtin External Tests

if WANT_DIPOLE
pkglib_LTLIBRARIES = HwMatchbox.la
endif

HwMatchbox_la_LDFLAGS = -module -version-info 3:0:0
HwMatchbox_la_LIBADD = \
Base/libHwMatchboxBase.la \
Utility/libHwMatchboxUtility.la \
Phasespace/libHwMatchboxPhasespace.la \
Dipoles/libHwMatchboxDipoles.la \
InsertionOperators/libHwMatchboxInsertionOperators.la \
Matching/libHwMatchboxMatching.la \
Builtin/Amplitudes/libHwMatchboxBuiltinAmplitudes.la \
External/BLHAGeneric/libHwMatchboxExternalBLHAGeneric.la \
Tests/libHwMatchboxTests.la

HwMatchbox_la_SOURCES =  \
MatchboxFactory.h MatchboxFactory.fh MatchboxFactory.cc
