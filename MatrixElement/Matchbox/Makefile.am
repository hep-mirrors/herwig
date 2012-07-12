SUBDIRS = \
Base Utility Phasespace \
Dipoles InsertionOperators Powheg \
Builtin

pkglib_LTLIBRARIES = HwMatchbox.la
HwMatchbox_la_LDFLAGS = -module -version-info 1:0:0
HwMatchbox_la_LIBADD = \
Base/libHwMatchboxBase.la \
Utility/libHwMatchboxUtility.la \
Phasespace/libHwMatchboxPhasespace.la \
Dipoles/libHwMatchboxDipoles.la \
InsertionOperators/libHwMatchboxInsertionOperators.la \
Powheg/libHwMatchboxPowheg.la \
Builtin/Processes/libHwMatchboxBuiltinProcesses.la \
Builtin/Amplitudes/libHwMatchboxBuiltinAmplitudes.la

HwMatchbox_la_SOURCES =  \
MatchboxFactory.h MatchboxFactory.cc
