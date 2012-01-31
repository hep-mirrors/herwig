SUBDIRS = \
Base Utility Phasespace \
Dipoles InsertionOperators Powheg \
Builtin

pkglib_LTLIBRARIES = HwMatchbox.la
HwMatchbox_la_LDFLAGS = -module -version-info 11:1:0
HwMatchbox_la_LIBADD = \
Base/libHwMatchboxBase.la \
Utility/libHwMatchboxUtility.la \
Phasespace/libHwMatchboxPhasespace.la \
Dipoles/libHwMatchboxDipoles.la \
InsertionOperators/libHwMatchboxInsertionOperators.la \
Powheg/libHwMatchboxPowheg.la \
Builtin/Processes/libHwMatchboxBuiltinProcesses.la

HwMatchbox_la_SOURCES =  \
MatchboxFactory.h MatchboxFactory.cc
