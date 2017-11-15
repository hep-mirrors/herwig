SUBDIRS = Base Kernels Kinematics Utility AlphaS Merging 

pkglib_LTLIBRARIES = HwDipoleShower.la

HwDipoleShower_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 8:2:0

HwDipoleShower_la_LIBADD = \
	Base/libHwDipoleShowerBase.la \
	Kernels/libHwDipoleShowerKernels.la \
	Kinematics/libHwDipoleShowerKinematics.la \
	Utility/libHwDipoleShowerUtility.la \
	Merging/libHwDipoleShowerMerging.la

HwDipoleShower_la_SOURCES =  \
	DipoleShowerHandler.h DipoleShowerHandler.fh DipoleShowerHandler.cc


pkglib_LTLIBRARIES += HwKrknloEventReweight.la
HwKrknloEventReweight_la_SOURCES = \
	KrkNLO/KrknloEventReweight.h KrkNLO/KrknloEventReweight.cc

HwKrknloEventReweight_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 1:0:0
