SUBDIRS = Base Kernels Kinematics Utility AlphaS Merging SpinCorrelations 

pkglib_LTLIBRARIES = HwDipoleShower.la

HwDipoleShower_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 11:0:0

HwDipoleShower_la_LIBADD = \
	Base/libHwDipoleShowerBase.la \
	Kernels/libHwDipoleShowerKernels.la \
	Kinematics/libHwDipoleShowerKinematics.la \
	Utility/libHwDipoleShowerUtility.la \
	Merging/libHwDipoleShowerMerging.la \
	SpinCorrelations/libHwDipoleShowerSpinCorrelations.la 

HwDipoleShower_la_SOURCES =  \
	DipoleShowerHandler.h DipoleShowerHandler.fh DipoleShowerHandler.cc


pkglib_LTLIBRARIES += HwKrknloEventReweight.la
HwKrknloEventReweight_la_SOURCES = \
	KrkNLO/KrknloEventReweight.h KrkNLO/KrknloEventReweight.cc

HwKrknloEventReweight_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 3:0:0
