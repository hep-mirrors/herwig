SUBDIRS = Base Kernels Kinematics Utility AlphaS

if WANT_DIPOLE
pkglib_LTLIBRARIES = HwDipoleShower.la
endif

HwDipoleShower_la_LDFLAGS = -module -version-info 1:0:0

HwDipoleShower_la_LIBADD = \
	Base/libHwDipoleShowerBase.la \
	Kernels/libHwDipoleShowerKernels.la \
	Kinematics/libHwDipoleShowerKinematics.la \
	Utility/libHwDipoleShowerUtility.la

HwDipoleShower_la_SOURCES =  \
	DipoleShowerHandler.h DipoleShowerHandler.cc

