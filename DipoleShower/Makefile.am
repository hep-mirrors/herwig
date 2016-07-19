SUBDIRS = Base Merging Kernels Kinematics Utility AlphaS

pkglib_LTLIBRARIES = HwDipoleShower.la


HwDipoleShower_la_LDFLAGS = $(AM_LDFLAGS) $(FASTJETINCLUDE) $(FASTJETLIBS)  -module -version-info 4:0:0

HwDipoleShower_la_LIBADD = \
	Base/libHwDipoleShowerBase.la \
	Kernels/libHwDipoleShowerKernels.la \
	Kinematics/libHwDipoleShowerKinematics.la \
        Merging/libHwDipoleShowerMerging.la \
	Utility/libHwDipoleShowerUtility.la \
        $(FASTJETLIBS)

HwDipoleShower_la_SOURCES =  \
	DipoleShowerHandler.h DipoleShowerHandler.cc
