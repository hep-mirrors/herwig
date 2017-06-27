SUBDIRS = CellGrids

pkglib_LTLIBRARIES = HwSampling.la

HwSampling_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 6:0:0
HwSampling_la_SOURCES = \
GeneralStatistics.h GeneralStatistics.cc \
MultiIterationStatistics.h MultiIterationStatistics.cc \
BinSampler.h BinSampler.cc \
GeneralSampler.h GeneralSampler.cc \
MonacoSampler.h MonacoSampler.cc \
Remapper.h Remapper.cc \
exsample/cell.icc \
exsample/generator.h \
exsample/linear_interpolator.icc \
exsample/config.h \
exsample/selectors.h \
exsample/adaption_info.h \
exsample/selectors.icc \
exsample/adaption_info.icc \
exsample/exponential_generator.h \
exsample/exponential_generator.icc \
exsample/generator.icc \
exsample/statistics.h \
exsample/linear_interpolator.h \
exsample/cell.h \
exsample/utility.h \
exsample/binary_tree.h \
exsample/statistics.icc

HwSampling_la_LIBADD = \
CellGrids/libHwCellGrids.la

EXTRA_DIST = herwig-mergegrids.py

install-exec-local:
	$(install_sh_SCRIPT) $(srcdir)/herwig-mergegrids.py $(DESTDIR)$(bindir)/herwig-mergegrids

uninstall-local:
	rm -f  $(DESTDIR)$(bindir)/herwig-mergegrids
