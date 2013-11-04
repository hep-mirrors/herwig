pkglib_LTLIBRARIES = HwExsample2.la

if WANT_DIPOLE
pkglib_LTLIBRARIES += HwExsampleMatchbox.la
endif

HwExsample2_la_LDFLAGS = -module -version-info 3:0:0
HwExsample2_la_SOURCES = \
GeneralStatistics.h GeneralStatistics.cc \
BinnedStatistics.h BinnedStatistics.cc \
MultiIterationStatistics.h MultiIterationStatistics.cc \
BinSampler.h BinSampler.cc \
ProjectingSampler.h ProjectingSampler.cc \
GeneralSampler.h GeneralSampler.cc \
ExSampler2.h ExSampler2.cc \
SamplingBias.h SamplingBias.cc \
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

HwExsampleMatchbox_la_LDFLAGS = -module -version-info 3:0:0

HwExsampleMatchbox_la_SOURCES =  \
EnhanceNLOContributions.h EnhanceNLOContributions.cc

