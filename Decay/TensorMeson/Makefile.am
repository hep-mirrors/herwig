pkglib_LTLIBRARIES =

if WANT_NEWDECAYERS
pkglib_LTLIBRARIES += HwTMDecay.la
endif

HwTMDecay_la_SOURCES = \
TensorMeson2PScalarDecayer.cc TensorMeson2PScalarDecayer.h TensorMeson2PScalarDecayer.icc \
TensorMesonVectorPScalarDecayer.cc TensorMesonVectorPScalarDecayer.h TensorMesonVectorPScalarDecayer.icc \
TensorMesonVectorVectorDecayer.cc TensorMesonVectorVectorDecayer.h TensorMesonVectorVectorDecayer.icc
HwTMDecay_la_LDFLAGS = -module -version-info 2:0:0
