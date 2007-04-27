SUBDIRS = Scalar Vector Tensor StandardModel RSModel Susy UED

noinst_LTLIBRARIES = libHwVertex.la
libHwVertex_la_SOURCES = VertexBase.cc VertexBase.h \
	 VertexBase.fh VertexBase.icc
libHwVertex_la_LIBADD = Scalar/libHwSVertex.la Vector/libHwVVertex.la \
	Tensor/libHwTVertex.la
