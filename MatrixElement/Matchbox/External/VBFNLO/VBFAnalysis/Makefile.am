
if WANT_LIBFASTJET

pkglib_LTLIBRARIES = VBFAnalysis.la
VBFAnalysis_la_LDFLAGS = -module -version-info 11:1:0
VBFAnalysis_la_CPPFLAGS = $(AM_CPPFLAGS) $(FASTJETINCLUDE) \
-I$(FASTJETPATH)

VBFAnalysis_la_CXXFLAGS = $(AM_CXXFLAGS)
#VBFAnalysis_la_CXXFLAGS += -Df2cFortran

VBFAnalysis_la_CFLAGS = $(AM_CFLAGS)
#VBFAnalysis_la_CFLAGS += -Df2cFortran

VBFAnalysis_la_LIBADD = $(FASTJETLIBS) 

VBFAnalysis_la_SOURCES = \
arHistogram.fh arHistogram.h arHistogram.cc arHistogram.icc \
BasicHistogramCollection.h BasicHistogramCollection.cc \
Cutter.h Cutter.cc \
NLOAnalysis.h NLOAnalysis.cc

endif