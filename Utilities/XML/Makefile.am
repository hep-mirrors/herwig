noinst_LTLIBRARIES = libHwXML.la
libHwXML_la_SOURCES = \
Element.cc \
Element.h \
ElementIO.cc \
ElementIO.h \
ElementIO.tcc
check_PROGRAMS = xml_test
xml_test_SOURCES = tests/xmlTests.cc 
xml_test_LDADD = $(BOOST_UNIT_TEST_FRAMEWORK_LIBS) $(BOOST_FILESYSTEM_LIBS) $(BOOST_SYSTEM_LIBS) libHwXML.la
xml_test_LDFLAGS = $(AM_LDFLAGS) -export-dynamic $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
xml_test_CPPFLAGS = $(AM_CPPFLAGS) $(BOOST_CPPFLAGS) -DHERWIG_PKGDATADIR="\"$(pkgdatadir)\"" -DHERWIG_PKGLIBDIR="\"$(pkglibdir)\"" -DTHEPEG_PKGLIBDIR="\"$(THEPEGLIBPATH)\""
TESTS = xml_test
