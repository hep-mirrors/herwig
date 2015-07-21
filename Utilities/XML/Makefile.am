noinst_LTLIBRARIES = libHwXML.la
libHwXML_la_SOURCES = \
Element.cc \
Element.h \
ElementIO.cc \
ElementIO.h \
ElementIO.tcc
check_PROGRAMS = xml_test
xml_test_SOURCES = tests/xmlTests.cc 
xml_test_LDADD = $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB) $(BOOST_UNIT_TEST_FRAMEWORK_LIB) libHwXML.la
TESTS = xml_test
