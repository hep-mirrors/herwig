DISTFILES = Makefile refman.conf

include ../Config/Makefile.common

doc: refman

distclean: clean

%: %.conf
	doxygen $<
# DO NOT DELETE THIS LINE -- g++dep uses it.
# DO NOT PUT ANYTHING AFTER THIS LINE, IT WILL GO AWAY.
