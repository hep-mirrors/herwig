SOURCES := $(wildcard *.cc)

INCLUDEDIR = Config

install-local:
	mkdir -p $(INSTALLSHARE)
	sed -e s:@:$(INSTALLSHARE):g Makefile.install > $(INSTALLSHARE)/Makefile
	sed -e s:-I\$$\(SRCDIR\)/..:-I$(INCLUDEINSTALL):g \
	    -e s:\$$\(SRCDIR\)/lib:$(INSTALLLIB):g \
	    -e s:\$$\(SRCDIR\)/Config:$(INSTALLINCLUDE)/Config:g \
	    Makefile.common > $(INSTALLSHARE)/Makefile.common
	mkdir -p $(INSTALLINCLUDE)/Config
	cp g++dep $(INSTALLINCLUDE)/Config/. 

include ../Config/Makefile.common


# DO NOT DELETE THIS LINE -- g++dep uses it.
# DO NOT PUT ANYTHING AFTER THIS LINE, IT WILL GO AWAY.
