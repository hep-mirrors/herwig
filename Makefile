SUBDIRS = Config Utilities Decay Hadronization Shower \
	  MatrixElement PDF src
#Interfaces \

DISTFILES = README Makefile configure configure.in

VERSION = 0.9

TAG = Herwig++-$(VERSION)

.PHONY: lib all depend clean nodebug setup

all: lib

check: lib
	cd src ;  $(MAKE) -$(MAKEFLAGS) check ; cd ..

lib:
	for dir in $(SUBDIRS) lib; do cd $$dir ; $(MAKE) -$(MAKEFLAGS) lib; cd .. ; done

clean: 
	for dir in $(SUBDIRS) lib src; do cd $$dir ; $(MAKE) -$(MAKEFLAGS) clean ; cd .. ; done

distclean: clean
	cd src ; $(MAKE) -$(MAKEFLAGS) distclean ; cd .. 
	rm -f config.cache config.status config.log Config/Makefile.common Config/config.h

depend: 
	for dir in $(SUBDIRS) src ; do cd $$dir ; $(MAKE) -$(MAKEFLAGS) depend ; cd .. ; done

install: lib
	for dir in $(SUBDIRS) lib src ; do cd $$dir ; $(MAKE) VERSION=$(VERSION) -$(MAKEFLAGS) install ; cd .. ; done

doc:
	for dir in $(SUBDIRS) ; do cd $$dir ; $(MAKE) -$(MAKEFLAGS) doc ; cd .. ; done

Doc/h2html: 
	cd Doc; dohtml.pl; cd ..

configure: configure.in
	autoconf

dist: doc
	rm -rf $(TAG)
	mkdir -p $(TAG)/Herwig
	cp $(DISTFILES) $(TAG)/Herwig
	for dir in $(SUBDIRS) src Doc lib Templates ; do cd $$dir ; $(MAKE) -$(MAKEFLAGS) TAGDIR=$(TAG)/Herwig VERSION=$(VERSION) dist ; cd .. ; done
	tar czf $(TAG).tgz $(TAG)
	rm -rf $(TAG)



