#SUBDIRS = Config Shower Hadronization Decay MatrixElement Utilities  
SUBDIRS = Config Shower Hadronization Decay Utilities  

.PHONY: lib all depend clean nodebug setup

all: current

current: lib
	cd src ;  $(MAKE) -$(MAKEFLAGS) current ; cd ..

lib: setup
	for dir in $(SUBDIRS) lib; do cd $$dir ; $(MAKE) -$(MAKEFLAGS) lib CXXDEBUGFLAG=-g CXXOPTFLAG=""; cd .. ; done

nodebug: setup
	for dir in $(SUBDIRS) lib; do cd $$dir ; $(MAKE) -$(MAKEFLAGS) lib ; cd .. ; done

clean:
	for dir in $(SUBDIRS) lib src; do cd $$dir ; $(MAKE) -$(MAKEFLAGS) clean ; cd .. ; done;
	cd Config/; rm Makefile.common configure config.h config.cache config.log config.status; cd ../src/; rm *.rpo *.run *.out *.log *.tex; cd ../Doc/; rm *.html; cd ..;

depend: setup
	for dir in $(SUBDIRS) src; do cd $$dir ; $(MAKE) -$(MAKEFLAGS) depend ; cd .. ; done

setup: Config/Makefile.common Config/config.h

Config/Makefile.common: Config/Makefile.common.in Config/configure
	cd Config/; ./configure; cd ../

Config/config.h: Config/config.h.in Config/configure
	cd Config/; ./configure; cd ../

Config/configure: Config/configure.in
	cd Config/; autoconf







