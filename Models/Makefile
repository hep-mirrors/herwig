SUBDIRS = RSModel StandardModel

all: lib

check: lib

lib:
	for dir in $(SUBDIRS); do cd $$dir ; $(MAKE) -$(MAKEFLAGS) lib; cd .. ; done

clean:
	for dir in $(SUBDIRS); do cd $$dir ; $(MAKE) -$(MAKEFLAGS) clean; cd .. ; done

depend: 
	for dir in $(SUBDIRS); do cd $$dir ; $(MAKE) -$(MAKEFLAGS) depend ; cd .. ; done
install: lib
	for dir in $(SUBDIRS); do cd $$dir ; $(MAKE) VERSION=$(VERSION) -$(MAKEFLAGS) install ; cd .. ; done
doc:
	for dir in $(SUBDIRS); do cd $$dir ; $(MAKE) -$(MAKEFLAGS) doc ; cd .. ; done
