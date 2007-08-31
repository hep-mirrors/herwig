CLEANFILES = AllInterfaces.h Herwig++-refman.tag HerwigDefaults.rpo

docdir = $(pkgdatadir)/Doc

dist_doc_DATA = refman.conf refman.h  tagfileThePEG.tag HerwigDefaults.in

doc_SCRIPTS = fixinterfaces.pl

clean-local: html-clean

html-clean:
	rm -rf refman-html refman-latex Herwig++-refman.tag

html-local: refman-html

pdf-local: refman-html
	if test -d refman-latex; then cd refman-latex; make refman.pdf; fi

AllInterfaces.h: HerwigDefaults.in
	$(top_builddir)/src/Herwig++ init -L$(top_builddir)/lib

refman-html: refman.conf refman.h AllInterfaces.h refmanThePEG.conf
	./fixinterfaces.pl AllInterfaces.h
	doxygen $<

tagfileThePEG.tag:
	@if test -e $(THEPEGPATH)/share/ThePEG/Doc/ThePEG-refman.tag; then $(LN_S) $(THEPEGPATH)/share/ThePEG/Doc/ThePEG-refman.tag tagfileThePEG.tag; fi

#install-data-local:
#	if test -d refman-html; then $(mkdir_p) $(DESTDIR)$(docdir);\
#	find refman-html -follow -name '*.*' \
#          -exec $(install_sh_DATA) \{\} $(DESTDIR)$(docdir)/\{\} \;; \
#	$(install_sh_DATA) ThePEG-refman.tag $(DESTDIR)$(docdir)/ThePEG-refman.tag; fi
#	if test -d refman-latex -a -f refman-latex/refman.pdf; then \
#	$(install_sh_DATA) refman-latex/refman.pdf $(DESTDIR)$(docdir)/refman.pdf; fi

#uninstall-local:
#	rm -rf $(DESTDIR)$(docdir)/refman-html \
#	$(DESTDIR)$(docdir)/ThePEG-refman.tag \
#	$(DESTDIR)$(docdir)/refman.pdf
