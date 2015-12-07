CLEANFILES = AllInterfaces.h Herwig-refman.tag HerwigDefaults.rpo tagfileThePEG.tag doxywarnings.log

docdir = $(pkgdatadir)/Doc

dist_doc_DATA = refman.conf refman.h HerwigDefaults.in BSMlibs.in
EXTRA_DATA = tagfileThePEG.tag

dist_noinst_SCRIPTS = fixinterfaces.pl

clean-local: html-clean

html-clean:
	rm -rf refman-html Herwig-refman.tag tagfileThePEG.tag

html-local: refman-html

AllInterfaces.h: HerwigDefaults.in BSMlibs.in
	$(top_builddir)/src/Herwig init -L$(top_builddir)/lib

refman-html: refman.conf refman.h AllInterfaces.h tagfileThePEG.tag
	$(builddir)/fixinterfaces.pl AllInterfaces.h
	doxygen $<

tagfileThePEG.tag:
	@if test -e $(THEPEGPATH)/share/ThePEG/Doc/ThePEG-refman.tag; then $(LN_S) $(THEPEGPATH)/share/ThePEG/Doc/ThePEG-refman.tag tagfileThePEG.tag; fi

install-data-local:
	if test -d refman-html; then $(mkdir_p) $(DESTDIR)$(docdir);\
	find refman-html -follow -name '*.*' \
          -exec $(install_sh_DATA) \{\} $(DESTDIR)$(docdir)/\{\} \;; \
	$(install_sh_DATA) tagfileThePEG.tag $(DESTDIR)$(docdir)/tagfileThePEG.tag; \
	$(install_sh_DATA) Herwig-refman.tag $(DESTDIR)$(docdir)/Herwig-refman.tag;fi

uninstall-local:
	rm -rf $(DESTDIR)$(docdir)/refman-html \
	$(DESTDIR)$(docdir)/tagfileThePEG.tag \
	$(DESTDIR)$(docdir)/Herwig-refman.tag
