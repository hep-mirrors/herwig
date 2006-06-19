BUILT_SOURCES = done-all-links
CLEANFILES = done-all-links \
Herwig++/Config \
Herwig++/Decay \
Herwig++/Hadronization \
Herwig++/Helicity \
Herwig++/Interfaces \
Herwig++/MatrixElement \
Herwig++/Models \
Herwig++/PDF \
Herwig++/PDT \
Herwig++/Shower \
Herwig++/UnderlyingEvent \
Herwig++/Utilities \
Herwig++/Analysis \
Herwig++/Looptools/clooptools.h

## replace Herwig++ with $(PACKAGE)
## replace find with explicit list of directories

my_directories=Config Decay Hadronization Helicity Interfaces \
MatrixElement Models PDF PDT Shower UnderlyingEvent Utilities Analysis
## Looptools is dealt with separately

done-all-links:
	mkdir -p Herwig++
	cd Herwig++; \
	for i in $(my_directories); \
	do \
	$(LN_S) -f ../$(top_srcdir)/$$i ; \
	done
	mkdir -p Herwig++/Looptools
	cd Herwig++/Looptools; \
	$(LN_S) -f ../../$(top_srcdir)/Looptools/include/clooptools.h
	echo "stamp" > done-all-links

install-data-local:
	find Herwig++ -follow -name '*.h' \
          -exec $(install_sh_DATA) \{\} $(DESTDIR)$(includedir)/\{\} \;
	find Herwig++ -follow -name '*.icc' \
          -exec $(install_sh_DATA) \{\} $(DESTDIR)$(includedir)/\{\} \;
	find Herwig++ -follow -name '*.fh'  \
          -exec $(install_sh_DATA) \{\} $(DESTDIR)$(includedir)/\{\} \;

uninstall-local:
	rm -rf $(DESTDIR)$(includedir)/Herwig++

