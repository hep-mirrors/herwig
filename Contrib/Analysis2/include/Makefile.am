BUILT_SOURCES = done-all-links

DIRLINKS = \
$(top_srcdir)/Histogram2 \
$(top_srcdir)/Analysis2

CLEANFILES = done-all-links 

done-all-links: $(DIRLINKS)
	mkdir -p Analysis2
	$(LN_S) -f $(addprefix ../, $(DIRLINKS)) Analysis2
	touch done-all-links

install-data-local:
	find Analysis2 -follow \( -name '*.h' -or -name '*.cc' -or -name '*.tcc' \) \
          -exec $(install_sh_DATA) \{\} $(DESTDIR)$(includedir)/\{\} \;

uninstall-local:
	rm -rf $(DESTDIR)$(includedir)/Analysis2

clean-local:
	rm -rf Analysis2
