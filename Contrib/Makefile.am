EXTRA_DIST = \
HiggsPair \
HiggsPairOL \
MultiWeight \
ShowerVeto

dist-hook:
	rm -rf `find $(distdir) -name '.svn' -or -name 'Makefile' -or -name '.hg'`

all:  
	bash make_makefiles.sh
