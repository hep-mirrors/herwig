#!/usr/bin/perl -w

# ------------------------------------------------------------------------
# 25-Mar-2002. File: dohtml.pl. 
#  Usage:  dohtml.pl  
#  Note:   it must be run in the same directory ( Herwig++/Doc ).
# Perl script that builds the main HTML file, main.html, containing all
# Herwig++ classes references, organized according to the directories
# where they are defined. Such references point to .html files (in the
# same Herwig++/Doc directory) that are created by this script, using 
# the awk script h2html.awk (which has been copied from Pythia7/Doc).  
# ------------------------------------------------------------------------


use strict;

my $i;
my $j;

if ( @ARGV != 0 ){
  die "---> Usage :   dohtml.pl"; 
} 

open( MAINHTML, ">main.html" ) || die "Cannot open : $! ";
print "\t Writing main.html ... \n";

print MAINHTML <<END_OF_FIRSTPART;
<HTML>
<HEAD>
 <TITLE>Herwig++ classes</TITLE>
</HEAD>
<BODY TEXT="#000000" BGCOLOR="#DDDDDD"  LINK="#A01010" VLINK="#808080" ALINK="#FF0000">
 <H1>Herwig++ classes</H1><BR>
Here the list of all classes, organized according to the directory where they are defined:
<BR>
END_OF_FIRSTPART

chdir "../.";
my @listDirs = <*>;
foreach $i (@listDirs) {
    if ( -d $i ){
	my $namePkg = $i;
	# We exclude some uninteresting directories: lib, my, herwig 
	# that start with a lower case character. Also Doc and Templates
	# are excluded because no class type definition can be found.
	if ( $namePkg =~ /^[A-Z]/  &&  
	     $namePkg ne "Doc"  &&  $namePkg ne "CVS" ) {
	    print MAINHTML '  <BR><BR><FONT COLOR="#0000F0"> <B><FONT SIZE=+2>',
	                   $namePkg.'/','</FONT></B> </FONT><BR><BR>',"\n";
	    chdir $namePkg;
	    print "\t \t Entering directory $namePkg ... \n";
	    my @listHeaders = <*.h>;
	    foreach $j (@listHeaders) {
		my $nameHeader = $j;
		if ( $nameHeader =~ /^[A-Z]/ ) {
		    my ($nameClass,$ignore) = split(/\./,$nameHeader);
		    my $nameHtml = $nameClass.".html";
                    # Use here the awk script, h2html.awk, copied from Pythia7/Doc/
		    system("gawk -f ../Doc/h2html.awk $nameHeader > ../Doc/$nameHtml");
		    print "\t \t \t Writing $nameHtml ... \n";
		    print MAINHTML '    <a href="http:',"$nameHtml",'">',
                                   $nameClass,'</a><BR>',"\n";
		}
	    }
	    chdir "../.";
	    print "\t \t Done directory $namePkg ! \n";
	}
    }
}

chdir "Doc";
print MAINHTML <<END_OF_LASTPART;
</BODY>
</HTML>
END_OF_LASTPART

close( MAINHTML );
print "\t Done main.html ! \n";

#--------------------------------------------------------------------------
