# h2html.awk
#
# This awk script tries to produce a web page from a C++ header file.
# The idea is to make the documentation more readable, without maintaining
# a separate file.  Typical usage:
#
#	gawk -f h2html.awk list.h > list.html 
#
# Although the result may not be exactly as intended, nothing in the file
# (except #if, #ifdef, #ifndef, #endif) is supposed to be thrown away.
# Major problem: comments associated with stuff that is saved (typedef, #define,
# etc.), is not saved but appears with some unrelated code.
# Improvemnts are gratefully accepted.
#
# Author: Leif Lonnblad. Modified version of classdoc by Dag Michael Bruck,
# Department of Automatic Control, Lund Institute
# of Technology, Box 118, S-221 00 Lund, SWEDEN.  E-mail: dag@control.lth.se
#
FNR == 1	{
                  COL1 = "<FONT COLOR=\"#0000F0\">";
		  COL2 = "<FONT COLOR=\"#000060\">";
		  COL3 = "<FONT COLOR=\"#000000\">";
		  COL4 = "<FONT COLOR=\"#008000\">";
		  COL5 = "<FONT COLOR=\"#A00000\">";
		  COL6 = "<FONT COLOR=\"#F00000\">";
		  COL7 = "<FONT COLOR=\"#505000\">";
		  COLOFF = "</FONT>"
                  print "<HTML>"
		  print "<TITLE>", FILENAME, "</TITLE>"
                  print "<BODY TEXT=\"#000000\" BGCOLOR=\"#DDDDDD\" ",
		    "LINK=\"#A01010\" VLINK=\"#808080\" ALINK=\"#FF0000\">"
		  print COL1, "<B><FONT SIZE=+2>File: </FONT></B>",
		    COLOFF, "<TT>", FILENAME,"</TT><UL>\n";
#		  section("Description", "DESCRIPTION:");
		  Section = "None";
		  SubSect = "None";
		  ContinueStat = 0;
		  SeenCode = 0;
		}
#
# Skip typical emacs mode specifiers
#
/\-\*\- C\+\+ \-\*\-/{next;}

#
# Convert < and > to &lt; and &gt; except in comments.
#
            {
              if ( ! (($1 ~ /^\/\// ) || ( $1 ~ /\#/ ) ) ) {
                gsub(">", "\\&gt;");
                gsub("<", "\\&lt;");
              }
            }

#
# Stop class ducumentation
#
/CLASSDOC[ \t]+OFF/	{
		  while ($0 !~ /CLASSDOC[ \t]+ON/ &&  getline );
		  getline;
		}
#
# Other classdoc directives
#
/CLASSDOC[ \t]+SECTION/	{
  print COL1,  "<P></UL>\n<H1>";
  for ( i=4;i<=NF;i++) print $i;
  print "</H1>\n<UL>\n", COLOFF;
  next;
}
/CLASSDOC[ \t]+SUBSECTION/	{
  print COL2,  "<P></UL>\n<H3>";
  for ( i=4;i<=NF;i++) print $i;
  print "</H3>\n<UL>\n", COLOFF;
  next;
}


#
# #define macro, placed at end of manual page
#
$1 == "#define" {
		  $0 = dropfirst($0);
		  match($0, /[A-Za-z0-9_]+(\([^\)]*\))?/);
		  save("define", substr($0, RSTART, RLENGTH));
		  eatcontinuation();
		  next;
		}
#
# #include files, placed at end
#
$1 == "#include" {
		  save("include", substr($2, 2, length($2)-2));
		  next;
		}
#
# Drop some cpp directives
#
$1 == "}"	{ next; }
$1 == "#if"	{ next; }
$1 == "#ifdef"	{ next; }
$1 == "#undef"	{ next; }
$1 == "#ifndef"	{ next; }
$1 == "#else"	{ next; }
$1 == "#elif"	{ next; }
$1 == "#endif"	{ next; }
$1 == "#error"	{ next; }
#
# Extrnal declarations, placed at end
#
$1 == "extern" {
		  save("extern", dropfirst($0));
		  next;
		}
#
# Typedef declaration, placed at end
#
$1 == "typedef" && Section != "Class" {
		  save("typedef", dropfirst($0));
		  next;
		}
#
# Class declarations, place at end
#
$1 ~ /^class$|^struct$/ && $0 ~ /;$/ {
		  save("classdecl", $1 " " substr($2, 1, length($2)-1));
		  next;
		}
#
# C++ comment (the /* ... */ variant is not recognized)
#
$1 ~ /\/\// && $2 ~ /#include/	{
  next;
}

$1 ~ /^\/\//	{ 
		  if (InComment == 0) print "<BR>&nbsp;";
		  if (NF == 1)
		    print "<P>";
		  else {
		    gsub("<!id>", "<TT><B>");
		    gsub("<!!id>", "</B></TT>");
		    xxx=gensub("<!class>([a-zA-z][a-zA-Z0-9_]*)<!!class>", "<a href=\"http:\\1.html\">\\1</a>", "g");
		    match(xxx, /\/\/[ ]*/);
		    print substr(xxx, RSTART+RLENGTH);
		  }
		  InComment = 1;
		  next;
		}
#
# Class definition
#
$1 ~ /^class$/	&& Section != "Class" {
		  section("Class",
			  COL4 "<P>CLASS <TT><B>" COLOFF COL6 $2 "</FONT></B></TT>");
		  if (NF > 4) {
		    subsect("Base", "Base class:");
		    for (i=4; i < NF; i++)
		      printcode("", $i, "");
		    print "<BR>";
		  };
		  next;
		}

#
# Struct definition - almost a class
#
$1 ~ /^struct$/	&& Section != "Class" {
		  section("Class",
			  COL4 "<P>STRUCT <TT><B>" COLOFF COL6 $2 "</FONT></B></TT>");
		  subsect("Public", "Public members:");
		  next;
		}

#
# Save template definition
#
$1 ~ /^template$/ && Section != "Class" {
      print "</UL>\n";
      printcode("", $0, "");
      print "\n<UL>\n";
      next;
}

# $1 ~ /^enum$/ {
#       print "\n<!code>", $0, "<BR>\n";
#       next;
# }

$1 ~ /^namespace$/ {
  print "</UL>";
  if ( !SeenCode ) print "<HR>";
  print "\n<H1>", COL4, "namespace ", COLOFF, COL6, $2, COLOFF, "</H1>\n<UL>\n";
  SeenCode = 1;
  next;
}

#
# End of class definition
#
$1 == "};"	{
                  print "</UL><HR><UL>";
		  Section = "None";
		  SubSect = "None";
		  InComment = 0;
		  next;
		}
#
# End of namespace definition
#
$1 == "}"	{
                  print "</UL><HR><UL>";
		  Section = "None";
		  SubSect = "None";
		  InComment = 0;
		  next;
		}
#
# Friend declaration, only in classes
#
$1 == "friend"	{
		  subsect("Friend", "Friends");
		  printcode("", dropfirst($0), "");
		  next;
		}
#
# Public, protected and private parts of a class
#
/public:/	{
		  subsect("Public", "Public members:");
		  next;
		}
/protected:/	{
		  subsect("Protected", "Protected members:");
		  next;
		}
/private:/	{
		  subsect("Private", "Private members:");
		  next;
		}
#
# Everything else inside a class (not comments)
#
Section == "Class" || $1 ~ /template/	{ 
		  if (InComment == 1)
		    print "<BR>";
		  if (ContinueStat == 1) {
		    if (NF > 0) {
		      printcode("<BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;",
				substr($0, index($0, $1)), "");
		    }
		  } else {
		    if (NF > 0) printcode("<BR>",substr($0, index($0, $1)), "");
		  }
		  if ( NF > 0 && !(match($0, /[\;\}][ \t]*$/))) {
		    ContinueStat = 1;
		  }
		  else if (ContinueStat == 1) {
		    ContinueStat = 0;
		  }
		  InComment = 0;
		  next;
		}
#
# Blank lines
#
NF == 0		{
		  print "<P>"
		  if (Section == "Code")
		    print "<P>";
		  next;
		}
#
# Everything else (outside classes)
#
		{
		  if (Section != "Code")
		    section("Code", "");
		  if (InComment == 1) {
		    InComment = 0;
		  };
		  if (ContinueStat == 1) {
		    if (NF > 0) printcode("<BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;",
				substr($0, index($0, $1)), "");
		  } else {
		    if (NF > 0) printcode("<BR>",
				substr($0, index($0, $1)), "");
		  }
		  if ( NF > 0 && !(match($0, /[\;\}][ \t]*$/))) {
		    ContinueStat = 1;
		  }
		  else if (ContinueStat == 1) {
		    ContinueStat = 0;
		  }
#		  print "<TT><B>", $0, "</TT></B>";
		  next;
		}
#
# Post-processing: list class declarations, external declarations,
# macro definitions and includec files.
#
END		{
		  print "</UL><HR><UL>";
		  if (issaved("classdecl")) {
		    print COL1, "</UL><H3>CLASS DECLARATIONS</H3><UL>",
		      COLOFF;
		    printsavedcode("classdecl");
		  }
		  if (issaved("extern")) {
		    print COL1, "</UL><H3>EXTERNAL DECLARATIONS</H3><UL>",
		      COLOFF;
		    printsavedcode("extern");
		  }
		  if (issaved("typedef")) {
		    print COL1, "</UL><H3>TYPE DEFINITIONS</H3><UL>",
		      COLOFF;
		    printsavedcode("typedef");
		  }
		  if (issaved("define")) {
		    print COL1, "</UL><H3>DEFINED MACROS</H3><UL>",
		      COLOFF;
		    printsaved("define");
		  }
		  if (issaved("include")) {
		    print COL1, "</UL><H3>INCLUDED FILES</H3><UL>",
		      COLOFF;
		    printsaved("include");
		  }
		  print "</UL>";
                  print "</BODY>"
                  print "</HTML>"
		}
#
# dropfirst(str) - returns string without its first field
#
function dropfirst(str)
{
  if (match(str, /^[ \t]*[^ \t]*[ \t]+/))
    return substr(str, RLENGTH+1);
  else
    return str;
}
#
# save, ... - functions for saving text in an area, printing it, etc.
#
function save(area, str)
{
  TextCount[area]++;
  TextArea[area, TextCount[area]] = str;
}

function issaved(area)
{
  return TextCount[area] > 0;
}

function printsaved(area,		i)
{
  for (i = 1; i <= TextCount[area]; i++)
    print "<TT>", TextArea[area, i], "</TT><BR>";
}
function printsavedcode(area,		i)
{
  for (i = 1; i <= TextCount[area]; i++)
    printcode("<BR>", TextArea[area, i], "");
}
#
# section(sect, heading) - change section, no subsection
#
function section(sect, heading)
{
  if (sect != Section) {
    Section = sect;
    print COL1, "</UL>\n<H1>", heading, "</H1>\n<UL>\n", COLOFF;
    InComment = 0;
  }
  SubSect = "None";
}
#
# subsect(subs, heading) - change subsection
#
function subsect(subs, heading)
{
  if (subs != SubSect) {
    SubSect = subs;
    print COL2,  "</UL>\n<H3>", heading, "</H3>\n<UL>", COLOFF;
    InComment = 0;
  }
}
#
# eatcontinuation() - eat continuation lines (preceding line ended in \)
#
function eatcontinuation()
{
  while ($0 ~ /\\$/)
    getline;
}

function printcode(bef, code, aft) {
  xxx = gensub("([a-zA-Z_][a-zA-Z0-9_]*)([ \t]*\\()",
	       COL5 "\\1" COLOFF "\\2", "g", code);
  yyy = gensub("\\<(auto|bool|char|class|const|double|enum|extern|float|friend|inline|int|long|register|short|signed|static|struct|namespace|typename|template|typedef|virtual|void|volatile|union|unsigned|operator|case|goto|break|catch|continue|delete|do|else|for|of|new|private|protected|public|return|sizeof|switch|this|throw|try|while)\\>",
	       COL4 "\\1" COLOFF, "g", xxx);
  zzz = gensub("\\<(string|vector|map|deque|set|multiset|multimap|iterator|const_iterator)\\>",
	       COL7 "\\1" COLOFF, "g", yyy);
  print bef, "<TT><B>", zzz, "</B></TT>", aft;
}

