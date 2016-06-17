// -*- C++ -*-
//
// DiagramDrawer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

#include "DiagramDrawer.h"

using namespace Herwig;

// merge two blocks by given line
vector<string> merge(tcPDPtr data, int id,
		     vector<string> block1,
		     vector<string> block2) {

  size_t l1 = block1.front().length();
  size_t l2 = block2.front().length();

  if ( l2 > l1 ) {
    unsigned int count = 0;
    for ( vector<string>::iterator b = block1.begin();
	  b != block1.end(); ++b, ++count ) {
      if ( count != block1.size()/2 ) {
	*b = string(l2-l1,' ') + *b;
      } else {
	*b = string(l2-l1,'-') + *b;
      }
    }
  } else if ( l1 > l2 ) {
    unsigned int count = 0;
    for ( vector<string>::iterator b = block2.begin();
	  b != block2.end(); ++b, ++count ) {
      if ( count != block2.size()/2 ) {
	*b = string(l1-l2,' ') + *b;
      } else {
	*b = string(l1-l2,'-') + *b;
      }
    }
  }

  for ( size_t k = 0; k < block1.size()/2; ++k ) {
    block1[k] = string(1,' ') + block1[k];
  }
  for ( size_t k = block1.size()/2; k < block1.size(); ++k ) {
    block1[k] = string(1,'|') + block1[k];
  }
  for ( size_t k = 0; k <= block2.size()/2; ++k ) {
    block2[k] = string(1,'|') + block2[k];
  }
  for ( size_t k = block2.size()/2+1; k < block2.size(); ++k ) {
    block2[k] = string(1,' ') + block2[k];
  }

  ostringstream lines("");
  lines << "--[" << data->PDGName() << "," << id << "]--";
  string line = lines.str();
  string pad(line.length(),' ');
  line += "|" + string(block1.front().length(),' ');

  vector<string> res;

  for ( vector<string>::iterator b = block1.begin();
	b != block1.end(); ++b )
    res.push_back(pad + *b);
  res.push_back(line);
  for ( vector<string>::iterator b = block2.begin();
	b != block2.end(); ++b )
    res.push_back(pad + *b);

  return res;

}


// draw timelike branch
vector<string> drawTimeLike(const Tree2toNDiagram& d,int line) {

  vector<int> children = d.children(line);
  if ( children.empty() ) {
    ostringstream lines("");
    lines << "--[" << d.allPartons()[line]->PDGName() << "," << line << "]--("
	  << d.externalId(line) << ")";
    return vector<string>(1,lines.str());
  }

  vector<string> left = drawTimeLike(d,children[0]);
  vector<string> right = drawTimeLike(d,children[1]);

  return merge(d.allPartons()[line],line,left,right);

}

// get all timelike blocks
vector<vector<string> > timeBlocks(const Tree2toNDiagram& d) {

  vector<vector<string> > res;
  vector<int> children;
  //children.push_back(0);
  //children.push_back(0);

  do {
    children = d.children(children[0]);
    res.push_back(drawTimeLike(d,children[1]));
  } while ( children[0] != d.nSpace() - 1 );

  size_t maxLength = 0;
  for ( vector<vector<string> >::const_iterator b = res.begin();
	b != res.end(); ++b )
    maxLength = max(maxLength,b->back().length());

  for ( vector<vector<string> >::iterator b = res.begin();
	b != res.end(); ++b ) {
    size_t bLength = b->back().length();
    if ( bLength < maxLength ) {
      unsigned int count = 0;
      for ( vector<string>::iterator c = b->begin();
	    c != b->end(); ++c, ++count ) {
	if ( count != b->size()/2 ) {
	  *c = string(maxLength-bLength,' ') + *c;
	} else {
	  *c = string(maxLength-bLength,'-') + *c;
	}
      }
    }
  }

  for ( vector<vector<string> >::iterator b = res.begin();
	b != res.end(); ++b ) {
    for ( vector<string>::iterator c = b->begin();
	  c != b->end(); ++c ) {
      *c = "  |" + *c;
    }
  }

  return res;

}

void DiagramDrawer::drawDiag(ostream& os,const Tree2toNDiagram& d) {

  os << d.partons()[0]->PDGName() << " "
     << d.partons()[1]->PDGName() << " -> ";
  for ( cPDVector::const_iterator p = d.partons().begin()+2;
	p != d.partons().end(); ++p )
    os << (**p).PDGName() << " ";
  os << "\n\n";

  vector<vector<string> > blocks = timeBlocks(d);

  os << " (0)\n";
  vector<int> children;
  children.push_back(0);
  //children.push_back(0);

  vector<vector<string> >::const_iterator b = blocks.begin();
  do {
    os << "  |\n"
       << "[" << d.allPartons()[children[0]]->PDGName() << "," << children[0] << "]\n"
       << "  |\n";
    for ( vector<string>::const_iterator l = b->begin();
	  l != b->end(); ++l )
      os << *l << "\n";
    children = d.children(children[0]);
    ++b;
  } while ( children[0] != d.nSpace() - 1 );

  os << "  |\n"
     << "[" << d.allPartons()[d.nSpace()-1]->PDGName() << "," << (d.nSpace()-1) << "]\n"
     << "  |\n" << " (1)\n\n" << flush;

}
