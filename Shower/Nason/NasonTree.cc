// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NasonTree class.
//

#include "NasonTree.h"

using namespace Herwig;

NasonTree::NasonTree(vector<NasonBranchingPtr> branchings) {
  for(unsigned int ix=0;ix<branchings.size();++ix) {
    cerr << ix << " " << *branchings[ix]->_particle << "\n";
    cerr << branchings[ix]->_particle->colourLine() << " " 
	 << branchings[ix]->_particle->antiColourLine() << "\n";
  }
}
