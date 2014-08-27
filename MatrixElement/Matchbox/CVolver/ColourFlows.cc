// -*- C++ -*-

//
// ColourFlowBasis.cpp is part of CVolver, (C) 2013 Simon Plätzer -- simon.plaetzer@desy.de
// CVolver is licenced under version 2 of the GPL, see COPYING for details.
//

#include "ColourFlows.h"

using namespace CVolver;
using namespace std;

set<ColourFlow> ColourFlow::allFlows(const size_t& n) {
  set<ColourFlow> res;
  vector<size_t> perm = identicalPermutation(n);
  do {
    res.insert(ColourFlow(perm));
  } while ( next_permutation(perm.begin(),perm.end()) );
  return res;
}

size_t ColourFlow::scalarProduct(const ColourFlow& other) const {
  assert(other.permutation().size() == thePermutation.size());
  map<size_t,size_t> initialProduct;
  for ( size_t k = 0; k < thePermutation.size(); ++k ) {
    initialProduct[k] = other.permutation()[thePermutation[k]];
  }
  size_t res = 0;
  while ( !initialProduct.empty() ) {
    if ( initialProduct.begin()->first == 
	 initialProduct.begin()->second ) {
      ++res;
      initialProduct.erase(initialProduct.begin());
      continue;
    }
    map<size_t,size_t>::iterator next = initialProduct.find(initialProduct.begin()->second);
    assert(next != initialProduct.end());
    initialProduct.begin()->second = next->second;
    initialProduct.erase(next);
  }
  return res;
}

pair<size_t,size_t> ColourFlow::getTranspositionOf(const ColourFlow& other) const {
  assert(other.permutation().size() == thePermutation.size());
  pair<size_t,size_t> differ(0,0);
  size_t diffCount = 0;
  for ( size_t k = 0; k < thePermutation.size(); ++k ) {
    if ( thePermutation[k] != other.permutation()[k] ) {
      if ( ++diffCount > 2 )
	return make_pair(0,0);
      if ( diffCount == 1 )
	differ.first = k;
      if ( diffCount == 2 )
	differ.second = k;
    }
  }
  if ( differ.first > differ.second )
    std::swap(differ.first,differ.second);
  return differ;
}

bool ColourFlow::isNonZero(const vector<size_t>& colours,
			   const vector<size_t>& antiColours) const {
  assert(colours.size() == antiColours.size() &&
	 colours.size() == thePermutation.size());
  for ( size_t k = 0; k < thePermutation.size(); ++k )
    if ( colours[k] != antiColours[antiColour(k)] )
      return false;
  return true;
}
