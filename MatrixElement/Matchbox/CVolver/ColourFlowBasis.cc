// -*- C++ -*-
//
// ColourFlowBasis.cc is a part of CVolver
// Copyright (C) 2013-2014 Simon Platzer
//
// CVolver is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ColourFlowBasis class.
//

#include "ColourFlowBasis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace CVolver;

ColourFlowBasis::ColourFlowBasis() {}

ColourFlowBasis::~ColourFlowBasis() {}

IBPtr ColourFlowBasis::clone() const {
  return new_ptr(*this);
}

IBPtr ColourFlowBasis::fullclone() const {
  return new_ptr(*this);
}

void ColourFlowBasis::clear() {
  ColourBasis::clear();
}

map<size_t,vector<vector<size_t> > > 
ColourFlowBasis::basisList(const vector<PDT::Colour>& basisId) const {

  assert(theCrossings.find(basisId) != theCrossings.end());
  const ColourFlowCrossing& crossing = theCrossings.find(basisId)->second;

  size_t n = crossing.nFlows();
  assert(theFlows.find(n) != theFlows.end());

  map<size_t,vector<vector<size_t> > > res;

  const vector<ColourFlow>& flows = theFlows.find(n)->second;

  for ( size_t n = 0; n < flows.size(); ++n ) {
    vector<size_t> fundamental;
    vector<size_t> antiFundamental;
    for ( size_t k = 0; k < flows[n].nLegs(); ++k ) {
      fundamental.push_back(crossing.colourLeg(k));
      size_t ac = flows[n].antiColour(k);
      antiFundamental.push_back(crossing.antiColourLeg(ac));
    }
    res[n].push_back(fundamental);
    res[n].push_back(antiFundamental);
  }

  return res;

}

size_t ColourFlowBasis::prepareBasis(const vector<PDT::Colour>& sub) {

  useMe();

  map<vector<PDT::Colour>,ColourFlowCrossing>::const_iterator cross =
    theCrossings.find(sub);

  if ( cross != theCrossings.end() )
    return cross->second.nFlows();

  ColourFlowCrossing newCrossing(sub,false);
  theCrossings[sub] = newCrossing;

  if ( theFlows.find(newCrossing.nFlows()) != theFlows.end() )
    return newCrossing.nFlows();

  set<ColourFlow> flows = ColourFlow::allFlows(newCrossing.nFlows());

  copy(flows.begin(),flows.end(),
       back_inserter(theFlows[newCrossing.nFlows()]));

  return newCrossing.nFlows();

}

void ColourFlowBasis::readBasisDetails(const vector<PDT::Colour>& sub) {
  prepareBasis(sub);
}

double ColourFlowBasis::scalarProduct(size_t i, size_t j,
				      const vector<PDT::Colour>& abBasis) const {

  if ( largeN() && i != j )
    return 0.;

  assert(theCrossings.find(abBasis) != theCrossings.end());

  size_t n = theCrossings.find(abBasis)->second.nFlows();

  assert(theFlows.find(n) != theFlows.end());

  const ColourFlow& iflow = theFlows.find(n)->second[i];
  const ColourFlow& jflow = theFlows.find(n)->second[j];

  return pow(3.,(double)(iflow.scalarProduct(jflow)));

}

double ColourFlowBasis::tMatrixElement(size_t m, size_t a, size_t b,
				       const vector<PDT::Colour>& aBasis,
				       const vector<PDT::Colour>& bBasis) const {

  assert(theCrossings.find(bBasis) != theCrossings.end());
  const ColourFlowCrossing& bcrossing = theCrossings.find(bBasis)->second;

  assert(theFlows.find(bcrossing.nFlows()) != theFlows.end());
  const ColourFlow& bflow = theFlows.find(bcrossing.nFlows())->second[b];

  assert(theCrossings.find(aBasis) != theCrossings.end());
  const ColourFlowCrossing& acrossing = theCrossings.find(aBasis)->second;

  assert(theFlows.find(acrossing.nFlows()) != theFlows.end());
  const ColourFlow& aflow = theFlows.find(acrossing.nFlows())->second[a];

  size_t bEmitterLine = 
    bBasis[m] == PDT::Colour3 || bBasis[m] == PDT::Colour8 ? 
    bcrossing.colourLine(m) : bcrossing.antiColourLine(m);
  size_t bSpectatorLine = 
    bBasis[m] == PDT::Colour3 || bBasis[m] == PDT::Colour8 ? 
    bflow.antiColour(bEmitterLine) : bflow.colour(bEmitterLine);

  size_t aEmitterLine = 
    bBasis[m] == PDT::Colour3 || bBasis[m] == PDT::Colour8 ? 
    acrossing.colourLine(m) : acrossing.antiColourLine(m);
  size_t aSpectatorLine = 
    bBasis[m] == PDT::Colour3 || bBasis[m] == PDT::Colour8 ? 
    aflow.antiColour(aEmitterLine) : aflow.colour(aEmitterLine);

  assert(aEmitterLine == bEmitterLine);

  if ( bBasis[m] == PDT::Colour3 ) {

    if ( bSpectatorLine != aSpectatorLine ) {
      return 1./2.;
    } else {
      return -1./2./3.;
    }

  }

  if ( bBasis[m] == PDT::Colour3bar ) {

    if ( bSpectatorLine != aSpectatorLine ) {
      return -1./2.;
    } else {
      return 1./2./3.;
    }

  }

  if ( bBasis[m] == PDT::Colour8 ) {

    if ( bSpectatorLine != aSpectatorLine ) {
      return 1.;
    } else {
      return -1.;
    }

  }

  return 0.0;

}

bool ColourFlowBasis::colourConnected(const cPDVector& sub,
				      const vector<PDT::Colour>& basisId,
				      const pair<int,bool>& first,
				      const pair<int,bool>& second, 
				      size_t tensor) const {
 
  assert(theCrossings.find(basisId) != theCrossings.end());
  const ColourFlowCrossing& crossing = theCrossings.find(basisId)->second;

  // translate process to basis ids
  map<cPDVector,map<size_t,size_t> >::const_iterator trans
    = indexMap().find(sub);
  assert(trans != indexMap().end());

  size_t idColoured = first.second ? second.first : first.first;
  idColoured = trans->second.find(idColoured)->second;
  size_t idAntiColoured = first.second ? first.first : second.first;
  idAntiColoured = trans->second.find(idAntiColoured)->second;

  size_t colourLine = crossing.colourLine(idColoured);
  size_t antiColourLine = crossing.antiColourLine(idAntiColoured);

  size_t n = crossing.nFlows();
  assert(theFlows.find(n) != theFlows.end());

  const ColourFlow& iflow = theFlows.find(n)->second[tensor];

  return antiColourLine == iflow.antiColour(colourLine);

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void ColourFlowBasis::persistentOutput(PersistentOStream &) const {}

void ColourFlowBasis::persistentInput(PersistentIStream & , int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<ColourFlowBasis,Herwig::ColourBasis>
  describeColourFlowBasis("CVolver::ColourFlowBasis", 
			  "HwCVolver.so");

void ColourFlowBasis::Init() {

  static ClassDocumentation<ColourFlowBasis> documentation
    ("ColourFlowBasis implements the colour flow basis.",
     "The colour algebra has been performed using CVolver \\cite{Platzer:2013fha}",
     "%\\cite{Platzer:2013fha}\n"
     "\\bibitem{Platzer:2013fha}\n"
     "S.~Platzer,\n"
     "``Summming Large-N Towers in Colour Flow Evolution,''\n"
     "arXiv:1312.2448 [hep-ph].\n"
     "%%CITATION = ARXIV:1312.2448;%%");

}

