// -*- C++ -*-
//
// NLOJetAmplitudeq2g4.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetAmplitudeq2g4 class.
//

#include "NLOJetAmplitudeq2g4.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/NLOJetRandomWrapper.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEgg2qqbgg.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEgq2qggg.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEgqb2qbggg.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEqg2qggg.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEqbg2qbggg.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEqqb2gggg.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEqbq2gggg.h"

using namespace Herwig;

NLOJetAmplitudeq2g4::NLOJetAmplitudeq2g4() 
  : NLOJetAmplitude<0,2,0>(), theAmplitude(0) {}

NLOJetAmplitudeq2g4::~NLOJetAmplitudeq2g4() {
  if ( theAmplitude ) {
    delete theAmplitude;
    theAmplitude = 0;
  }
}


bool NLOJetAmplitudeq2g4::canHandle(const PDVector& pd) const {

  pair<unsigned int, unsigned int> n = countColoured(pd);
  return n.first == 2 && n.second == 4;

}

Ptr<MatchboxMEBase>::ptr NLOJetAmplitudeq2g4::makeME(const vector<PDVector>& pd) const {

  if ( pd.front()[0]->id() == 21 && pd.front()[1]->id() == 21 )
    return new_ptr(NLOJetMEgg2qqbgg());

  if ( pd.front()[0]->id() == 21 && abs(pd.front()[1]->id()) < 7 ) {
    if ( pd.front()[1]->id() > 0 )
      return new_ptr(NLOJetMEgq2qggg());
    if ( pd.front()[1]->id() < 0 )
      return new_ptr(NLOJetMEgqb2qbggg());
  }

  if ( pd.front()[1]->id() == 21 && abs(pd.front()[0]->id()) < 7 ) {
    if ( pd.front()[0]->id() > 0 )
      return new_ptr(NLOJetMEqg2qggg());
    if ( pd.front()[0]->id() < 0 )
      return new_ptr(NLOJetMEqbg2qbggg());
  }

  if ( pd.front()[2]->id() == 21 && pd.front()[3]->id() == 21 ) {
    if ( pd.front()[0]->id() > 0 )
      return new_ptr(NLOJetMEqqb2gggg());
    if ( pd.front()[0]->id() < 0 )
      return new_ptr(NLOJetMEqbq2gggg());
  }

  assert(false);
  return new_ptr(MatchboxMEBase());

}

double NLOJetAmplitudeq2g4::colourOrdered2(const int* c, size_t) const {
  complex<double> mhv[8];
  complex<double> nmhv[7];
  double res=0.0;

  theAmplitude->matrix_tree_MHV(c[0], c[1], c[2], c[3], c[4], c[5], mhv);
  for(int i=0;i< 8 ;i++){
    res+=norm(mhv[i]);
  }
  theAmplitude->matrix_tree_NMHV(c[0], c[1], c[2], c[3], c[4], c[5], nmhv);
  for(int i=1;i< 7 ;i++){
    res+=norm(nmhv[i]);
  }

  return 2.*res;
}

void NLOJetAmplitudeq2g4::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {
  if ( !theAmplitude ) {
    static NLOJetRandomWrapper dummy;
    theAmplitude = new nlo::ampq2g4(lastInvariants(),dummy);
  }
  NLOJetAmplitude<0,2,0>::prepareAmplitudes(me);
}


double NLOJetAmplitudeq2g4::treeLevel2(const vector<int>& c) const {
  // for one crossed fermion the sign is already included in nlojet
  // as opposed to the 4 and 5 parton processes.
  double sign = 1.; 
  if ( crossingSign() == -1. )
    sign = -1.;
  return sign*theAmplitude->su3_tree(c[0], c[1], c[2], c[3], c[4], c[5]);
}

double NLOJetAmplitudeq2g4::treeOneLoop(const vector<int>&) const {
  throw Exception() << name() << " does not provide one-loop matrix elements"
		    << Exception::abortnow;
  return 0.;
}

double NLOJetAmplitudeq2g4::treeLevelCC(pair<int,int>,
					const vector<int>&) const {
  throw Exception() << name() << " does not provide colour-correlated matrix elements"
		    << Exception::abortnow;
  return 0.;
}

pair<double,Complex> NLOJetAmplitudeq2g4::treeLevelSCC(pair<int,int>,
						       const vector<int>&) const {
  throw Exception() << name() << " does not provide colour-correlated matrix elements"
		    << Exception::abortnow;
  return make_pair(0.,0.);
}

IBPtr NLOJetAmplitudeq2g4::clone() const {
  return new_ptr(*this);
}

IBPtr NLOJetAmplitudeq2g4::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void NLOJetAmplitudeq2g4::persistentOutput(PersistentOStream &) const {}

void NLOJetAmplitudeq2g4::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<NLOJetAmplitudeq2g4,Herwig::NLOJetAmplitude<0,2,0> >
  describeHerwigNLOJetAmplitudeq2g4("Herwig::NLOJetAmplitudeq2g4", "HwMatchboxNLOJet.so");

void NLOJetAmplitudeq2g4::Init() {

  static ClassDocumentation<NLOJetAmplitudeq2g4> documentation
    ("NLOJetAmplitudeq2g4");

}

