// -*- C++ -*-
//
// NLOJetAmplitudeq2g2.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetAmplitudeq2g2 class.
//

#include "NLOJetAmplitudeq2g2.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/NLOJetRandomWrapper.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEgg2qqb.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEgq2qg.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEgqb2qbg.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEqg2qg.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEqbg2qbg.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEqqb2gg.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEqbq2gg.h"


using namespace Herwig;

NLOJetAmplitudeq2g2::NLOJetAmplitudeq2g2() 
  : NLOJetAmplitude<0,2,0>(), theAmplitude(0) {}

NLOJetAmplitudeq2g2::~NLOJetAmplitudeq2g2() {
  if ( theAmplitude ) {
    delete theAmplitude;
    theAmplitude = 0;
  }
}

bool NLOJetAmplitudeq2g2::canHandle(const PDVector& pd) const {

  pair<unsigned int, unsigned int> n = countColoured(pd);
  return n.first == 2 && n.second == 2;

}

Ptr<MatchboxMEBase>::ptr NLOJetAmplitudeq2g2::makeME(const vector<PDVector>& pd) const {

  if ( pd.front()[0]->id() == 21 && pd.front()[1]->id() == 21 )
    return new_ptr(NLOJetMEgg2qqb());

  if ( pd.front()[0]->id() == 21 && abs(pd.front()[1]->id()) < 7 ) {
    if ( pd.front()[1]->id() > 0 )
      return new_ptr(NLOJetMEgq2qg());
    if ( pd.front()[1]->id() < 0 )
      return new_ptr(NLOJetMEgqb2qbg());
  }

  if ( pd.front()[1]->id() == 21 && abs(pd.front()[0]->id()) < 7 ) {
    if ( pd.front()[0]->id() > 0 )
      return new_ptr(NLOJetMEqg2qg());
    if ( pd.front()[0]->id() < 0 )
      return new_ptr(NLOJetMEqbg2qbg());
  }

  if ( pd.front()[2]->id() == 21 && pd.front()[3]->id() == 21 ) {
    if ( pd.front()[0]->id() > 0 )
      return new_ptr(NLOJetMEqqb2gg());
    if ( pd.front()[0]->id() < 0 )
      return new_ptr(NLOJetMEqbq2gg());
  }

  assert(false);
  return new_ptr(MatchboxMEBase());

}

double NLOJetAmplitudeq2g2::colourOrdered2(const int* c, size_t) const {
  complex<double> tmp[2];
  double res=0.0;
  
  theAmplitude->matrix_tree_ppmm(c[0], c[1], c[2], c[3], tmp);
  res+=norm(tmp[0]);
  theAmplitude->matrix_tree_pmpm(c[0], c[1], c[2], c[3], tmp);
  res+=norm(tmp[0]);

  return 2.*res;
}

void NLOJetAmplitudeq2g2::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {
  if ( !theAmplitude ) {
    static NLOJetRandomWrapper dummy;
    theAmplitude = new nlo::ampq2g2(lastInvariants(),dummy);
  }
  NLOJetAmplitude<0,2,0>::prepareAmplitudes(me);
}

double NLOJetAmplitudeq2g2::treeLevel2(const vector<int>& c) const {
  return theAmplitude->su3_tree(c[0], c[1], c[2], c[3]);
}

double NLOJetAmplitudeq2g2::treeOneLoop(const vector<int>& c) const {
  // sign conventions different here, again
  double sign = 1.; 
  if ( crossingSign() == -1. )
    sign = -1.;
  return sign*theAmplitude->su3_1loop(c[0], c[1], c[2], c[3]);
}

double NLOJetAmplitudeq2g2::treeLevelCC(pair<int,int> ij,
					const vector<int>& c) const {
  return theAmplitude->su3_cc(ij.first,ij.second,c[0], c[1], c[2], c[3]);
}

pair<double,Complex> NLOJetAmplitudeq2g2::treeLevelSCC(pair<int,int> ij,
						       const vector<int>& c) const {
  return make_pair(theAmplitude->su3_cc(ij.first,ij.second,c[0], c[1], c[2], c[3]),0.);
}

IBPtr NLOJetAmplitudeq2g2::clone() const {
  return new_ptr(*this);
}

IBPtr NLOJetAmplitudeq2g2::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void NLOJetAmplitudeq2g2::persistentOutput(PersistentOStream &) const {}

void NLOJetAmplitudeq2g2::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<NLOJetAmplitudeq2g2,Herwig::NLOJetAmplitude<0,2,0> >
  describeHerwigNLOJetAmplitudeq2g2("Herwig::NLOJetAmplitudeq2g2", "HwMatchboxNLOJet.so");

void NLOJetAmplitudeq2g2::Init() {

  static ClassDocumentation<NLOJetAmplitudeq2g2> documentation
    ("NLOJetAmplitudeq2g2");

}

