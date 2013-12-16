// -*- C++ -*-
//
// NLOJetAmplitudeq4g2.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetAmplitudeq4g2 class.
//

#include "NLOJetAmplitudeq4g2.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/NLOJetRandomWrapper.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEqbq2kkbgg.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEqqb2kkbgg.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEkbq2qkbgg.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEkqb2kqbgg.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEkbqb2kbqbgg.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEkq2kqgg.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEgq2kqkbg.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEgqb2kkbqbg.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEqg2kqkbg.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEqbg2kkbqbg.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEgg2kqkbqb.h"

using namespace Herwig;

NLOJetAmplitudeq4g2::NLOJetAmplitudeq4g2() 
  : NLOJetAmplitude<0,2,0>(), theAmplitude(0) {}

NLOJetAmplitudeq4g2::~NLOJetAmplitudeq4g2() {
  if ( theAmplitude ) {
    delete theAmplitude;
    theAmplitude = 0;
  }
}

bool NLOJetAmplitudeq4g2::canHandle(const PDVector& pd) const {

  pair<unsigned int, unsigned int> n = countColoured(pd);
  return n.first == 4 && n.second == 2;

}

Ptr<MatchboxMEBase>::ptr NLOJetAmplitudeq4g2::makeME(const vector<PDVector>& pd) const {

  if ( pd.front()[0]->id() == 21 &&
       pd.front()[1]->id() == 21 )
    return new_ptr(NLOJetMEgg2kqkbqb());

  if ( pd.front()[0]->id() == 21 ) {
    if ( pd.front()[1]->id() > 0 )
      return new_ptr(NLOJetMEgq2kqkbg());
    if ( pd.front()[1]->id() < 0 )
      return new_ptr(NLOJetMEgqb2kkbqbg());
  }

  if ( pd.front()[1]->id() == 21 ) {
    if ( pd.front()[0]->id() > 0 )
      return new_ptr(NLOJetMEqg2kqkbg());
    if ( pd.front()[0]->id() < 0 )
      return new_ptr(NLOJetMEqbg2kkbqbg());
  }

  if ( pd.front()[0]->id() + pd.front()[1]->id() == 0 ) {
    if ( pd.front()[0]->id() < 0 )
      return new_ptr(NLOJetMEqbq2kkbgg());
    if ( pd.front()[0]->id() > 0 )
      return new_ptr(NLOJetMEqqb2kkbgg());
  }

  if ( pd.front()[0]->id() < 0 &&
       pd.front()[1]->id() > 0 ) {
    return new_ptr(NLOJetMEkbq2qkbgg());
  }

  if ( pd.front()[0]->id() > 0 &&
       pd.front()[1]->id() < 0 ) {
    return new_ptr(NLOJetMEkqb2kqbgg());
  }

  if ( pd.front()[0]->id() < 0 &&
       pd.front()[1]->id() < 0 ) {
    return new_ptr(NLOJetMEkbqb2kbqbgg());
  }

  if ( pd.front()[0]->id() > 0 &&
       pd.front()[1]->id() > 0 ) {
    return new_ptr(NLOJetMEkq2kqgg());
  }

  assert(false);
  return new_ptr(MatchboxMEBase());

}

double NLOJetAmplitudeq4g2::colourOrdered2(const int* c, size_t) const {
  static const int LEND= -999;

  std::complex<double> tmp[16];
  double res=0.0;
 
  if(c[2] == LEND){
    theAmplitude->matrix_tree_02(c[0], c[1], c[3], c[6], c[4], c[5], tmp);
  }
  else if(c[3] == LEND){
    theAmplitude->matrix_tree_11(c[0], c[2], c[4], c[6], c[1], c[5], tmp);
  }
  else if(c[4] == LEND){
    theAmplitude->matrix_tree_20(c[0], c[3], c[5], c[6], c[1], c[2], tmp);
  }
  else{
    assert(false && "not proper values of colourOrdered2 argument");
    return 0.;
  }
  for(int i=0;i<16;i++)  
    res +=norm(tmp[i]); 
  
  return 2.*res;
}

void NLOJetAmplitudeq4g2::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {
  if ( !theAmplitude ) {
    static NLOJetRandomWrapper dummy;
    theAmplitude = new nlo::ampq4g2(lastInvariants(),dummy);
  }
  NLOJetAmplitude<0,2,0>::prepareAmplitudes(me);
}


double NLOJetAmplitudeq4g2::treeLevel2(const vector<int>& c) const {
  static double res[2];
  theAmplitude->su3_tree(c[0], c[1], c[2], c[3], c[4], c[5], res);
  // for one crossed fermion the sign is already included in nlojet
  // as opposed to the 4 and 5 parton processes.
  double sign = 1.; 
  if ( crossingSign() == -1. )
    sign = -1.;
  return identicalQuarks() == 2 ? sign*res[1] : sign*res[0];
}

double NLOJetAmplitudeq4g2::treeOneLoop(const vector<int>&) const {
  throw Exception() << name() << " does not provide one-loop matrix elements"
		    << Exception::abortnow;
  return 0.;
}

double NLOJetAmplitudeq4g2::treeLevelCC(pair<int,int>,
					   const vector<int>&) const {
  throw Exception() << name() << " does not provide colour-correlated matrix elements"
		    << Exception::abortnow;
  return 0.;
}

pair<double,Complex> NLOJetAmplitudeq4g2::treeLevelSCC(pair<int,int>,
						       const vector<int>&) const {
  throw Exception() << name() << " does not provide colour-correlated matrix elements"
		    << Exception::abortnow;
  return make_pair(0.,0.);
}

IBPtr NLOJetAmplitudeq4g2::clone() const {
  return new_ptr(*this);
}

IBPtr NLOJetAmplitudeq4g2::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void NLOJetAmplitudeq4g2::persistentOutput(PersistentOStream & ) const {
}

void NLOJetAmplitudeq4g2::persistentInput(PersistentIStream &, int) {
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<NLOJetAmplitudeq4g2,Herwig::NLOJetAmplitude<0,2,0> >
  describeHerwigNLOJetAmplitudeq4g2("Herwig::NLOJetAmplitudeq4g2", "HwMatchboxNLOJet.so");

void NLOJetAmplitudeq4g2::Init() {

  static ClassDocumentation<NLOJetAmplitudeq4g2> documentation
    ("NLOJetAmplitudeq4g2");

}

