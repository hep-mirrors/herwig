// -*- C++ -*-
//
// NLOJetAmplitudeq4.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetAmplitudeq4 class.
//

#include "NLOJetAmplitudeq4.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/NLOJetRandomWrapper.h"

#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEqbq2kkb.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEqqb2kkb.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEkbq2qkb.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEkqb2kqb.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEkbqb2kbqb.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEkq2kq.h"

using namespace Herwig;

NLOJetAmplitudeq4::NLOJetAmplitudeq4() 
  : NLOJetAmplitude<0,2,0>(), theAmplitude(0) {}

NLOJetAmplitudeq4::~NLOJetAmplitudeq4() {
  if ( theAmplitude ) {
    delete theAmplitude;
    theAmplitude = 0;
  }
}

bool NLOJetAmplitudeq4::canHandle(const PDVector& pd) const {

  pair<unsigned int, unsigned int> n = countColoured(pd);
  return n.first == 4 && n.second == 0;

}

Ptr<MatchboxMEBase>::ptr NLOJetAmplitudeq4::makeME(const vector<PDVector>& pd) const {

  if ( pd.front()[0]->id() + pd.front()[1]->id() == 0 ) {
    if ( pd.front()[0]->id() < 0 )
      return new_ptr(NLOJetMEqbq2kkb());
    if ( pd.front()[0]->id() > 0 )
      return new_ptr(NLOJetMEqqb2kkb());
  }

  if ( pd.front()[0]->id() < 0 &&
       pd.front()[1]->id() > 0 ) {
    return new_ptr(NLOJetMEkbq2qkb());
  }

  if ( pd.front()[0]->id() > 0 &&
       pd.front()[1]->id() < 0 ) {
    return new_ptr(NLOJetMEkqb2kqb());
  }

  if ( pd.front()[0]->id() < 0 &&
       pd.front()[1]->id() < 0 ) {
    return new_ptr(NLOJetMEkbqb2kbqb());
  }

  if ( pd.front()[0]->id() > 0 &&
       pd.front()[1]->id() > 0 ) {
    return new_ptr(NLOJetMEkq2kq());
  }

  assert(false);
  return new_ptr(MatchboxMEBase());

}

double NLOJetAmplitudeq4::colourOrdered2(const int* c, size_t) const {
  complex<double> tmp[1];
  double res=0.0;
  
  theAmplitude->matrix_tree_pmpm(c[0], c[1], c[3], c[4], tmp);
  res+=norm(tmp[0]);
  theAmplitude->matrix_tree_pmmp(c[0], c[1], c[3], c[4], tmp);
  res+=norm(tmp[0]);

  return 2.*res;
}

void NLOJetAmplitudeq4::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {
  if ( !theAmplitude ) {
    static NLOJetRandomWrapper dummy;
    theAmplitude = new nlo::ampq4(lastInvariants(),dummy);
  }
  NLOJetAmplitude<0,2,0>::prepareAmplitudes(me);
}

double NLOJetAmplitudeq4::treeLevel2(const vector<int>& c) const {
  static double res[2];
  theAmplitude->su3_tree(c[0], c[1], c[2], c[3], res);
  return identicalQuarks() == 2 ? res[1] : res[0];
}

double NLOJetAmplitudeq4::treeOneLoop(const vector<int>& c) const {
  static double res[2];
  theAmplitude->su3_1loop(nLight(), c[0], c[1], c[2], c[3], res);
  return identicalQuarks() == 2 ? res[1] : res[0];
}

double NLOJetAmplitudeq4::treeLevelCC(pair<int,int> ij,
				      const vector<int>& c) const {
  static double res[2];
  theAmplitude->su3_cc(ij.first,ij.second,c[0], c[1], c[2], c[3], res);
  return identicalQuarks() == 2 ? res[1] : res[0];
}

pair<double,Complex> NLOJetAmplitudeq4::treeLevelSCC(pair<int,int>,
						     const vector<int>&) const {
  throw Exception() << name() << " does not provide spin-colour-correlated matrix elements"
		    << Exception::abortnow;
  return make_pair(0.,0.);
}

IBPtr NLOJetAmplitudeq4::clone() const {
  return new_ptr(*this);
}

IBPtr NLOJetAmplitudeq4::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void NLOJetAmplitudeq4::persistentOutput(PersistentOStream &) const {
}

void NLOJetAmplitudeq4::persistentInput(PersistentIStream &, int) {
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<NLOJetAmplitudeq4,Herwig::NLOJetAmplitude<0,2,0> >
  describeHerwigNLOJetAmplitudeq4("Herwig::NLOJetAmplitudeq4", "HwMatchboxNLOJet.so");

void NLOJetAmplitudeq4::Init() {

  static ClassDocumentation<NLOJetAmplitudeq4> documentation
    ("NLOJetAmplitudeq4");

}

