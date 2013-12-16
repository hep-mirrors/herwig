// -*- C++ -*-
//
// NLOJetAmplitudeg4.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetAmplitudeg4 class.
//

#include "NLOJetAmplitudeg4.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/NLOJetRandomWrapper.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEgg2gg.h"

using namespace Herwig;

NLOJetAmplitudeg4::NLOJetAmplitudeg4() 
  : NLOJetAmplitude<0,2,0>(), theAmplitude(0) {}

NLOJetAmplitudeg4::~NLOJetAmplitudeg4() {
  if ( theAmplitude ) {
    delete theAmplitude;
    theAmplitude = 0;
  }
}

bool NLOJetAmplitudeg4::canHandle(const PDVector& pd) const {

  pair<unsigned int, unsigned int> n = countColoured(pd);
  return n.first == 0 && n.second == 4;

}

Ptr<MatchboxMEBase>::ptr NLOJetAmplitudeg4::makeME(const vector<PDVector>&) const {
  return new_ptr(NLOJetMEgg2gg());
}

double NLOJetAmplitudeg4::colourOrdered2(const int* c, size_t) const {
  complex<double> tmp;
  double res=0.0;
  static const int pr[]={0,1, 0,2, 0,3, 1,2, 1,3, 2,3   };
  static int nopr=sizeof(pr)/sizeof(int)/2 ;

  for(int i=0;i< nopr ;i++){
    theAmplitude->matrix_tree(c[pr[2*i]], c[pr[2*i+1]], c[0], c[1], c[2], c[3], &tmp);
    res+=norm(tmp);
  }

  return res;
}

void NLOJetAmplitudeg4::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {
  if ( !theAmplitude ) {
    static NLOJetRandomWrapper dummy;
    theAmplitude = new nlo::ampg4(lastInvariants(),dummy);
  }
  NLOJetAmplitude<0,2,0>::prepareAmplitudes(me);
}

double NLOJetAmplitudeg4::treeLevel2(const vector<int>& c) const {
  return theAmplitude->su3_tree(c[0], c[1], c[2], c[3]);
}

double NLOJetAmplitudeg4::treeOneLoop(const vector<int>& c) const {
  return theAmplitude->su3_1loop(nLight(),c[0], c[1], c[2], c[3]);
}

double NLOJetAmplitudeg4::treeLevelCC(pair<int,int> ij,
				      const vector<int>& c) const {
  return theAmplitude->su3_cc(ij.first,ij.second,c[0], c[1], c[2], c[3]);
}

pair<double,Complex> NLOJetAmplitudeg4::treeLevelSCC(pair<int,int> ij,
						     const vector<int>& c) const {
  return make_pair(theAmplitude->su3_cc(ij.first,ij.second,c[0], c[1], c[2], c[3]),0.);
}

IBPtr NLOJetAmplitudeg4::clone() const {
  return new_ptr(*this);
}

IBPtr NLOJetAmplitudeg4::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void NLOJetAmplitudeg4::persistentOutput(PersistentOStream &) const {}

void NLOJetAmplitudeg4::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<NLOJetAmplitudeg4,Herwig::NLOJetAmplitude<0,2,0> >
  describeHerwigNLOJetAmplitudeg4("Herwig::NLOJetAmplitudeg4", "HwMatchboxNLOJet.so");

void NLOJetAmplitudeg4::Init() {

  static ClassDocumentation<NLOJetAmplitudeg4> documentation
    ("NLOJetAmplitudeg4");

}

