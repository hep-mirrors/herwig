// -*- C++ -*-
//
// NLOJetAmplitudeg5.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetAmplitudeg5 class.
//

#include "NLOJetAmplitudeg5.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/NLOJetRandomWrapper.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEgg2ggg.h"

using namespace Herwig;

NLOJetAmplitudeg5::NLOJetAmplitudeg5() 
  : NLOJetAmplitude<0,2,0>(), theAmplitude(0) {}

NLOJetAmplitudeg5::~NLOJetAmplitudeg5() {
  if ( theAmplitude ) {
    delete theAmplitude;
    theAmplitude = 0;
  }
}

bool NLOJetAmplitudeg5::canHandle(const PDVector& pd) const {

  pair<unsigned int, unsigned int> n = countColoured(pd);
  return n.first == 0 && n.second == 5;

}

Ptr<MatchboxMEBase>::ptr NLOJetAmplitudeg5::makeME(const vector<PDVector>&) const {
  return new_ptr(NLOJetMEgg2ggg());
}

double NLOJetAmplitudeg5::colourOrdered2(const int* c, size_t) const {
  complex<double> tmp;
  double res=0.0;
  static const int pr[]={0,1, 0,2, 0,3, 0,4, 1,2, 1,3, 1,4, 2,3, 2,4, 3,4  };
  static int nopr=sizeof(pr)/sizeof(int)/2 ;

  for(int i=0;i< nopr ;i++){
    theAmplitude->matrix_tree(c[pr[2*i]], c[pr[2*i+1]], c[0], c[1], c[2], c[3], c[4], &tmp);
    res+=norm(tmp);
  }

  return 2.*res;
}

void NLOJetAmplitudeg5::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {
  if ( !theAmplitude ) {
    static NLOJetRandomWrapper dummy;
    theAmplitude = new nlo::ampg5(lastInvariants(),dummy);
  }
  NLOJetAmplitude<0,2,0>::prepareAmplitudes(me);
}

double NLOJetAmplitudeg5::treeLevel2(const vector<int>& c) const {
  return theAmplitude->su3_tree(c[0], c[1], c[2], c[3], c[4]);
}

double NLOJetAmplitudeg5::treeOneLoop(const vector<int>& c) const {
  return theAmplitude->su3_1loop(nLight(),c[0], c[1], c[2], c[3], c[4]);
}

double NLOJetAmplitudeg5::treeLevelCC(pair<int,int> ij,
				      const vector<int>& c) const {
  return theAmplitude->su3_cc(ij.first,ij.second,c[0], c[1], c[2], c[3], c[4]).first;
}

pair<double,Complex> NLOJetAmplitudeg5::treeLevelSCC(pair<int,int> ij,
						     const vector<int>& c) const {
  return theAmplitude->su3_cc(ij.first,ij.second,c[0], c[1], c[2], c[3], c[4]);
}

IBPtr NLOJetAmplitudeg5::clone() const {
  return new_ptr(*this);
}

IBPtr NLOJetAmplitudeg5::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void NLOJetAmplitudeg5::persistentOutput(PersistentOStream &) const {}

void NLOJetAmplitudeg5::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<NLOJetAmplitudeg5,Herwig::NLOJetAmplitude<0,2,0> >
  describeHerwigNLOJetAmplitudeg5("Herwig::NLOJetAmplitudeg5", "HwMatchboxNLOJet.so");

void NLOJetAmplitudeg5::Init() {

  static ClassDocumentation<NLOJetAmplitudeg5> documentation
    ("NLOJetAmplitudeg5");

}

