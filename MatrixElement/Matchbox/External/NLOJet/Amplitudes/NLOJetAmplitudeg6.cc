// -*- C++ -*-
//
// NLOJetAmplitudeg6.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetAmplitudeg6 class.
//

#include "NLOJetAmplitudeg6.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/NLOJetRandomWrapper.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEgg2gggg.h"

using namespace Herwig;

NLOJetAmplitudeg6::NLOJetAmplitudeg6() 
  : NLOJetAmplitude<0,2,0>(), theAmplitude(0) {}

NLOJetAmplitudeg6::~NLOJetAmplitudeg6() {
  if ( theAmplitude ) {
    delete theAmplitude;
    theAmplitude = 0;
  }
}

bool NLOJetAmplitudeg6::canHandle(const PDVector& pd) const {

  pair<unsigned int, unsigned int> n = countColoured(pd);
  return n.first == 0 && n.second == 6;

}

Ptr<MatchboxMEBase>::ptr NLOJetAmplitudeg6::makeME(const vector<PDVector>&) const {
  return new_ptr(NLOJetMEgg2gggg());
}

double NLOJetAmplitudeg6::colourOrdered2(const int* c, size_t) const {
  complex<double> mhv;
  complex<double> nmhv[11];
  double res=0.0;
  static const int pr[]={0,1, 0,2, 0,3, 0,4, 0,5, 1,2, 1,3, 1,4, 1,5, 2,3, 2,4, 2,5, 3,4, 3,5, 4,5   };
  static int nopr=sizeof(pr)/sizeof(int)/2 ;

  for(int i=0;i< nopr ;i++){
    theAmplitude->matrix_tree_MHV(c[pr[2*i]], c[pr[2*i+1]], c[0], c[1], c[2], c[3], c[4], c[5], &mhv);
    res+=norm(mhv);
  }
  theAmplitude->matrix_tree_NMHV(c[0], c[1], c[2], c[3], c[4], c[5], nmhv);
  for(int i=1;i< 11 ;i++){
    res+=norm(nmhv[i]);
  }

  return 2.*res;
}

void NLOJetAmplitudeg6::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {
  if ( !theAmplitude ) {
    static NLOJetRandomWrapper dummy;
    theAmplitude = new nlo::ampg6(lastInvariants(),dummy);
  }
  NLOJetAmplitude<0,2,0>::prepareAmplitudes(me);
}

double NLOJetAmplitudeg6::treeLevel2(const vector<int>& c) const {
  return theAmplitude->su3_tree(c[0], c[1], c[2], c[3], c[4], c[5]);
}

double NLOJetAmplitudeg6::treeOneLoop(const vector<int>&) const {
  throw Exception() << name() << " does not provide one-loop matrix elements"
		    << Exception::abortnow;
  return 0.;
}

double NLOJetAmplitudeg6::treeLevelCC(pair<int,int>,
				      const vector<int>&) const {
  throw Exception() << name() << " does not provide colour-correlated matrix elements"
		    << Exception::abortnow;
  return 0.;
}

pair<double,Complex> NLOJetAmplitudeg6::treeLevelSCC(pair<int,int>,
						      const vector<int>&) const {
  throw Exception() << name() << " does not provide colour-correlated matrix elements"
		    << Exception::abortnow;
  return make_pair(0.,0.);
}

IBPtr NLOJetAmplitudeg6::clone() const {
  return new_ptr(*this);
}

IBPtr NLOJetAmplitudeg6::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void NLOJetAmplitudeg6::persistentOutput(PersistentOStream &) const {}

void NLOJetAmplitudeg6::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<NLOJetAmplitudeg6,Herwig::NLOJetAmplitude<0,2,0> >
  describeHerwigNLOJetAmplitudeg6("Herwig::NLOJetAmplitudeg6", "HwMatchboxNLOJet.so");

void NLOJetAmplitudeg6::Init() {

  static ClassDocumentation<NLOJetAmplitudeg6> documentation
    ("NLOJetAmplitudeg6");

}

