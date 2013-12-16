// -*- C++ -*-
//
// NLOJetAmplitudeq6.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetAmplitudeq6 class.
//

#include "NLOJetAmplitudeq6.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/NLOJetRandomWrapper.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEqqb2rkrbkb.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEqbq2rkrbkb.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEkq2rkqrb.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEkqb2rkrbqb.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEkbq2rqrbkb.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEkbqb2rrbkbqb.h"

using namespace Herwig;

NLOJetAmplitudeq6::NLOJetAmplitudeq6() 
  : NLOJetAmplitude<0,2,0>(), theAmplitude(0) {}

NLOJetAmplitudeq6::~NLOJetAmplitudeq6() {
  if ( theAmplitude ) {
    delete theAmplitude;
    theAmplitude = 0;
  }
}

bool NLOJetAmplitudeq6::canHandle(const PDVector& pd) const {

  pair<unsigned int, unsigned int> n = countColoured(pd);
  return n.first == 6 && n.second == 0;

}

Ptr<MatchboxMEBase>::ptr NLOJetAmplitudeq6::makeME(const vector<PDVector>& pd) const {

  if ( pd.front()[0]->id() + pd.front()[1]->id() == 0 ) {
    if ( pd.front()[0]->id() > 0 )
      return new_ptr(NLOJetMEqqb2rkrbkb());
    if ( pd.front()[0]->id() < 0 )
      return new_ptr(NLOJetMEqbq2rkrbkb());
  }

  if ( pd.front()[0]->id() > 0 ) {
    if ( pd.front()[1]->id() > 0 )
      return new_ptr(NLOJetMEkq2rkqrb());
    if ( pd.front()[1]->id() < 0 )
      return new_ptr(NLOJetMEkqb2rkrbqb());
  }

  if ( pd.front()[0]->id() < 0 ) {
    if ( pd.front()[1]->id() > 0 )
      return new_ptr(NLOJetMEkbq2rqrbkb());
    if ( pd.front()[1]->id() < 0 )
      return new_ptr(NLOJetMEkbqb2rrbkbqb());
  }

  assert(false);
  return new_ptr(MatchboxMEBase());

}

double NLOJetAmplitudeq6::colourOrdered2(const int* c, size_t) const {
  std::complex<double> tmp[4];
  double res=0.;

  theAmplitude->matrix_tree(c[0], c[1], c[3], c[4], c[6], c[7], tmp);

  for(int i=0;i<4;i++)
    res+=norm(tmp[i]);
  
  return 2.*res;
}

void NLOJetAmplitudeq6::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {
  if ( !theAmplitude ) {
    static NLOJetRandomWrapper dummy;
    theAmplitude = new nlo::ampq6(lastInvariants(),dummy);
  }
  NLOJetAmplitude<0,2,0>::prepareAmplitudes(me);
}

double NLOJetAmplitudeq6::treeLevel2(const vector<int>& c) const {

  // work out the magic strings
  string what = "";
  size_t id = 0;

  if ( amplitudePartonData()[0]->id() != amplitudePartonData()[1]->id() &&
       amplitudePartonData()[1]->id() != amplitudePartonData()[2]->id() &&
       amplitudePartonData()[0]->id() != amplitudePartonData()[2]->id() ) {
    id = 0; what = "10000";
  } else if ( amplitudePartonData()[0]->id() == amplitudePartonData()[1]->id() &&
	      amplitudePartonData()[1]->id() != amplitudePartonData()[2]->id() &&
	      amplitudePartonData()[0]->id() != amplitudePartonData()[2]->id() ) {
    id = 1; what = "01000";
  } else if ( amplitudePartonData()[0]->id() != amplitudePartonData()[1]->id() &&
	      amplitudePartonData()[1]->id() != amplitudePartonData()[2]->id() &&
	      amplitudePartonData()[0]->id() == amplitudePartonData()[2]->id() ) {
    id = 2; what = "00100";
  } else if ( amplitudePartonData()[0]->id() != amplitudePartonData()[1]->id() &&
	      amplitudePartonData()[1]->id() == amplitudePartonData()[2]->id() &&
	      amplitudePartonData()[0]->id() != amplitudePartonData()[2]->id() ) {
    id = 3; what = "00010";
  } else if ( amplitudePartonData()[0]->id() == amplitudePartonData()[1]->id() &&
	      amplitudePartonData()[1]->id() == amplitudePartonData()[2]->id() &&
	      amplitudePartonData()[0]->id() == amplitudePartonData()[2]->id() ) {
    id = 4; what = "00001";
  }

  assert(what != "");

  static double res[5];

  theAmplitude->su3_tree(c[0],c[1],c[2],c[3],c[4],c[5],what.c_str(),res);

  return res[id];

}

double NLOJetAmplitudeq6::treeOneLoop(const vector<int>&) const {
  throw Exception() << name() << " does not provide one-loop matrix elements"
		    << Exception::abortnow;
  return 0.;
}

double NLOJetAmplitudeq6::treeLevelCC(pair<int,int>,
				      const vector<int>&) const {
  throw Exception() << name() << " does not provide colour-correlated matrix elements"
		    << Exception::abortnow;
  return 0.;
}

pair<double,Complex> NLOJetAmplitudeq6::treeLevelSCC(pair<int,int>,
						     const vector<int>&) const {
  throw Exception() << name() << " does not provide colour-correlated matrix elements"
		    << Exception::abortnow;
  return make_pair(0.,0.);
}

IBPtr NLOJetAmplitudeq6::clone() const {
  return new_ptr(*this);
}

IBPtr NLOJetAmplitudeq6::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void NLOJetAmplitudeq6::persistentOutput(PersistentOStream &) const {}

void NLOJetAmplitudeq6::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<NLOJetAmplitudeq6,Herwig::NLOJetAmplitude<0,2,0> >
  describeHerwigNLOJetAmplitudeq6("Herwig::NLOJetAmplitudeq6", "HwMatchboxNLOJet.so");

void NLOJetAmplitudeq6::Init() {

  static ClassDocumentation<NLOJetAmplitudeq6> documentation
    ("NLOJetAmplitudeq6");

}

