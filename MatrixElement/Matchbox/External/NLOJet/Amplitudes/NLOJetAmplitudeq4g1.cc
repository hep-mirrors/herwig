// -*- C++ -*-
//
// NLOJetAmplitudeq4g1.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetAmplitudeq4g1 class.
//

#include "NLOJetAmplitudeq4g1.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/NLOJetRandomWrapper.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEqbq2kkbg.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEqqb2kkbg.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEkbq2qkbg.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEkqb2kqbg.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEkbqb2kbqbg.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEkq2kqg.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEgq2kqkb.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEgqb2kkbqb.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEqg2kqkb.h"
#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/PP2Jets/NLOJetMEqbg2kkbqb.h"


using namespace Herwig;

NLOJetAmplitudeq4g1::NLOJetAmplitudeq4g1() 
  : NLOJetAmplitude<0,2,0>(), theAmplitude(0) {}

NLOJetAmplitudeq4g1::~NLOJetAmplitudeq4g1() {
  if ( theAmplitude ) {
    delete theAmplitude;
    theAmplitude = 0;
  }
}

bool NLOJetAmplitudeq4g1::canHandle(const PDVector& pd) const {

  pair<unsigned int, unsigned int> n = countColoured(pd);
  return n.first == 4 && n.second == 1;

}

Ptr<MatchboxMEBase>::ptr NLOJetAmplitudeq4g1::makeME(const vector<PDVector>& pd) const {

  if ( pd.front()[0]->id() == 21 ) {
    if ( pd.front()[1]->id() > 0 )
      return new_ptr(NLOJetMEgq2kqkb());
    if ( pd.front()[1]->id() < 0 )
      return new_ptr(NLOJetMEgqb2kkbqb());
  }

  if ( pd.front()[1]->id() == 21 ) {
    if ( pd.front()[0]->id() > 0 )
      return new_ptr(NLOJetMEqg2kqkb());
    if ( pd.front()[0]->id() < 0 )
      return new_ptr(NLOJetMEqbg2kkbqb());
  }

  if ( pd.front()[0]->id() + pd.front()[1]->id() == 0 ) {
    if ( pd.front()[0]->id() < 0 )
      return new_ptr(NLOJetMEqbq2kkbg());
    if ( pd.front()[0]->id() > 0 )
      return new_ptr(NLOJetMEqqb2kkbg());
  }

  if ( pd.front()[0]->id() < 0 &&
       pd.front()[1]->id() > 0 ) {
    return new_ptr(NLOJetMEkbq2qkbg());
  }

  if ( pd.front()[0]->id() > 0 &&
       pd.front()[1]->id() < 0 ) {
    return new_ptr(NLOJetMEkqb2kqbg());
  }

  if ( pd.front()[0]->id() < 0 &&
       pd.front()[1]->id() < 0 ) {
    return new_ptr(NLOJetMEkbqb2kbqbg());
  }

  if ( pd.front()[0]->id() > 0 &&
       pd.front()[1]->id() > 0 ) {
    return new_ptr(NLOJetMEkq2kqg());
  }

  assert(false);
  return new_ptr(MatchboxMEBase());

}

double NLOJetAmplitudeq4g1::colourOrdered2(const int* c, size_t) const {
  static const int LEND= -999;

  nlo::ampq4g1::amp_tree tmp[2];
  double res=0.0;
 
  if(c[2] == LEND){
    theAmplitude->matrix_tree_pmpmp(c[0], c[1], c[3], c[5], c[4], tmp);
    res+=norm(tmp[0].A12);
    theAmplitude->matrix_tree_pmmpp(c[0], c[1], c[3], c[5], c[4], tmp);
    res+=norm(tmp[0].A12);
    theAmplitude->matrix_tree_mppmp(c[0], c[1], c[3], c[5], c[4], tmp);
    res+=norm(tmp[0].A12);
    theAmplitude->matrix_tree_mpmpp(c[0], c[1], c[3], c[5], c[4], tmp);
    res+=norm(tmp[0].A12);
  }
  else if(c[3] == LEND){
    theAmplitude->matrix_tree_pmpmp(c[0], c[2], c[4], c[5], c[1], tmp);
    res+=norm(tmp[0].A34);
    theAmplitude->matrix_tree_pmmpp(c[0], c[2], c[4], c[5], c[1], tmp);
    res+=norm(tmp[0].A34);
    theAmplitude->matrix_tree_mppmp(c[0], c[2], c[4], c[5], c[1], tmp);
    res+=norm(tmp[0].A34);
    theAmplitude->matrix_tree_mpmpp(c[0], c[2], c[4], c[5], c[1], tmp);
    res+=norm(tmp[0].A34);
  }
  else{
    assert(false && "not proper values of colourOrdered2 argument");
  }

  return 2.*res;
}

void NLOJetAmplitudeq4g1::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {
  if ( !theAmplitude ) {
    static NLOJetRandomWrapper dummy;
    theAmplitude = new nlo::ampq4g1(lastInvariants(),dummy);
  }
  NLOJetAmplitude<0,2,0>::prepareAmplitudes(me);
}

double NLOJetAmplitudeq4g1::treeLevel2(const vector<int>& c) const {
  static double res[2];
  theAmplitude->su3_tree(c[0], c[1], c[2], c[3], c[4], res);
  return identicalQuarks() == 2 ? res[1] : res[0];
}

double NLOJetAmplitudeq4g1::treeOneLoop(const vector<int>& c) const {
  static double res[2];
  theAmplitude->su3_1loop(nLight(), c[0], c[1], c[2], c[3], c[4], res);
  return identicalQuarks() == 2 ? res[1] : res[0];
}

double NLOJetAmplitudeq4g1::treeLevelCC(pair<int,int> ij,
					const vector<int>& c) const {
  static pair<double,complex<double> > res[2];
  theAmplitude->su3_cc(ij.first,ij.second,c[0], c[1], c[2], c[3], c[4], res);
  return identicalQuarks() == 2 ? res[1].first : res[0].first;
}

pair<double,Complex> NLOJetAmplitudeq4g1::treeLevelSCC(pair<int,int> ij,
						       const vector<int>& c) const {
  static pair<double,complex<double> > res[2];
  theAmplitude->su3_cc(ij.first,ij.second,c[0], c[1], c[2], c[3], c[4], res);
  return identicalQuarks() == 2 ? res[1] : res[0];
}

IBPtr NLOJetAmplitudeq4g1::clone() const {
  return new_ptr(*this);
}

IBPtr NLOJetAmplitudeq4g1::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void NLOJetAmplitudeq4g1::persistentOutput(PersistentOStream &) const {
}

void NLOJetAmplitudeq4g1::persistentInput(PersistentIStream &, int) {
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<NLOJetAmplitudeq4g1,Herwig::NLOJetAmplitude<0,2,0> >
  describeHerwigNLOJetAmplitudeq4g1("Herwig::NLOJetAmplitudeq4g1", "HwMatchboxNLOJet.so");

void NLOJetAmplitudeq4g1::Init() {

  static ClassDocumentation<NLOJetAmplitudeq4g1> documentation
    ("NLOJetAmplitudeq4g1");

}

